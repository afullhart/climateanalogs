
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess as sub
import os
import pandas as pd
from classes.Formatting import Formatting
formatting_obj = Formatting()
import shutil
from subprocess import Popen, PIPE
from osgeo import osr
from osgeo import gdal
from scipy import optimize
import numpy as np
from datetime import datetime

cliDIR = '/home/afullhart/Downloads/cligen_53004_Linux'
parentparDIR = '/home/afullhart/Downloads/GDBs/RCP45/pars'
ptFILE = '/home/afullhart/Downloads/Points.csv'
outDIR = '/home/afullhart/Downloads/GDBs/RCP45/maps'


parFolders = [os.path.join(parentparDIR, d) for d in os.listdir(parentparDIR)]

n_workers = 100
REC_LEN = 30
eo = 0.29; a = 0.72; io = 12.195

ndays_ref = [31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

trans = (-120.0, 0.0083333333, 0.0, 31.3333333333, 0.0, 0.0083333333)
spatial_ref = osr.SpatialReference()
spatial_ref.ImportFromEPSG(4326)

driver = gdal.GetDriverByName('GTiff')
data_type = gdal.GDT_Float32

with open(ptFILE) as f:
  point_list = []
  next(f)
  for line in f:
    row = line.split(',')
    i = int(row[0])
    x = float(row[1])
    y = float(row[2])
    point_list.append([i, x, y])

def coord2pixel(x, y, trans):
  xpt = (((x - trans[0]) / trans[1]))
  ypt = (((y - trans[3]) / trans[5]))
  return(xpt, ypt)

def newton(b_, ip_):
  return (ip_ * (1 - np.exp(-b_))) - (b_)

def i30(p_, ip_, d_, b_):
  return ((2*p_*ip_)/(b_)) * (1 - np.exp(-((b_)/(2*d_))))

def energy(p_, ip_, lp_, b_, eo_, a_, io_,):
  inside = np.exp( -(lp_ / io_)*np.exp( -b_ ) ) - np.exp( -lp_ / io_ )
  middle = ((a_*ip_)/(b_)) * (io_/lp_)
  outside = p_*eo_
  return (outside)*((1)-(middle*inside))

def make_dtime(yr, mo, dy):
  dtime = datetime(yr.astype(int), mo.astype(int), dy.astype(int))
  return dtime

def daylength(doy, lat):
  latInRad = np.deg2rad(lat)
  dec = 23.45*np.sin(np.deg2rad(360.0*(283.0+doy)/365.0))
  if -np.tan(latInRad) * np.tan(np.deg2rad(dec)) <= -1.0:
    return 24.0
  elif -np.tan(latInRad) * np.tan(np.deg2rad(dec)) >= 1.0:
    return 0.0
  else:
    hourAngle = np.rad2deg(np.arccos(-np.tan(latInRad) * np.tan(np.deg2rad(dec))))
    return 2.0*hourAngle/15.0


def thornthwaite(lines, lat):
  
  data = []
  for line in lines:
    line = line.replace('*****', ' 0.0 ')
    row = [float(x) for x in line.strip('\n').split()]
    data.append([row[0], row[1], row[2]+2000, row[3], row[7], row[8]])
  
  df = pd.DataFrame(data, columns=['day', 'mo', 'yr', 'precip', 'tmax', 'tmin'])
  df['dtime'] = df.apply(lambda row: make_dtime(row['yr'], row['mo'], row['day']), axis=1)
  df['doy'] = df.apply(lambda row: row['dtime'].timetuple().tm_yday, axis=1)
  df['dayhrs'] = df.apply(lambda row: daylength(row['doy'], lat), axis=1)
  dyhrs_monthly = df.groupby('mo').mean()['dayhrs']
  df['tmean'] = (df['tmax'] + df['tmin']) / 2
  tmean_monthly = df.groupby('mo').mean()['tmean']
  
  I = sum([(x/5)**1.514 for x in tmean_monthly])
  a = (6.75e-7)*I**3 - (7.71e-5)*I**2 + (1.792e-2)*I + 0.49239
  pet_monthly = []
  for i in range(12):
    L = dyhrs_monthly[i+1]
    N = ndays_ref[i]
    Td = tmean_monthly[i+1]
    pet_mo = 16*(L/12.0)*(N/30.0)*(10*Td/I)**a 
    pet_monthly.append(pet_mo)
  
  pet = sum(pet_monthly)

  return pet


#cd /home/afullhart/Downloads/cligen_53004_Linux && parallel --progress "script -q -c './cligen_53004_Linux -b1 -y30 -t5 -i{} -o{.}.txt' /dev/null > /dev/null 2>&1" ::: *.par


def run_cligen(lbl):
  output_file = os.path.join(cliDIR, f'{lbl}.txt')
  command = f"""script -q -c 'cd {cliDIR} && ./cligen_53004_Linux -b1 -y{REC_LEN} -t5 -i{lbl}.par -o{lbl}.txt' /dev/null > /dev/null 2>&1"""
  status = os.system(command)
  if status == 0 and os.path.exists(output_file):
    good = 1
  else:
    good = 0

def make_Zgrid(xyz_df):
  x_data = np.array(xyz_df['x'])
  y_data = np.array(xyz_df['y'])
  z_data = np.array(xyz_df['z'])
  unique_x = np.unique(x_data)
  unique_y = np.unique(y_data)
  X_grid, Y_grid = np.meshgrid(unique_x, unique_y)
  Z_grid = np.zeros(X_grid.shape)
  Z_grid = np.empty(Z_grid.shape)
  Z_grid.fill(np.nan)
  for i in range(len(x_data)):
    x_val = x_data[i]
    y_val = y_data[i]
    z_val = z_data[i]
    x_idx = np.where(unique_x == x_val)[0][0]
    y_idx = np.where(unique_y == y_val)[0][0]
    Z_grid[y_idx, x_idx] = z_val
  return Z_grid

def make_gtif(Z_grid, tif_f):
  n_rows, n_cols = Z_grid.shape
  mem_raster = driver.Create(tif_f, n_cols, n_rows, 1, data_type)
  mem_raster.SetGeoTransform(trans)
  mem_raster.SetProjection(spatial_ref.ExportToWkt())
  mem_raster.GetRasterBand(1).WriteArray(Z_grid)
  mem_path = '/vsimem/mem_raster.vrt'
  gdal.Warp(mem_path, mem_raster)
  translate_options = gdal.TranslateOptions(
    format='GTiff', 
    outputSRS='EPSG:4326',
    noData=np.nan
  )
  gdal.Translate(tif_f, mem_path, options=translate_options)
  


def main(point, parDIR):

  fname = str(point[0])
  shutil.copyfile(os.path.join(parDIR, fname + '.par'), os.path.join(cliDIR, fname + '.par'))
  
  run_cligen(fname)

  with open(os.path.join(cliDIR, fname + '.txt')) as f:
    lines = f.readlines()[15:-1]
  
  os.remove(os.path.join(cliDIR, fname + '.txt'))
  os.remove(os.path.join(cliDIR, fname + '.par'))

  ei_sum = 0.0
  swe_sum = 0.0
  accum = 0.0
  ndays_ct = 0
  consecutive_ct = 0
  yr_last = 1
  consecutive_list = []
  consecutive_sublist = []
  for line in lines:

    row = line.split()
    yr = int(row[2])
    mo = int(row[1])

    if yr != yr_last:
      consecutive_sublist.append(consecutive_ct)
      consecutive_list.append(max(consecutive_sublist))
      consecutive_sublist = []
      consecutive_ct = 0

    if row[3] != '*****':
      p = float(row[3])
      tavg = (float(row[7]) + float(row[8])) / 2.

      if p >= 12.7 and tavg > 5.0:
        ip = float(row[6])
        b = optimize.newton(func=newton, x0=ip, args=(ip,))
        d = float(row[4])
        lp = ip * (p/d)
        if d > 0.5:
            l30 = i30(p, ip, d, b)
        else:
            l30 = 2*p
        e = energy(p, ip, lp, b, eo, a, io)
        ei = e*l30
        ei_sum += ei

      else:
        ei = 0.
        ei_sum += ei

      if p > 0.0 and tavg <= 5.0:
        swe_sum += p

      if p < 0.3:
        consecutive_ct += 1
      else:
        consecutive_ct = 0
        ndays_ct += 1

      yr_last = yr

  consecutive_sublist.append(consecutive_ct)
  consecutive_list.append(max(consecutive_sublist))

  ann_ero = ei_sum / REC_LEN
  ann_swe = swe_sum / REC_LEN
  ann_ndays = float(ndays_ct) / REC_LEN
  ann_consec = np.mean(consecutive_list)
  ann_pet = thornthwaite(lines, point[2])

  return([point[1], point[2], ann_ero, ann_swe, ann_ndays, ann_consec, ann_pet])



if __name__ == '__main__':

  for pdir in parFolders:

    xyz_results = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
      pool = {executor.submit(main, p, pdir): p for p in point_list}
      ct = 0
      for future in as_completed(pool):
        res = future.result()
        xyz_results.append([res[0], res[1], res[2], res[3], res[4], res[5], res[6]])
        ct += 1
        if ct % 100000 == 0:
          print(ct)

    results_df = pd.DataFrame(xyz_results, columns=['x','y','z1','z2', 'z3', 'z4', 'z5'])
    results_df.sort_values(by=['y', 'x'], ascending=[True, True], axis=0, inplace=True)

    ero_df = results_df[['x', 'y', 'z1']]
    ero_df = ero_df.rename(columns={'x':'x', 'y':'y', 'z1':'z'})
    Z_grid_ero = make_Zgrid(ero_df)

    swe_df = results_df[['x', 'y', 'z2']]
    swe_df = swe_df.rename(columns={'x':'x', 'y':'y', 'z2':'z'})
    Z_grid_swe = make_Zgrid(swe_df)

    ndays_df = results_df[['x', 'y', 'z3']]
    ndays_df = ndays_df.rename(columns={'x':'x', 'y':'y', 'z3':'z'})
    Z_grid_ndays = make_Zgrid(ndays_df)

    consec_df = results_df[['x', 'y', 'z4']]
    consec_df = consec_df.rename(columns={'x':'x', 'y':'y', 'z4':'z'})
    Z_grid_consec = make_Zgrid(consec_df)

    pet_df = results_df[['x', 'y', 'z5']]
    pet_df = pet_df.rename(columns={'x':'x', 'y':'y', 'z5':'z'})
    Z_grid_pet = make_Zgrid(pet_df)

    scenario = os.path.split(pdir)[-1]

    ero_tif_f = os.path.join(outDIR, scenario + '_ero.tif')
    swe_tif_f = os.path.join(outDIR, scenario + '_swe.tif')
    ndays_tif_f = os.path.join(outDIR, scenario + '_ndays.tif')
    consec_tif_f = os.path.join(outDIR, scenario + '_consec.tif')

    make_gtif(Z_grid_ero, ero_tif_f)
    make_gtif(Z_grid_swe, swe_tif_f)
    make_gtif(Z_grid_ndays, ndays_tif_f)
    make_gtif(Z_grid_consec, consec_tif_f)

print('ALL DONE')


