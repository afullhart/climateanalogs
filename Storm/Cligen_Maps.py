
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


cliDIR = '/home/afullhart/Downloads/cligen_53004_Linux'
parDIR = '/home/afullhart/Downloads/GDBs/RCP45/pars/CCSM4_2000_2029'
ptFILE = '/home/afullhart/Downloads/Points.csv'
outDIR = '/home/afullhart/Downloads/GDBs/RCP45/maps'


var_labels = ['mean', 'sdev', 'skew', 'pww', 'pwd', 'tmax', 'tmin', 'txsd', 'tnsd', 'srad', 'srsd', 'mx5p', 'tdew', 'timepk']
historical_var_labels = ['timepk']

yr_str = '1974_2013'
n_workers = 100
REC_LEN = 30
eo = 0.29; a = 0.72; io = 12.195

trans = (-120.0, 0.0083333333, 0.0, 31.3333333333, 0.0, 0.0083333333)
spatial_ref = osr.SpatialReference()
spatial_ref.ImportFromWkt('WGS84')
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
  Z_grid = np.flipud(Z_grid)
  return Z_grid

def make_gtif(Z_grid, tif_f):
  n_rows, n_cols = Z_grid.shape
  mem_raster = driver.Create(tif_f, n_cols, n_rows, 1, data_type)
  mem_raster.SetGeoTransform(trans)
  mem_raster.SetProjection(spatial_ref.ExportToWkt())
  mem_raster.GetRasterBand(1).WriteArray(Z_grid)
  #will this collide if multiprocessing?
  mem_path = '/vsimem/mem_raster.vrt'
  gdal.Warp(mem_path, mem_raster)
  translate_options = gdal.TranslateOptions(
    format='GTiff', 
    outputSRS='EPSG:4326', 
    outputGeotransform=trans,
    noData=np.nan
  )
  gdal.Translate(tif_f, mem_path, options=translate_options)
  


def main(point):

  fname = str(point[0])
  shutil.copyfile(os.path.join(parDIR, fname + '.par'), os.path.join(cliDIR, fname + '.par'))
  
  run_cligen(fname)

  with open(os.path.join(cliDIR, fname + '.txt')) as f:
    lines = f.readlines()[15:]
  
  os.remove(os.path.join(cliDIR, fname + '.txt'))
  os.remove(os.path.join(cliDIR, fname + '.par'))

  ei_sum = 0.0
  swe_sum = 0.0
  accum = 0.0
  for line in lines[:-1]:
    row = line.split()
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

  ann_ero = ei_sum / REC_LEN
  ann_swe = swe_sum / REC_LEN

  return([point[1], point[2], ann_ero, ann_swe])



if __name__ == '__main__':
  xyz_results = []
  with ProcessPoolExecutor(max_workers=n_workers) as executor:
    pool = {executor.submit(main, p): p for p in point_list}
    ct = 0
    for future in as_completed(pool):
      res = future.result()
      xyz_results.append([res[0], res[1], res[2], res[3]])
      ct += 1
      if ct % 100000 == 0:
        print(ct)

  results_df = pd.DataFrame(xyz_results, columns=['x','y','z1','z2'])
  results_df.sort_values(by=['y', 'x'], ascending=[True, True], axis=0, inplace=True)

  ero_df = results_df[['x', 'y', 'z1']]
  ero_df = ero_df.rename(columns={'x':'x', 'y':'y', 'z1':'z'})
  Z_grid_ero = make_Zgrid(ero_df)

  
  swe_df = results_df[['x', 'y', 'z2']]
  swe_df = swe_df.rename(columns={'x':'x', 'y':'y', 'z2':'z'})
  Z_grid_swe = make_Zgrid(swe_df)

  scenario = os.path.split(parDIR)[-1]

  ero_tif_f = os.path.join(outDIR, scenario + '_ero.tif')
  swe_tif_f = os.path.join(outDIR, scenario + '_swe.tif')

  make_gtif(Z_grid_ero, ero_tif_f)
  make_gtif(Z_grid_swe, swe_tif_f)



