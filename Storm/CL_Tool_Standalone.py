
####
'GDAL requires Python 3.9 minimum.'
'Same as CL_Tool.py script except it doesnt create an exe.'
'Exe or py script need to be run from top of directory structure.'
'Modify list.txt with query info with columns: id, lat, lon, year window, gcm name, wind option.'
'VALID YEAR WINDOW STRINGS: 1974_2013, 2000_2029, 2010_2039, 2020_2049, 2030_2059, 2040_2069, 2050_2079, 2060_2089, 2070_2099.'
'VALID GCM NAME STRINGS: CCSM4, CanESM2, MIROC5.'
'id can be string of own choosing, but valid year windows and gcms must be used.'
'There are two wind options to get wind parameters, since wind data isnt included in the GDBs.'
'Wind option takes data from nearest CLIGEN ground net station with wind data when set to Search.'
'Wind option takes data from formatted wind string copied to txt file in /wind-strings when set to the name of the txt file (without .txt extension).'
####

from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import pandas as pd
from classes.Formatting import Formatting
formatting_obj = Formatting()
from osgeo import ogr
from osgeo import gdal


cwd = os.getcwd()

parDIR = os.path.join(cwd, 'pars')
parnetDIR = os.path.join(cwd, 'parfiles-2015')
windDIR = os.path.join(cwd, 'wind-strings')
listFILE = os.path.join(cwd, 'list.txt')
searchFILE = os.path.join(cwd, 'search_stations.txt')


def readlist(listFILE):
  with open(listFILE) as f:
    next(f)
    lines = f.readlines()
  points = []
  for line in lines:
    row = line.strip('\n').split(',')
    iD = str(row[0])
    y = float(row[2])
    x = float(row[1])
    yrs = str(row[3])
    gcm = str(row[4])
    wind = str(row[5])
    points.append([iD, y, x, yrs, gcm, wind])
  return points

def readsearchlist(searchFILE):  
  with open(os.path.join(cwd, searchFILE)) as f:
    next(f)
    slines = f.readlines()
    spoints = []
    for line in slines:
      row = line.strip('\n').split(',')
      iD = str(row[0])
      y = float(row[2])
      x = float(row[1])
      spoints.append([iD, y, x])
  return spoints

def readwindstring(parname):
  with open(os.path.join(parnetDIR, parname + '.par')) as f:
    lines = f.readlines()
  lines = [l.strip('\n') for l in lines[17:]]
  wind_lines = []
  for i, line in enumerate(lines):
    wind_lines.append(line)
    if '---Wind Stations---' in line:
      wind_lines.append(lines[i+1])
      break
  return wind_lines

def readcustomwindstring(windname):
  with open(os.path.join(windDIR, windname + '.txt')) as f:
    wlines = f.readlines()
  wind_lines = [l.strip('\n') for l in wlines]
  return wind_lines

def coord2pixel(x, y, trans):
  xpt = (((x - trans[0]) / trans[1]))
  ypt = (((y - trans[3]) / trans[5]))
  return(xpt, ypt)

def find_closest_point(lat, lon, spoints):
  point1 = ogr.Geometry(ogr.wkbPoint)
  point1.AddPoint(lon, lat)
  point2 = ogr.Geometry(ogr.wkbPoint)
  min_distance = 99999999.9
  for spoint in spoints:
    slat, slon = spoint[1], spoint[2]
    point2.AddPoint(slon, slat)
    distance = point1.Distance(point2)
    if distance < min_distance:
      min_distance = distance
      closest_point = spoint[0]
  return closest_point



def main(point):

  var_labels = ['mean', 'sdev', 'skew', 'pww', 'pwd', 'tmax', 'tmin', 'txsd', 'tnsd', 'srad', 'srsd', 'mx5p', 'tdew', 'timepk']
  historical_var_labels = ['timepk']
  fname = str(point[0])
  yr_str = str(point[3])
  lat, lon = point[2], point[1]
  wind = str(point[5])

  closest_point = find_closest_point(lat, lon, spoints)

  if wind == 'Search':
    wlines = readwindstring(closest_point)
    wind_str = '\n'.join(l for l in wlines) + '\n'
  else:
    wlines = readcustomwindstring(wind)
    wind_str = '\n'.join(l for l in wlines) + '\n'

  par_df = pd.DataFrame(columns=range(1, 12), index=var_labels)

  for varlb in var_labels:
    
    if varlb not in historical_var_labels:
      raster_name = varlb + '_' + yr_str
    else:
      raster_name = varlb + '_1974_2013'
  
    raster = gdal.Open(f'OpenFileGDB:{gdbDIR}:{raster_name}')
    transform = raster.GetGeoTransform()
    xp, yp = coord2pixel(lon, lat, transform)
  
    for mo in range(1, 13):
      band = raster.GetRasterBand(mo)
      pixel_value = band.ReadAsArray(int(xp), int(yp), 1, 1)[0, 0] 
      par_df.at[varlb, mo] = pixel_value  
      band = None
    
    raster.FlushCache()
    raster = None
  
  raster_name = 'DEM'
  raster = gdal.Open(f'OpenFileGDB:{gdbDIR}:{raster_name}')
  band = raster.GetRasterBand(1)
  pixel_value = band.ReadAsArray(int(xp), int(yp), 1, 1)[0, 0] 
  dem_elev = pixel_value*3.28084
  band = None
  raster.FlushCache()
  raster = None
  
  par_df.loc['tmax'] = par_df.loc['tmax'].apply(lambda x: x*(9.0/5.0) + 32.0)
  par_df.loc['tmin'] = par_df.loc['tmin'].apply(lambda x: x*(9.0/5.0) + 32.0)
  par_df.loc['tdew'] = par_df.loc['tdew'].apply(lambda x: x*(9.0/5.0) + 32.0)
  par_df.loc['txsd'] = par_df.loc['txsd'].apply(lambda x: x*(9.0/5.0))
  par_df.loc['tnsd'] = par_df.loc['tnsd'].apply(lambda x: x*(9.0/5.0))
  par_df.loc['srad'] = par_df.loc['srad'].apply(lambda x: x*86400.0/41868.0)
  par_df.loc['srsd'] = par_df.loc['srsd'].apply(lambda x: x*86400.0/41868.0)
  
  meanP_list = par_df.loc['mean'].clip(lower=0.01)
  sdevP_list = par_df.loc['sdev'].clip(lower=0.01)
  skewP_list = par_df.loc['skew']
  ww_list = par_df.loc['pww'].clip(lower=0.01, upper=0.99)
  wd_list = par_df.loc['pwd'].clip(lower=0.01, upper=0.99)
  tmax_list = par_df.loc['tmax']
  tmin_list = par_df.loc['tmin']
  sdtmax_list = par_df.loc['txsd'].clip(lower=0.01)
  sdtmin_list = par_df.loc['tnsd'].clip(lower=0.01)
  solrad_list = par_df.loc['srad'].clip(lower=0.0)
  solsdev_list = par_df.loc['srsd'].clip(lower=0.01)
  mx5p_list = par_df.loc['mx5p'].clip(lower=0.01)
  dewpt_list = par_df.loc['tdew']
  timepk_list = par_df.loc['timepk'].sort_values().clip(lower=0.01, upper=1.00)
  
  meanP, sdevP, skewP, ww, wd = [], [], [], [], []
  tmax, tmin, sdtmax, sdtmin = [], [], [], []
  solrad, solsdev, mx5p, dewpt, timepk = [], [], [], [], []
  
  for mo in range(1, 13):
    meanP.append('{:.2f}'.format(meanP_list[mo]).lstrip('0'))
    sdevP.append('{:.2f}'.format(sdevP_list[mo]).lstrip('0'))
    skewP.append('{:.2f}'.format(skewP_list[mo]).lstrip('0'))
    ww.append('{:.2f}'.format(ww_list[mo]).lstrip('0'))
    wd.append('{:.2f}'.format(wd_list[mo]).lstrip('0'))
    tmax.append('{:.2f}'.format(tmax_list[mo]).lstrip('0'))
    tmin.append('{:.2f}'.format(tmin_list[mo]).lstrip('0'))
    sdtmax.append('{:.2f}'.format(sdtmax_list[mo]).lstrip('0'))
    sdtmin.append('{:.2f}'.format(sdtmin_list[mo]).lstrip('0')) 
    solrad.append('%.0f' %round(float(solrad_list[mo]), 0) + '.')
    solsdev.append('%.1f' %round(float(solsdev_list[mo]), 1))
    mx5p.append(('%.2f' %round(float(mx5p_list[mo]), 2)).lstrip('0'))
    dewpt.append('%.2f' %round(float(dewpt_list[mo]), 2))
    timepk.append(('%.3f' %round(float(timepk_list[mo]), 3)).lstrip('0'))
    
  with open(os.path.join(parDIR, fname + '.par'), 'w') as f_out:
  
    title = ' Grid point'
    yrsStr = 'YEARS= 30. TYPE= 1'
    tpStr = 'TP5 = 1.23 TP6= 2.34' 
  
    lat, lon = '%.2f' % round(float(lat), 2), '%.2f' % round(float(lon), 2)
    elev = '%.0f' % round(float(dem_elev), 0) + '.'
  
    f_out.write(title + '\n' )
    f_out.write(' LATT=' + formatting_obj.spacing_lat_lon(lat) + lat)
    f_out.write(' LONG=' + formatting_obj.spacing_lat_lon(lon) + lon)
    f_out.write(' ' + yrsStr + '\n')
    f_out.write(' ELEVATION =' + formatting_obj.spacing_elev(elev) + elev)
    f_out.write(' ' + tpStr + '\n')
  
    gen = formatting_obj.spacing_gen(meanP)
    f_out.write(' MEAN P ' + ''.join([next(gen) + str(x) + next(gen) if i < 1 else str(x) if i == 11 else str(x) + next(gen) for i, x in enumerate(meanP)]) + '\n')
    gen = formatting_obj.spacing_gen(sdevP)
    f_out.write(' S DEV P' + ''.join([next(gen) + str(x) + next(gen) if i < 1 else str(x) if i == 11 else str(x) + next(gen) for i, x in enumerate(sdevP)]) + '\n')
    gen = formatting_obj.spacing_gen(skewP)
    f_out.write(' SKEW  P' + ''.join([next(gen) + str(x) + next(gen) if i < 1 else str(x) if i == 11 else str(x) + next(gen) for i, x in enumerate(skewP)]) + '\n')
    gen = formatting_obj.spacing_gen(ww)
    f_out.write(' P(W/W) ' + ''.join([next(gen) + str(x) + next(gen) if i < 1 else str(x) if i == 11 else str(x) + next(gen) for i, x in enumerate(ww)]) + '\n')
    gen = formatting_obj.spacing_gen(wd)
    f_out.write(' P(W/D) ' + ''.join([next(gen) + str(x) + next(gen) if i < 1 else str(x) if i == 11 else str(x) + next(gen) for i, x in enumerate(wd)]) + '\n')
    gen = formatting_obj.spacing_gen(tmax)
    f_out.write(' TMAX AV' + ''.join([next(gen) + str(x) + next(gen) if i < 1 else str(x) if i == 11 else str(x) + next(gen) for i, x in enumerate(tmax)]) + '\n')
    gen = formatting_obj.spacing_gen(tmin)
    f_out.write(' TMIN AV' + ''.join([next(gen) + str(x) + next(gen) if i < 1 else str(x) if i == 11 else str(x) + next(gen) for i, x in enumerate(tmin)]) + '\n')
    gen = formatting_obj.spacing_gen(sdtmax)
    f_out.write(' SD TMAX' + ''.join([next(gen) + str(x) + next(gen) if i < 1 else str(x) if i == 11 else str(x) + next(gen) for i, x in enumerate(sdtmax)]) + '\n')
    gen = formatting_obj.spacing_gen(sdtmin)
    f_out.write(' SD TMIN' + ''.join([next(gen) + str(x) + next(gen) if i < 1 else str(x) if i == 11 else str(x) + next(gen) for i, x in enumerate(sdtmin)]) + '\n')
    gen = formatting_obj.spacing_gen(solrad)
    f_out.write( ' SOL.RAD' + ''.join([next(gen) + str(x) + next(gen) if i < 1 else str(x) if i == 11 else str(x) + next(gen) for i, x in enumerate(solrad)]) + '\n')
    gen = formatting_obj.spacing_gen(solsdev)
    f_out.write( ' SD SOL ' + ''.join([next(gen) + str(x) + next(gen) if i < 1 else str(x) if i == 11 else str(x) + next(gen) for i, x in enumerate(solsdev)]) + '\n')
    gen = formatting_obj.spacing_gen(mx5p)
    f_out.write( ' MX .5 P' + ''.join([next(gen) + str(x) + next(gen) if i < 1 else str(x) if i == 11 else str(x) + next(gen) for i, x in enumerate(mx5p)]) + '\n')
    gen = formatting_obj.spacing_gen(dewpt)
    f_out.write( ' DEW PT ' + ''.join([next(gen) + str(x) + next(gen) if i < 1 else str(x) if i == 11 else str(x) + next(gen) for i, x in enumerate(dewpt)]) + '\n')
    gen = formatting_obj.spacing_gen(timepk)
    f_out.write( 'Time Pk ' + ''.join([next(gen) + str(x) + next(gen) if i < 1 else str(x) if i == 11 else str(x) + next(gen) for i, x in enumerate(timepk)]) + '\n')
    
    f_out.write(wind_str)
                 
    if wind != 'Search':
      f_out.write('---\n')


points = readlist(listFILE)
spoints = readsearchlist(searchFILE)
n_workers = 100
spoints = readsearchlist(searchFILE)    
gcm_name = points[0][4]
gdbDIR = os.path.join(cwd, '{}.gdb'.format(gcm_name))

if __name__ == '__main__':

  with ProcessPoolExecutor(max_workers=n_workers) as executor:
    pool = {executor.submit(main, p): p for p in points}
    ct = 0
    for future in as_completed(pool):
      res = future.result()
      ct += 1
      if ct % 100000 == 0:
        print(ct)

