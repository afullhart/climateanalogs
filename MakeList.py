import os

points_list = r'/home/afullhart/Downloads/GDBs/RCP45/Points.csv'
par_folder = r'/home/afullhart/Downloads/GDBs/RCP45/pars-lists'
list_hdr = 'id,lat_dd,lon_dd,yr_window,gcm_name,wind_str_option'

yr_windows = ['1974_2013', '2000_2029', '2010_2039', '2020_2049', '2030_2059', '2040_2069', '2050_2079', '2060_2089', '2070_2099']
gcm_names = ['CCSM4', 'CanESM2', 'MIROC5']



with open(points_list) as f:
  lineOne = f.readline()
  lineTwo = f.readline()

print(lineOne)
print(lineTwo)

with open(points_list) as f:
  next(f)
  listLines = f.readlines()

p_data = []
for l in listLines:
  row = l.strip('\n').split(',')
  iD = row[0]
  lon = row[1]
  lat = row[2]
  p_data.append([iD, lon, lat])



for i in gcm_names:

  for j in yr_windows:

    with open(os.path.join(par_folder, i + '_' + j + '_list.txt'), 'w') as f_out:
      f_out.write(list_hdr + '\n')
      for data in p_data:
        iD = data[0]
        lon = data[1]
        lat = data[2]
        f_out.write(','.join([iD, lat, lon, j, i, 'Search']) + '\n')
      
