#this script will read all the WRF radout and surfout files and will merge 6 variables and concatenate along the time dimension
#need to update the names of the WRF model and output file for each run of the script

#Import Python libraries

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import xarray as xray
import xarray
#import seaborn as sn
from datetime import datetime
from dask.diagnostics import ProgressBar
import warnings
import os
import pandas as pd
from datetime import datetime
import re

path = '/Users/carina/Desktop/WRF_data/'
#filePos = 6
# Anthony's path
#path = 'c:\\work\\datadrive\\WRF\\'
#filePos = 3

# Create a file list of all the netCDF files
import glob
fileList = glob.glob(path + '*.nc')
fileList.sort()

files_dict = {}

for file in fileList:
    radout = re.search('radout', file) is not None
    surfout = re.search('surfout', file) is not None
    if not (radout or surfout):
        continue
    #NARR = re.search('NARR', file) != None
    #Morr = re.search('Morr', file) != None

    year = file[len(file)-10:len(file)-6]
    month = file[len(file)-5:len(file)-3]

    index = "{}{}".format(year, month)
    if index in files_dict:
        file_dict = files_dict.get(index)
    else:
        file_dict = {}

    if radout:
        file_dict['radout'] = file
    if surfout:
        file_dict['surfout'] = file

    files_dict[index] = file_dict


def decode(d):
    decoded = datetime.strptime(d.decode(encoding='UTF-8'),  "%Y-%m-%d_%H:%M:%S")
    return decoded


def process_month(radFile, surfFile, month):
    radDS = xarray.open_dataset(radFile, engine = 'scipy') # may not need engine = 'netcdf4' others are using engine = 'scipy'
    radDates = radDS.Times.to_series().apply(decode) # this creates a list of the varaibles Times in the netCDF file ds and then applies the decode function to each element
    dsTotal_rad = xarray.Dataset({'SWdown': (['time','x','y'], radDS.swdnbhourly.values),
                                  'LWdown':(['time','x','y'], radDS.lwdnbhourly.values)},
                                 coords = {'longitude': (['x','y'], radDS.lon.values),
                                           'latitude': (['x','y'], radDS.lat.values),
                                           'time': xarray.DataArray.from_series(radDates).values})

    surfDS = xarray.open_dataset(surfFile , engine = 'scipy') # may not need engine = 'netcdf4' others are using engine = 'scipy'
    surfDates = surfDS.Times.to_series().apply(decode) # this creates a list of the varaibles Times in the netCDF file ds and then applies the decode function to each element
    #calculate 2m wind speed from the 10 m u10 and v10
    wind2m = 0.54 * np.sqrt(surfDS.u10.values * surfDS.u10.values + surfDS.v10.values * surfDS.v10.values)
    dsTotal_surf = xarray.Dataset({'prec': (['time','x','y'], surfDS.prehourly.values),
                              'temp2m':(['time','x','y'], surfDS.t2.values),
                              'rh2m':(['time','x','y'], surfDS.rh2.values),
                              'wind2m':(['time','x','y'], wind2m)},
                             coords = {'longitude': (['x','y'], surfDS.lon.values),
                                       'latitude': (['x','y'], surfDS.lat.values),
                                       'time': xarray.DataArray.from_series(surfDates).values})
    dsTotal_surf = dsTotal_surf.merge(dsTotal_rad)
    # first the prec and temp2m variables are created
    # add attributes to the netCDF
    dsTotal_surf.attrs['prec'] = 'Precipitation Hourly [mm]'
    dsTotal_surf.attrs['temp2m'] = 'Two Meter Temperature [deg K]'
    dsTotal_surf.attrs['rh2m'] = 'Two Meter Relative Humidity [%?]'
    dsTotal_surf.attrs['wind2m'] = 'Two Meter Wind Speed [m/s]'
    dsTotal_surf.attrs['SWdown'] = 'Shortwave radiation [W/m2]'
    dsTotal_surf.attrs['LWdown'] = 'Incoming longwave radiation [W/m2]'

    # write the netcdf files back to the drive, using the year-month as the name

    #dsTotal_surf.to_netcdf(path + 'temp1/' + month + '.nc', mode = 'w')

    #print('done cleaning files')

    return dsTotal_surf

ds_total = None

for k in sorted(files_dict):
    print("Processing month {}".format(k))
    v = files_dict[k]
    radout_file = v['radout']
    surfout_file = v['surfout']
    ds = process_month(radout_file, surfout_file, k)

    if ds_total:
        ds_total = xray.concat([ds_total, ds], dim = 'time')
    else:
        ds_total = ds


# ds_total.to_netcdf(path + 'ds_total_NARR_Morr')

#subsetting the WRF dataset for the upper Tuolumne area

bb = {'minLong':-119.40, 'maxLong':-119.20, 'minLat': 37.73, 'maxLat':37.96}

long = ds_total.coords['longitude'].values
lat = ds_total.coords['latitude'].values

xycord = np.where( (long > bb['minLong'] ) & (long < bb['maxLong']) & (lat > bb['minLat']) & (lat < bb['maxLat']))
xcord = xycord[:][0]
ycord = xycord[:][1]

reduced_ds = ds_total.isel(x=xcord, y=ycord)
# df = reduced_ds.to_dataframe()
# df.to_csv(path + 'reduced.csv')
#reduced_ds.to_netcdf(path + 'ds_reduced_NARR_Morr.nc', format='NETCDF4', mode='w')
reduced_ds.to_netcdf(path + 'ds_reduced_NARR_Morr.nc')


'''
#this script will read all the WRF radout and surfout files and will merge 6 variables and concatenate along the time dimension
#need to update the names of the WRF model and output file for each run of the script

#Import Python libraries

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import xarray as xray
import xarray
#import seaborn as sn
from datetime import datetime
from dask.diagnostics import ProgressBar
import warnings
import os
import pandas as pd
from datetime import datetime
import re

path = '/Users/carina/desktop/WRF_data/'
#filePos = 6
# Anthony's path
#path = 'c:\\work\\datadrive\\WRF\\'
#filePos = 3

# Create a file list of all the netCDF files
import glob
fileList = glob.glob(path + '*.nc')
# fileList.sort()

files_dict = {}

for file in fileList:
    radout = re.search('radout', file) != None
    surfout = re.search('surfout', file) != None
    #NARR = re.search('NARR', file) != None
    #Morr = re.search('Morr', file) != None

    year = file[len(file)-10:len(file)-6]
    month = file[len(file)-5:len(file)-3]

    index = "{}{}".format(year, month)
    if index in files_dict:
        file_dict = files_dict.get(index)
    else:
        file_dict = {}

    if radout:
        file_dict['radout'] = file
    if surfout:
        file_dict['surfout'] = file

    files_dict[index] = file_dict


def decode(d):
    decoded = datetime.strptime(d.decode(encoding='UTF-8'),  "%Y-%m-%d_%H:%M:%S")
    return decoded


def process_month(radFile, surfFile):
    radDS = xarray.open_dataset(radFile, engine = 'scipy') # may not need engine = 'netcdf4' others are using engine = 'scipy'
    radDates = radDS.Times.to_series().apply(decode) # this creates a list of the varaibles Times in the netCDF file ds and then applies the decode function to each element
    dsTotal_rad = xarray.Dataset({'SWdown': (['time','x','y'], radDS.swdnbhourly.values),
                                  'LWdown':(['time','x','y'], radDS.lwdnbhourly.values)},
                                 coords = {'longitude': (['x','y'], radDS.lon.values),
                                           'latitude': (['x','y'], radDS.lat.values),
                                           'time': xarray.DataArray.from_series(radDates).values})

    surfDS = xarray.open_dataset(file , engine = 'scipy') # may not need engine = 'netcdf4' others are using engine = 'scipy'
    surfDates = surfDS.Times.to_series().apply(decode) # this creates a list of the varaibles Times in the netCDF file ds and then applies the decode function to each element
    #calculate 2m wind speed from the 10 m u10 and v10
    wind2m = 0.54 * np.sqrt(surfDS.u10.values * surfDS.u10.values + surfDS.v10.values * surfDS.v10.values)
    dsTotal_surf = xarray.Dataset({'prec': (['time','x','y'], surfDS.prehourly.values),
                                   'temp2m':(['time','x','y'], surfDS.t2.values),
                                   'rh2m':(['time','x','y'], surfDS.rh2.values),
                                   'wind2m':(['time','x','y'], wind2m)},
                                  coords = {'longitude': (['x','y'], surfDS.lon.values),
                                            'latitude': (['x','y'], surfDS.lat.values),
                                            'time': xarray.DataArray.from_series(surfDates).values})
    dsTotal_surf = dsTotal_surf.merge(dsTotal_rad)
    # first the prec and temp2m variables are created
    # add attributes to the netCDF
    dsTotal_surf.attrs['prec'] = 'Precipitation Hourly [mm]'
    dsTotal_surf.attrs['temp2m'] = 'Two Meter Temperature [deg K]'
    dsTotal_surf.attrs['rh2m'] = 'Two Meter Relative Humidity [%?]'
    dsTotal_surf.attrs['wind2m'] = 'Two Meter Wind Speed [m/s]'
    dsTotal_surf.attrs['SWdown'] = 'Shortwave radiation [W/m2]'
    dsTotal_surf.attrs['LWdown'] = 'Incoming longwave radiation [W/m2]'

    # write the netcdf files back to the drive, using the year-month as the name
    #dsTotal.to_netcdf(path + 'temp1/' + file.split('_')[filePos], mode = 'w')
    #print('done cleaning files')

    return dsTotal_surf

ds_total = None

for k in sorted(files_dict):
    v = files_dict[k]
    radout_file = v['radout']
    surfout_file = v['surfout']
    ds = process_month(radout_file, surfout_file)
    if ds_total:
        ds_total = xray.concat([ds_total, ds], dim = 'time')
    else:
        ds_total = ds

# ds_total.to_netcdf(path + 'ds_total_NARR_Morr')

#subsetting the WRF dataset for the upper Tuolumne area
bb = {'minLong':-119.40, 'maxLong':-119.20, 'minLat': 37.73, 'maxLat':37.96}

long = ds_total.coords['longitude'].values
lat = ds_total.coords['latitude'].values

xycord = np.where( (long > bb['minLong'] ) & (long < bb['maxLong']) & (lat > bb['minLat']) & (lat < bb['maxLat']))
xcord = xycord[:][0]
ycord = xycord[:][1]

reduced_ds = ds_total.isel(x=xcord, y=ycord)
reduced_ds.to_netcdf(path + 'ds_reduced_NARR_Morr')'''