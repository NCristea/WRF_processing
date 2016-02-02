# Import a bunch of python sub packages
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import xarray
#import seaborn as sn
from datetime import datetime
from dask.diagnostics import ProgressBar
import warnings
import os
import pandas as pd
from datetime import datetime

warnings.filterwarnings('ignore')

# Nicoleta path on local computer: /Users/carina/desktop/WRF_data
# Create a file list of all the netCDF files
import glob
fileList = glob.glob('/Users/carina/desktop/WRF_data/*.nc')
fileList.sort()
#fileList

def clean_netCDF(fileList):

    def decode(d):
        decoded = datetime.strptime(d.decode(encoding='UTF-8'),  "%Y-%m-%d_%H:%M:%S")
        return decoded

# Stores the creates a new netCDF files with coordinat
    for file in fileList:
        print(file)   
        ds = xarray.open_dataset(file , engine = 'scipy') # may not need engine = 'netcdf4' others are using engine = 'scipy'
        dates = ds.Times.to_series().apply(decode) # this creates a list of the varaibles Times in the netCDF file ds and then applies the decode function to each element
        dsTotal = xarray.Dataset({'prec': (['time','x','y'], ds.prehourly.values),
                  'temp2m':(['time','x','y'], ds.t2.values)},
                   coords = {'longitude': (['x','y'], ds.lon.values),
                          'latitude': (['x','y'], ds.lat.values),
                          'time': xarray.DataArray.from_series(dates).values})    
    # add attributes to the netCDF
        dsTotal.attrs['prec'] = 'Precipitation Hourly [mm]'
        dsTotal.attrs['temp2m'] = 'Two Meter Temperature [deg K]'
    # write the netcdf files back to the drive, using the year-month as the name
        dsTotal.to_netcdf('/Users/carina/desktop/WRF_data/temp/' + file.split('_')[6], mode = 'w')
    print('done cleaning files')


# call this once to create temporary cleaned netCDF files

clean_netCDF(fileList)

dsTotal = xarray.open_mfdataset('/Users/carina/desktop/WRF_data/temp/*.nc')
dsTotal.chunk({'time':400,'x':50,'y':50})

#this loops through the lon list and back calculates the x and y indices needed to plot or extract data

bb = {'minLong':-119.40, 'maxLong':-119.20, 'minLat': 37.73, 'maxLat':37.96}

long = dsTotal.coords['longitude'].values
lat = dsTotal.coords['latitude'].values

xycord = np.where( (long > bb['minLong'] ) & (long < bb['maxLong']) & (lat > bb['minLat']) & (lat < bb['maxLat'])) 
xcord = xycord[:][0] 
ycord = xycord[:][1]

all_Temp = np.zeros((dsTotal.time.shape[0], xcord.shape[0]))

for i in range(xcord.shape[0]):
    selectTemp = dsTotal.sel_points(x = [xcord[i]], y = [ycord[i]]).temp2m
    all_Temp[:, i] = selectTemp

# convert to pandas
df_Temp = pd.DataFrame(all_Temp)
# add time variable to data frame
df_Time = pd.DataFrame(dsTotal.time.values)

df_TimeTemp = pd.concat([df_Time, df_Temp], axis = 1)

df = pd.read_csv('DHSVM_example.txt', sep = '\t')
names = ['date_time', 'temp2m', 'wind2m', 'RH', 'SW', 'LW', 'Precip']
df.columns = names
decoded_test = datetime.strptime("10/01/2002-01", "%m/%d/%Y-%H")

#test changes