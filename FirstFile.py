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

path = '/Users/carina/desktop/WRF_data/'
filePos = 6
# Anthony's path
#path = 'c:\\work\\datadrive\\WRF\\'
#filePos = 3

# Create a file list of all the netCDF files
import glob
fileList = glob.glob(path + '*.nc')
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
        #calculate 2m wind speed from the 10 m u10 and v10
        wind2m = 0.74 * np.sqrt(ds.u10.values * ds.u10.values + ds.v10.values * ds.v10.values)
        dsTotal = xarray.Dataset({'prec': (['time','x','y'], ds.prehourly.values),
                  'temp2m':(['time','x','y'], ds.t2.values),
                  'rh2m':(['time','x','y'], ds.rh2.values),
                  'wind2m':(['time','x','y'], wind2m)}
        coords = {'longitude': (['x','y'], ds.lon.values),
                          'latitude': (['x','y'], ds.lat.values),
                          'time': xarray.DataArray.from_series(dates).values})

    # add attributes to the netCDF
        dsTotal.attrs['prec'] = 'Precipitation Hourly [mm]'
        dsTotal.attrs['temp2m'] = 'Two Meter Temperature [deg K]'
        dsTotal.attrs['rh2m'] = 'Two Meter Relative Humidity [%?]'
        dsTotal.attrs['wind2m'] = 'Two Meter Wind Speed [m/s]'
    # write the netcdf files back to the drive, using the year-month as the name
        dsTotal.to_netcdf(path + 'temp/' + file.split('_')[filePos], mode = 'w')
    print('done cleaning files')


# call this once to create temporary cleaned netCDF files

clean_netCDF(fileList)

dsTotal = xarray.open_mfdataset(path + 'temp/*.nc')
dsTotal.chunk({'time':400,'x':50,'y':50})

#this loops through the lon list and back calculates the x and y indices needed to plot or extract data

bb = {'minLong':-119.40, 'maxLong':-119.20, 'minLat': 37.73, 'maxLat':37.96}

long = dsTotal.coords['longitude'].values
lat = dsTotal.coords['latitude'].values

xycord = np.where( (long > bb['minLong'] ) & (long < bb['maxLong']) & (lat > bb['minLat']) & (lat < bb['maxLat'])) 
xcord = xycord[:][0] 
ycord = xycord[:][1]

all_Temp = np.zeros((dsTotal.time.shape[0], xcord.shape[0]))
all_Prec = np.zeros((dsTotal.time.shape[0], xcord.shape[0]))

columnNames = ['date_time']

for i in range(xcord.shape[0]):
    selectTemp = dsTotal.sel_points(x = [xcord[i]], y = [ycord[i]]).temp2m
    selectPrec = dsTotal.sel_points(x = [xcord[i]], y = [ycord[i]]).prec
    all_Temp[:, i] = selectTemp
    all_Prec[:, i] = selectPrec
    columnNames.append(str(selectTemp.coords['latitude'].values[0]) + '_' + str(selectTemp.coords['longitude'].values[0]))

# convert to pandas
df_Temp = pd.DataFrame(all_Temp)
df_Prec = pd.DataFrame(all_Prec)
# add time variable to data frame
df_Time = pd.DataFrame(dsTotal.time.values)

df_TimeTemp = pd.concat([df_Time, df_Temp], axis = 1)
df_TimePrec = pd.concat([df_Time, df_Prec], axis = 1)
# add names to the columns
df_TimeTemp.columns = columnNames
df_TimePrec.columns = columnNames

df_TimeTemp.set_index('date_time', inplace = True)
df_TimePrec.set_index('date_time', inplace = True)

df = pd.read_csv(path + 'DHSVM_example.txt', sep = '\t')
names = ['date_time', 'temp2m', 'wind2m', 'RH', 'SW', 'LW', 'Precip']
df.columns = names
df['date_time'] = pd.to_datetime(df['date_time'], format = "%m/%d/%Y-%H")
df.set_index('date_time', inplace = True)

for i in range(xcord.shape[0]):
    node_prec = df_TimePrec[df_TimePrec.columns[i]]
    names_node_prec = ['date_time', 'Precip']
    node_prec.columns = names_node_prec
    df.update(node_prec, join='left', overwrite=True)
    fileName = 'DHSVM_prec_Input_' + columnNames[i+1]
    df.to_csv(path + fileName, sep='\t')

#concatenate station data to the wrf data frame to explore
#all_T = pd.concat([df_TimeTemp-273.15, df['temp2m']], axis = 1, join_axes = [df_TimeTemp.index])
#all_P = pd.concat([df_TimePrec/1000, df['Precip']], axis = 1, join_axes = [df_TimeTemp.index])

#all_P.plot()

#example of how to explore the data
#all_T.resample('D', how = 'mean').plot()
#df['temp2m'].groupby(df.index.year).mean().plot()

#node_prec = df_TimePrec['37.734_-119.348'];
#names_node_prec = ['date_time', 'Precip']
#node_prec.columns = names_node_prec