# Import a bunch of python sub packages
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import xray
import xray.ufuncs as xu
import dask
#import seaborn as sn
from datetime import datetime
from dask.diagnostics import ProgressBar
import warnings
import os
import pandas as pd
from datetime import datetime

warnings.filterwarnings('ignore')

# Create a file list of all the netCDF files
import glob
fileList = glob.glob('/datadrive/WRF/radoutd02.nc')
fileList.sort()
fileList

# This cell converts concatenates the netCDF files and creates new ones that are better organized

# Create a function that will take the way the time is stored as a variable in the
# WRF netCDF files and decode it
# Build a function to go from numpy.bytes_ to a string - WRF files were compressed
def decode(d):
    decoded = datetime.strptime(d.decode(encoding='UTF-8'),  "%Y-%m-%d_%H:%M:%S")
    return decoded

first = True # create a boolean for the loop below


# Stores the creates a new netCDF files with coordinat
for file in fileList:
    print(file)
    ds = xray.open_dataset(file , engine = 'scipy') # may not need enginer = 'netcdf4' others are using engine = 'scipy'
    dates = ds.Times.to_series().apply(decode) # this creates a list of the varaibles Times in the netCDF file ds and then applies the decode function to each element
    dsTotal = xray.Dataset({'prec': (['time','x','y'], ds.prehourly.values),
                            'temp2m':(['time','x','y'], ds.t2.values)},
                           coords = {'longitude': (['x','y'], ds.lon.values),
                                     'latitude': (['x','y'], ds.lat.values),
                                     'time': xray.DataArray.from_series(dates).values})
    # add attributes to the netCDF
    dsTotal.attrs['prec'] = 'Precipitation Hourly [mm]'
    dsTotal.attrs['temp2m'] = 'Two Meter Temperature [deg K]'
    # write the netcdf files back to the drive, using the year-month as the name
    dsTotal.to_netcdf('/datadrive/WRF/temp/' + file.split('_')[2], mode = 'w')
print('done')


df = pd.read_csv('DHSVM_example.txt', sep = '\t')
names = ['date_time', 'temp2m', 'wind2m', 'RH', 'SW', 'LW', 'Precip']
df.columns = names
decoded_test = datetime.strptime("10/01/2002-01", "%m/%d/%Y-%H")
