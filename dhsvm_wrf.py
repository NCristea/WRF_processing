import pandas as pd
import xarray
import numpy as np
from dask.diagnostics import ProgressBar
from pytz import timezone

western = timezone('US/Pacific')
utc = timezone('UTC')

path = '/Users/carina/desktop/WRF_data/'
#current WRF runs - combination of microphysics and boundary conditions

#open the reduced dataset
# dsTotal = xarray.open_dataset('/Users/carina/desktop/WRF_data/ds_reduced_NARR_Morr.nc', engine = 'scipy')
dsTotal = xarray.open_dataset('/Users/carina/desktop/WRF_data/ds_reduced_NARR_Morr.nc', engine = 'netcdf4')
#dsTotal.chunk({'time':400,'x':50,'y':50})
#converts the dataset into pandas and numpy arrays

xcord = dsTotal.x.to_series().values
ycord = dsTotal.y.to_series().values

#initialize numpy arrays
all_Temp = np.zeros((dsTotal.time.shape[0], xcord.shape[0]))
all_Prec = np.zeros((dsTotal.time.shape[0], xcord.shape[0]))
all_Wind = np.zeros((dsTotal.time.shape[0], xcord.shape[0]))
all_RH = np.zeros((dsTotal.time.shape[0], xcord.shape[0]))
all_SW = np.zeros((dsTotal.time.shape[0], xcord.shape[0]))
all_LW = np.zeros((dsTotal.time.shape[0], xcord.shape[0]))

columnNames = ['date_time']

for i in range(xcord.shape[0]):
    selectTemp = dsTotal.isel_points(x = [i], y = [i]).temp2m
    selectPrec = dsTotal.isel_points(x = [i], y = [i]).prec
    selectWind = dsTotal.isel_points(x = [i], y = [i]).wind2m
    selectRH = dsTotal.isel_points(x = [i], y = [i]).rh2m
    selectSW = dsTotal.isel_points(x = [i], y = [i]).SWdown
    selectLW = dsTotal.isel_points(x = [i], y = [i]).LWdown
    all_Temp[:, i] = selectTemp - 273.15 #convert to Celsius
    all_Prec[:, i] = selectPrec/1000 # convert to m
    all_Wind[:, i] = selectWind
    all_RH[:, i] = selectRH
    all_SW[:, i] = selectSW
    all_LW[:, i] = selectLW
    #assign names to different columns as a function of the WRF node location (lat-long is in the name)
    columnNames.append(str(selectTemp.coords['latitude'].values[0]) + '_' + str(selectTemp.coords['longitude'].values[0]))

# convert the numpy arrays to pandas dataframes
df_Temp = pd.DataFrame(all_Temp)
df_Prec = pd.DataFrame(all_Prec)
df_Wind = pd.DataFrame(all_Wind)
df_RH = pd.DataFrame(all_RH)
df_SW = pd.DataFrame(all_SW)
df_LW = pd.DataFrame(all_LW)


#get the time variable from the original netCDF file
df_Time = pd.DataFrame(dsTotal.time.values)
dtIndex = pd.DatetimeIndex(dsTotal.time.values)
freq = pd.infer_freq(dtIndex);
#df_Time = df_Time.set_index(dtIndex)
# this is slow
#df_Time.apply(tz_update_utc)

df_TimeTemp = pd.concat([df_Time, df_Temp], axis = 1)
df_TimePrec = pd.concat([df_Time, df_Prec], axis = 1)
df_TimeWind = pd.concat([df_Time, df_Wind], axis = 1)
df_TimeRH = pd.concat([df_Time, df_RH], axis = 1)
df_TimeSW = pd.concat([df_Time, df_SW], axis = 1)
df_TimeLW = pd.concat([df_Time, df_LW], axis = 1)

# add time variable to data frame
# add names to the columns
df_TimeTemp.columns = columnNames
df_TimePrec.columns = columnNames
df_TimeWind.columns = columnNames
df_TimeRH.columns = columnNames
df_TimeSW.columns = columnNames
df_TimeLW.columns = columnNames

#the newly created dataframes can't read the time values-> need to set the time as index
df_TimeTemp.set_index(pd.DatetimeIndex(df_TimeTemp['date_time'], tz=utc, freq=freq).tz_convert(western), inplace=True)
df_TimePrec.set_index(pd.DatetimeIndex(df_TimePrec['date_time'], tz=utc, freq=freq).tz_convert(western), inplace=True)
df_TimeWind.set_index(pd.DatetimeIndex(df_TimeWind['date_time'], tz=utc, freq=freq).tz_convert(western), inplace=True)
df_TimeRH.set_index(pd.DatetimeIndex(df_TimeRH['date_time'], tz=utc, freq=freq).tz_convert(western), inplace=True)
df_TimeSW.set_index(pd.DatetimeIndex(df_TimeSW['date_time'], tz=utc, freq=freq).tz_convert(western), inplace=True)
df_TimeLW.set_index(pd.DatetimeIndex(df_TimeLW['date_time'], tz=utc, freq=freq).tz_convert(western), inplace=True)


#open the station data file and stich WRF data
df = pd.read_csv(path + 'DHSVM_example.txt', sep = '\t')
names = ['date_time', 'temp2m', 'wind2m', 'RH', 'SW', 'LW', 'Precip']
df.columns = names
'''
The dates in DHSVM are not taking into account daylight savings. Because of that
we can't use the timezone for the Datetime index.
To work around this we convert the UTC timezone
used in netCDF to US/Pacific...
'''
df['date_time'] = pd.to_datetime(df['date_time'], format="%m/%d/%Y-%H")
df.set_index(pd.DatetimeIndex(df['date_time']), inplace=True)


#open the station data file and stich WRF data
#df = pd.read_csv(path + 'DHSVM_example.txt', sep = '\t')
#names = ['date_time', 'temp2m', 'wind2m', 'RH', 'SW', 'LW', 'Precip']
#df.columns = names
#df['date_time'] = pd.to_datetime(df['date_time'], format = "%m/%d/%Y-%H")
#df.set_index('date_time', inplace = True)

#take a quick look at all the T and P data - how do they compare with the station data

#concatenate station data to the wrf data frame to explore
all_T_with_station = pd.concat([df_TimeTemp, df['temp2m']], axis = 1, join_axes = [df_TimeTemp.index])
all_P_with_station = pd.concat([df_TimePrec, df['Precip']], axis = 1, join_axes = [df_TimeTemp.index])


#scenario 1 - update the precipitation only and save DHSVM input files
for i in range(xcord.shape[0]):
    node_prec = df_TimePrec[df_TimePrec.columns[i]]
    #names of the columns need to be the same in both dataframes to be able to use the pandas update fcn
    names_node_prec = ['date_time', 'Precip']
    node_prec.columns = names_node_prec
    df.update(node_prec, join='left', overwrite=True)
    fileName = 'DHSVM_prec_Input_NARR_Morr_' + columnNames[i+1]
    df.to_csv(path + fileName, sep='\t')

#scenario2 - update all variables and save DHSVM input files
for i in range(xcord.shape[0]):
    #read all the WRF locations sequentially
    node_prec = df_TimePrec[df_TimePrec.columns[i]]
    node_temp = df_TimeTemp[df_TimeTemp.columns[i]]
    node_wind = df_TimeWind[df_TimeWind.columns[i]]
    node_rh = df_TimeRH[df_TimeRH.columns[i]]
    node_sw = df_TimeSW[df_TimeSW.columns[i]]
    node_lw = df_TimeLW[df_TimeLW.columns[i]]
    names_node_prec = ['date_time', 'Precip']
    names_node_temp = ['date_time', 'temp2m']
    names_node_wind = ['date_time', 'wind2m']
    names_node_rh = ['date_time', 'RH']
    names_node_sw = ['date_time', 'SW']
    names_node_lw = ['date_time', 'LW']
    node_prec.columns = names_node_prec
    node_temp.columns = names_node_temp
    node_wind.columns = names_node_wind
    node_rh.columns = names_node_rh
    node_sw.columns = names_node_sw
    node_lw.columns = names_node_lw
    df.update(node_prec, join='left', overwrite=True)
    df.update(node_temp, join='left', overwrite=True)
    df.update(node_wind, join='left', overwrite=True)
    df.update(node_rh, join='left', overwrite=True)
    df.update(node_sw, join='left', overwrite=True)
    df.update(node_lw, join='left', overwrite=True)
    fileName = 'DHSVM_all_var_Input_NARR_Morr_' + columnNames[i+1]
    df.to_csv(path + fileName, sep='\t')

# test the output

