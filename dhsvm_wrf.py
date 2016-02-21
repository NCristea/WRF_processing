import pandas as pd
import xarray
import numpy as np
from dask.diagnostics import ProgressBar

path = '/Users/carina/desktop/WRF_data/'
#current WRF runs - combination of microphysics and boundary conditions
#NARR_Morr
#Narr_WSM6
#Narr_Thomm

#open the reduced dataset
dsTotal = xarray.open_dataset('/Users/carina/desktop/WRF_data/ds_reduced_NARR_Morr.nc', engine = 'scipy')
#dsTotal.chunk({'time':400,'x':50,'y':50})
xcord = x[:][0]
ycord = y[:][0]

#initialize numpy arrays
all_Temp = np.zeros((dsTotal.time.shape[0], xcord.shape[0]))
all_Prec = np.zeros((dsTotal.time.shape[0], xcord.shape[0]))
all_Wind = np.zeros((dsTotal.time.shape[0], xcord.shape[0]))
all_RH = np.zeros((dsTotal.time.shape[0], xcord.shape[0]))
all_SW = np.zeros((dsTotal.time.shape[0], xcord.shape[0]))
all_LW = np.zeros((dsTotal.time.shape[0], xcord.shape[0]))

columnNames = ['date_time']

for i in range(xcord.shape[0]):
    selectTemp = dsTotal.sel_points(x = [xcord[i]], y = [ycord[i]]).temp2m
    selectPrec = dsTotal.sel_points(x = [xcord[i]], y = [ycord[i]]).prec
    selectWind = dsTotal.sel_points(x = [xcord[i]], y = [ycord[i]]).wind2m
    selectRH = dsTotal.sel_points(x = [xcord[i]], y = [ycord[i]]).rh2m
    selectSW = dsTotal.sel_points(x = [xcord[i]], y = [ycord[i]]).SWdown
    selectLW = dsTotal.sel_points(x = [xcord[i]], y = [ycord[i]]).LWdown
    # append all WRF nodes in a dataframe for each variable
    all_Temp[:, i] = selectTemp
    all_Prec[:, i] = selectPrec
    all_Wind[:, i] = selectWind
    all_RH[:, i] = all_RH
    all_SW[:, i] = all_SW
    all_LW[:, i] = all_LW
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
df_TimeTemp.set_index('date_time', inplace = True)
df_TimePrec.set_index('date_time', inplace = True)
df_TimeWind.set_index('date_time', inplace = True)
df_TimeRH.set_index('date_time', inplace = True)
df_TimeSW.set_index('date_time', inplace = True)
df_TimeLW.set_index('date_time', inplace = True)

#open the station data file and stich WRF data
df = pd.read_csv(path + 'DHSVM_example.txt', sep = '\t')
names = ['date_time', 'temp2m', 'wind2m', 'RH', 'SW', 'LW', 'Precip']
df.columns = names
df['date_time'] = pd.to_datetime(df['date_time'], format = "%m/%d/%Y-%H")
df.set_index('date_time', inplace = True)

#scenario 1 - update the precipitation only and save DHSVM input files
for i in range(xcord.shape[0]):
    node_prec = df_TimePrec[df_TimePrec.columns[i]]
    #names of the columns need to be the same in both dataframes to be able to use the pandas update fcn
    names_node_prec = ['date_time', 'Precip']
    node_prec.columns = names_node_prec
    df.update(node_prec, join='left', overwrite=True)
    fileName = 'DHSVM_prec_Input_NARR_Morr' + columnNames[i+1]
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
    fileName = 'DHSVM_all_var_Input_NARR_Morr' + columnNames[i+1]
    df.to_csv(path + fileName, sep='\t')