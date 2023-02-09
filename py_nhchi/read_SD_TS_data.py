from random import lognormvariate
import numpy as np
import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
import datetime
import os
import cftime

def load_SD_qcTS_arctic2019(year, SD_name):
    ### load saildrone Arctic 2019 T, S (with quick QC) trajectory (in netcdf) data files to np array
    ### units of u, v in cm/s
    ### time format is datetime
    ### Outputs: yday of the year specified, longitude, latitude, T, S
    platf_name = ['1033','1034','1035','1036','1037','1041']
    pathSD = os.path.expanduser('~/Documents/Data/Arctic'+year+'/nhchi/')
    filenames = ['SD'+platf_name[i]+'_TS_SBE_flag_TstdGT0.1_SstdGT0.05.txt' \
        for i in range(len(platf_name))]
    print(filenames)
    ### find the index of SD_name
    for i in range(len(platf_name)):
        if platf_name[i] == SD_name:
            sd = i
            break
    ### load SD data
    filename = filenames[sd]
    dataIN = np.genfromtxt(pathSD+filename,delimiter=' ',skip_header=True,dtype=float)
    yday = dataIN[:,0]
    longitude = dataIN[:,1]
    latitude = dataIN[:,2]
    T = dataIN[:,3]
    S = dataIN[:,4]
    ### remove 1e+34
    # largenum = -1e+10
    # yday[yday< largenum] = np.nan
    # longitude[longitude< largenum] = np.nan
    # latitude[latitude< largenum] = np.nan
    # T[T< largenum] = np.nan
    # S[S< largenum] = np.nan
    # put masked array to list, then to array
    print('Successfully load',SD_name)
    print('shape of lon, lat: ',longitude.shape,latitude.shape,'shape of T: ',T.shape)

    return yday, longitude, latitude, T, S

def load_SD_qcTS_arctic2018(year, SD_name):
    ### load saildrone Arctic 2018 T, S (with quick QC) trajectory (in netcdf) data files to np array
    ### units of u, v in cm/s
    ### time format is datetime
    ### Outputs: yday of the year specified, longitude, latitude, T, S
    platf_name = ['1020','1021','1022','1023']
    pathSD = os.path.expanduser('~/Documents/Data/Arctic'+year+'/nhchi/')
    filenames = ['SD'+platf_name[i]+'_TS_SBE_flag_TstdGT0.1_SstdGT0.05.txt' \
        for i in range(len(platf_name))]
    print(filenames)
    ### find the index of SD_name
    for i in range(len(platf_name)):
        if platf_name[i] == SD_name:
            sd = i
            break
    ### load SD data
    filename = filenames[sd]
    dataIN = np.genfromtxt(pathSD+filename,delimiter=' ',skip_header=True,dtype=float)
    yday = dataIN[:,0]
    longitude = dataIN[:,1]
    latitude = dataIN[:,2]
    T = dataIN[:,3]
    S = dataIN[:,4]
    # put masked array to list, then to array
    print('Successfully load',SD_name)
    print('shape of lon, lat: ',longitude.shape,latitude.shape,'shape of T: ',T.shape)

    return yday, longitude, latitude, T, S