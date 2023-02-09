import numpy as np
import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
import datetime
import os
import cftime

def load_SD_nc_arctic2019(year, SD_name):
    ### load saildrone Arctic 2019 ADCP trajectory (in netcdf) data files to masked np arrays
    ### units of u, v in cm/s
    ### time format is datetime
    ### Outputs: time, longitude, latitude, u, v
    platf_name = ['1035','1036','1037']
    pathSD = os.path.expanduser('~/Documents/Data/Arctic'+year+'/')
    pathSD_name = [pathSD,pathSD+'nhchi/',pathSD+'nhchi/']
    fext_name = ['.cdf','.nc','.nc']
    time_name = ['MYNEWT','time','time']
    lon_name = ['LONADCP','longitude','longitude']
    lat_name = ['LATADCP','latitude','latitude']
    u_name = ['U','u','u']
    v_name = ['V','v','v']
    depth_name = ['DEPTH_ADCP','cell_depth','cell_depth']

    ### find the index of SD_name
    for i in range(len(platf_name)):
        if platf_name[i] == SD_name:
            sd = i
            break
    ### load SD data
    filename = 'SD'+platf_name[sd]+'uvadcp'+year+fext_name[sd]
    ds = nc.Dataset(pathSD_name[sd]+filename)
    time = ds.variables[time_name[sd]] # nc class variable
    dates = cftime.num2date(time,units=time.units,calendar=time.calendar)
    # put masked array to list, then to array
    dtime_list = np.ma.masked_array.tolist(dates)
    dtime_arr = np.asarray(dtime_list)
    dtime_arr = np.reshape(dtime_arr,(dtime_arr.size,))
    for i in range(0, dtime_arr.size): # convert to datetime
        dtime_arr[i] = datetime.datetime.strptime(str(dtime_arr[i]),'%Y-%m-%d %H:%M:%S')
    longitude = ds.variables[lon_name[sd]][:]
    latitude = ds.variables[lat_name[sd]][:]
    depth = ds.variables[depth_name[sd]][:]
    u = np.squeeze(ds.variables[u_name[sd]][:])
    v = np.squeeze(ds.variables[v_name[sd]][:])
    if sd> 0:
        u = u*100 # convert 1036 & 1037 to cm/s
        v = v*100
    print('Successfully load',SD_name)
    print('shape of lon, lat, u: ',longitude.shape,latitude.shape,'shape of u: ',u.shape)

    return dtime_arr, longitude, latitude, depth, u, v


def load_tide_current_arctic(year,SD_name):
    # load the respective saildrone barotropic tidal data (have the same time dimension as their nc data file)
    # estimated by AVISO FES2014 data. 

    # load saildrone tidal current data
    pathTide = os.path.expanduser('~/Documents/Data/Arctic'+year+'/nhchi/')
    dataTide = np.genfromtxt(pathTide+'SD'+SD_name+'uvadcp'+year+'_tidal.txt',delimiter=' ',skip_header=False,dtype=float)
    u_tide = dataTide[:,3]
    v_tide = dataTide[:,4]
    print('Succesfully load tidal current of',SD_name)
    return u_tide, v_tide


def detide_current_arctic(u,v,u_tide,v_tide):
    # perform u[index_trajectory, index_depth] - u[index_trajectory] for each depth
    # --> remove barotropic tidal current
    if u.shape[0] == u_tide.size:
        u_detide = np.nan*np.ma.ones(u.shape)
        v_detide = np.nan*np.ma.ones(u.shape)
        for i in range(u.shape[1]):
            u_detide[:,i] = u[:,i]-u_tide[:]
            v_detide[:,i] = v[:,i]-v_tide[:]
        print('Successfully detide')
        return u_detide, v_detide
    else:
        print('The length of the trajectory velocity data and its tidal current do not match!')


def load_SD_nc_arctic2018(year, SD_name):
    ### load saildrone Arctic 2018 ADCP trajectory (in netcdf) data files to masked np arrays
    ### units of u, v in cm/s
    ### time format is datetime
    ### Outputs: time, longitude, latitude, u, v
    platf_name = ['1020','1021']
    pathSD = os.path.expanduser('~/Documents/Data/Arctic'+year+'/')
    pathSD_name = [pathSD,pathSD]
    fext_name = ['.cdf','.cdf']
    time_name = ['MYNEWT1','MYNEWT1']
    lon_name = ['LON','LON']
    lat_name = ['LAT','LAT']
    u_name = ['U','U']
    v_name = ['V','V']
    depth_name = ['DEPTH_ADCP','DEPTH_ADCP']
    ### find the index of SD_name
    for i in range(len(platf_name)):
        if platf_name[i] == SD_name:
            sd = i
            break
    ### load SD data
    filename = 'ADCPsd'+platf_name[sd]+'-'+year+fext_name[sd]
    ds = nc.Dataset(pathSD_name[sd]+filename)
    time = ds.variables[time_name[sd]] # nc class variable
    dates = cftime.num2date(time,units=time.units,calendar=time.calendar)
    # put masked array to list, then to array
    dtime_list = np.ma.masked_array.tolist(dates)
    dtime_arr = np.asarray(dtime_list)
    dtime_arr = np.reshape(dtime_arr,(dtime_arr.size,))
    for i in range(0, dtime_arr.size): # convert to datetime
        dtime_arr[i] = datetime.datetime.strptime(str(dtime_arr[i]),'%Y-%m-%d %H:%M:%S')
    longitude = ds.variables[lon_name[sd]][:]
    latitude = ds.variables[lat_name[sd]][:]
    depth = ds.variables[depth_name[sd]][:]
    u = np.squeeze(ds.variables[u_name[sd]][:])
    v = np.squeeze(ds.variables[v_name[sd]][:])
    print('shape of lon, lat, u: ',longitude.shape,latitude.shape,'shape of u: ',u.shape)
    ### remove some data
    # if SD_name == '1021':
    #     iremove = (dtime_arr<= datetime.datetime(2018, 6, 30, 22, 15)) | \
    #         ((dtime_arr>= datetime.datetime(2018, 8, 21, 0, 4)) & (dtime_arr<= datetime.datetime(2018, 8, 22, 3, 25)))\
    #             | (dtime_arr>= datetime.datetime(2018, 10, 6, 22, 49))
    #     dtime_arr = np.delete(dtime_arr, iremove)
    #     longitude = np.delete(longitude, iremove)
    #     latitude = np.delete(latitude, iremove)
    #     u = np.delete(u, iremove, axis=0)
    #     v = np.delete(v, iremove, axis=0)
    ###
    print('Successfully load',SD_name)
    print('shape of lon, lat, u: ',longitude.shape,latitude.shape,'shape of u: ',u.shape)

    return dtime_arr, longitude, latitude, depth, u, v
