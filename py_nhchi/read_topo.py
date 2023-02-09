import numpy as np
import netCDF4 as nc
import os

def load_topo_arctic():
    pathB = os.path.expanduser('~/Documents/Data/Topography/')
    filename = 'usgsCeSrtm30v6_SDArctic2019.nc'
    ds = nc.Dataset(pathB+filename)
    latTOPO = ds.variables['latitude'][:]
    lonTOPO = ds.variables['longitude'][:]
    lonTOPO = (360+(lonTOPO%360))%360
    sort_index = np.argsort(lonTOPO)
    lonTOPO = np.sort(lonTOPO)
    topo = ds.variables['topo'][:]
    # sort
    topo = topo[:,sort_index]

    return latTOPO, lonTOPO, topo

def load_topo_arctic1():
    pathB = os.path.expanduser('~/Documents/Data/Topography/')
    filename = 'usgsCeSrtm30v6_Arctic.nc'
    ds = nc.Dataset(pathB+filename)
    latTOPO = ds.variables['latitude'][:]
    lonTOPO = ds.variables['longitude'][:]
    lonTOPO = (360+(lonTOPO%360))%360
    sort_index = np.argsort(lonTOPO)
    lonTOPO = np.sort(lonTOPO)
    topo = ds.variables['topo'][:]
    # sort
    topo = topo[:,sort_index]

    return latTOPO, lonTOPO, topo

def load_topo_arctic_HiRes(lonmin,lonmax,latmin,latmax):
    pathB = os.path.expanduser('~/Documents/Data/Topography/')
    filename = 'usgsCeSrtm30v6_90N_50N.nc'
    ds = nc.Dataset(pathB+filename)
    latTOPO = ds.variables['latitude'][:]
    lonTOPO = ds.variables['longitude'][:]
    lonTOPO = (360+(lonTOPO%360))%360
    sort_index = np.argsort(lonTOPO)
    lonTOPO = np.sort(lonTOPO)
    topo = ds.variables['topo'][:]
    # sort
    topo = topo[:,sort_index]
    # slice the data specified
    lonmin = (360+(lonmin%360))%360
    lonmax = (360+(lonmax%360))%360
    ilon = np.where( (lonTOPO>= lonmin) & (lonTOPO<= lonmax) )[0]
    ilat = np.where( (latTOPO>= latmin) & (latTOPO<= latmax) )[0]
    latTOPO = latTOPO[ilat]
    lonTOPO = lonTOPO[ilon]
    topo = topo[np.ix_(ilat,ilon)]

    return latTOPO, lonTOPO, topo