import numpy as np
import netCDF4 as nc
import datetime
import os

def load_nc2dic(filename, sticker, basetime, timedelta):
    ''' Read netcdf files and put data in a dictionary 
        - input:    (1) path+filename (str)
                    (2) a sticker name (str, could be empty str "") to the new variable name in the dictionary
                    (3) if there is "time" variable in the netcdf file, need 
                        - base time (datetime)
                        - timedelta (datetime.timedelta): timedelta(seconds=1), timedelta(days=1), ...
        -  output:  a dictionary with variables (all fill_values are changed to np.nan) '''
    dic = {}
    ds = nc.Dataset(filename)
    keys = list( ds.variables.keys() )
    ''' put data to dictionary '''
    for key in keys:
        data_ma = np.squeeze( ds.variables[key][:] )
        fillvalue = data_ma.get_fill_value()
        shape = data_ma.shape
        if len(shape) == 0:
            print(key,'shape is (). Will skip.')
        else:
            print(key, shape, 'fillvalue=',fillvalue)
            # np.ma.set_fill_value( data_ma, np.nan )
            # data = np.squeeze( np.ma.getdata( data_ma ) )
            # ifill = np.where( (data_ma== fillvalue) | (np.abs(data_ma>=1e20)) )[0]
            ifill = (np.ma.getmask(data_ma)) | (np.abs(data_ma>= 1e20))
            data = data_ma
            data[ifill] = np.nan
            ''' make new variable names (1) lower cases  
                                        (3) vel_east to u, (4) vel_north to v
                                        (5) lonXXX to longitude (5) latXXX to latitude '''
            var = key.lower()
            if var == 'vel_east':
                var = 'u'
            if var == 'vel_north':
                var = 'v'
            if 'lon' in var:
                var = 'longitude'
            if 'lat' in var:
                var = 'latitude'
            # i1 = var.find('_')
            # if i1> 0:
            #     var = var[:i1]
            if sticker != '':
                var = sticker+'-'+var
            # print('variable name= '+var)
            dic[var] = data
            ''' add datetime to dictionary '''
            if 'time' in key.lower():
                dtime = np.array([basetime+item*timedelta for item in dic[var]])
                dic[sticker+'-dtime'] = dtime
                print('dtime rangen:',np.min(dtime),np.max(dtime))
    print('load '+filename+' to dictionary')
    print('variables in dictionary:',list(dic.keys()))
    return dic