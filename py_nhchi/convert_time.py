import numpy as np
import datetime

def datetime_to_yearday(datetimeIN):
    ### Convert 1 datetime value or 1D array/ list of year datetime to yearday (include fractional day) of that year
    ### Output year, yearday in nd-array

    ### first, convert into np.array
    if isinstance(datetimeIN, datetime.date):
        year = int( datetimeIN.timetuple().tm_year )
        yday_int = datetimeIN.timetuple().tm_yday
        day_frac = (datetimeIN-datetime.datetime.strptime(str(year)+'-'+str(yday_int),'%Y-%j')).total_seconds()/86400
        yearday = yday_int + day_frac
        return np.array(year), np.array(yearday)
    elif isinstance(datetimeIN, np.ndarray) | isinstance(datetimeIN, list):
        ### conver to 1d list
        datetimeIN = np.squeeze( np.ndarray.flatten(np.array(datetimeIN)) )
        datetimeIN = datetimeIN.tolist()
        if isinstance(datetimeIN, datetime.date):
            year = int( datetimeIN.timetuple().tm_year )
            yday_int = datetimeIN.timetuple().tm_yday
            day_frac = (datetimeIN-datetime.datetime.strptime(str(year)+'-'+str(yday_int),'%Y-%j')).total_seconds()/86400
            yearday = yday_int + day_frac
            return np.array(year), np.array(yearday)
        elif len(datetimeIN)> 1:
            yearday = []
            year = []
            for i in range( len(datetimeIN) ):
                year.append( int( datetimeIN[i].timetuple().tm_year ) )
                yday_int = datetimeIN[i].timetuple().tm_yday
                day_frac = (datetimeIN[i]-datetime.datetime.strptime(str(year[i])+'-'+str(yday_int),'%Y-%j')).total_seconds()/86400
                yearday.append(yday_int + day_frac)
            return np.array(year), np.array(yearday)
    else: 
        print('datetime_to_yearday function only deals with 1D array or list.')
        return None

# test = datetime.datetime(2019,1,15,12,0,5)
# test = [datetime.datetime(2019,1,15,12,0,5)]
# test = np.array([datetime.datetime(2019,1,15,12,0,5)])
# test = np.array([datetime.datetime(2019,1,15,12,0,5),datetime.datetime(2019,1,16)])
# print(type(test))
# print(datetime_to_yearday(test))

def jday_to_datetime(jday, ref_dtime):
    ### convert jday (days after a specified reference datetime, ref_dtime.). This can be either a number or a 1D np.array/ list. 
    ### output: datetime
    if isinstance(ref_dtime, datetime.date):
        # print('input ref_dtime ok')
        pass
    else:
        print('input reference datetime needs to be a datetime.date object')
        return 'input reference datetime needs to be a datetime.date object'
    if isinstance(jday, float) or isinstance(jday, int) or isinstance(jday, np.int64):
        sec_after_ref = float(jday)*86400 # in sec
        dtime_out = ref_dtime+datetime.timedelta(seconds=sec_after_ref)
        return dtime_out
    else:
        dtime_out = [np.nan]*len(jday)
        for i in range( len(jday) ):
            sec_after_ref = jday[i]*86400
            dtime_out[i] = ref_dtime+datetime.timedelta(seconds=sec_after_ref)
        return dtime_out
        # day_frac = (jday-datetime.datetime.strptime(str(year)+'-'+str(yday_int),'%Y-%j')).total_seconds()/86400
        # yearday = yday_int + day_frac
        # return np.array(year), np.array(yearday)

# print(1)
# print(jday_to_datetime([1], datetime.datetime(2018,12,31)))
# print(2)
# print(jday_to_datetime(250000.5, 100))
# print(3)
# print( jday_to_datetime([250000.5,2500.2], datetime.datetime(1950,1,1)) )
