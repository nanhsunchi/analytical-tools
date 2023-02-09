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