import numpy as np
import math

def distance(lat1,lon1,lat2,lon2):
    # calculate Distance Between Two Points on Earth
    # input: lat1,lon1,lat2,lon2
    # output: distance in kilometer (km)
    # The math module contains a function named
    # radians which converts from degrees to radians.
    lon1 = math.radians(lon1)
    lon2 = math.radians(lon2)
    lat1 = math.radians(lat1)
    lat2 = math.radians(lat2)

    # Haversine formula
    dlon = lon2 - lon1 
    dlat = lat2 - lat1
    a = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.acos
    c = 2 * math.asin(math.sqrt(a)) 
    # Radius of earth in kilometers. Use 3956 for miles
    r = 6371
     # calculate the result
    return(c * r)