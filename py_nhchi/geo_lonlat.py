import numpy as np
import math

def lonlat_to_perpendicular_vector(lon1,lat1,lon2,lat2):
    # follow these steps:
    # https://gis.stackexchange.com/questions/29664/given-a-line-on-the-earths-surface-how-do-i-plot-a-line-perpendicular-to-it
    # step 1
    lon1_adj = lon1*math.cos(math.radians((lat1+lat2)/2))
    lon2_adj = lon2*math.cos(math.radians((lat1+lat2)/2))
    # print(lon1_adj, lat1, lon2_adj, lat2)
    # print(lon2_adj-lon1_adj, lat2-lat1)
    # step 2
    # The direction vector for the path is the displacement from its beginning to its end
    v = np.array([lon2_adj,lat2])-np.array([lon1_adj,lat1])
    # Turning any vector (x,y) at right angles clockwise produces (y,-x), whence a perpendicular direction to the right is
    w = np.array([v[1],-v[0]])
    return v,w
    
def vectors_to_angle(a,b):
    # input: a & b have to be np.array
    a = np.array(a)
    b = np.array(b)
    # compute the angle between two vectors in numpy
    # compute angle between two vectors in numpy
    inner = np.inner(a,b)
    norms = np.linalg.norm(a) * np.linalg.norm(b)
    cos = inner / norms
    rad = np.arccos(np.clip(cos,-1,1))
    degree = np.rad2deg(rad)
    return rad, degree, cos