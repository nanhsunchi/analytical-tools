import numpy as np
import math

def vector_dir_to_from(u,v):
    ### input: u and v are vector component in zonal and meridional direction.
    ###        Take one number at a time. 
    ### output: vector (ex: wind or current) to and from (in degrees)
    ### https://stackoverflow.com/questions/21484558/how-to-calculate-wind-direction-from-u-and-v-wind-components-in-r/21486028#21486028
    
    if isinstance(u,int) or isinstance(u,float):
        ### indices of calm winds
        is_calm_wind = ((u==0) & (v==0))
        if is_calm_wind:
            ### assign u=v=0 to 0 vector direction
            vector_dir_cardinal_to = 0
            vector_dir_cardinal_from = 0
        else:
            vector_dir_trig_to = math.atan2(v, u)
            vector_dir_trig_to_degrees = vector_dir_trig_to * 180/np.pi
            ### convert vector to the meteorological convention of the direction the wind is coming from:
            vector_dir_trig_from_degrees = (vector_dir_trig_to_degrees + 180)%360
            vector_dir_cardinal_to = ((90-vector_dir_trig_from_degrees)+180)%360
            vector_dir_cardinal_from = (vector_dir_cardinal_to + 180)%360
        
        # print(vector_dir_cardinal_to,'deg', vector_dir_cardinal_from,'deg')
        return vector_dir_cardinal_to, vector_dir_cardinal_from
    else:
        print('vector_dir_to_from only takes integer or float inputs')


def rotate_origin_only(xy, radians):
    """Only rotate a point around the origin (0, 0).
    https://gist.github.com/LyleScott/e36e08bfb23b1f87af68c9051f985302
    
    Note: The input of each item in tuple "xy" should be a number or an array, not list. 
    radians is a number, not an array or list. 
    """
    x, y = xy
    xx = x * math.cos(radians) + y * math.sin(radians)
    yy = -x * math.sin(radians) + y * math.cos(radians)

    return xx, yy