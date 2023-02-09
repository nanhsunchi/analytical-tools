import numpy as np
# import matplotlib
from matplotlib import cm

def values_to_rgba(values,ColormapName,vmin,vmax):
    ### Getting/ converting individual colors from a color map in matplotlib
    ### The input "values" has to be eitehr a value, list/ array. ColormapName is string. vmin & vmax are min & max to normalize colors for values.
    ### The output is rgba collapsed into 1D (use np.ndarray.flatten)

    norm = cm.colors.Normalize(vmin=vmin,vmax=vmax,clip=True)
    cmap = cm.get_cmap(ColormapName)

    ### first, convert into np.array
    if isinstance(values, (float, int)):
        rgba = np.array(cmap(norm(values)))
        return rgba, norm, cmap
    elif isinstance(values, list) | isinstance(values, np.ndarray):
        values = np.ndarray.flatten( np.array(values) )
        if values.size == 1:
            rgba = np.array(cmap(norm(values[0])))
            return rgba, norm, cmap
        else:
            rgba = np.array([None]*len(values))
            for i in range( len(values) ):
                rgba[i] = cmap(norm(values[i]))
            return rgba, norm, cmap
    else: 
        print('function only deals with single value or 1D array or list.')

# test = [0,1,2,3]
# out = values_to_rgba(test,'gray',0,2)
# print(out, type(out))