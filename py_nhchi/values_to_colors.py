import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib

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

def categorical_cmap(nc, nsc, cmap="tab10", continuous=False):
    '''You may use the HSV system to obtain differently saturated and luminated colors for the same hue. 
    Suppose you have at most 10 categories, then the tab10 map can be used to get a certain number of base colors. 
    From those you can choose a couple of lighter shades for the subcategories.
    
    The following would be a function categorical_cmap, 
    which takes as input the number of categories (nc) and the number of subcategories (nsc) and 
    returns a colormap with nc*nsc different colors, where for each category there are nsc colors of same hue.
    From https://stackoverflow.com/questions/47222585/matplotlib-generic-colormap-from-tab10
    '''
    if nc > plt.get_cmap(cmap).N:
        raise ValueError("Too many categories for colormap.")
    if continuous:
        ccolors = plt.get_cmap(cmap)(np.linspace(0,1,nc))
    else:
        ccolors = plt.get_cmap(cmap)(np.arange(nc, dtype=int))
    cols = np.zeros((nc*nsc, 3))
    for i, c in enumerate(ccolors):
        chsv = matplotlib.colors.rgb_to_hsv(c[:3])
        arhsv = np.tile(chsv,nsc).reshape(nsc,3)
        arhsv[:,1] = np.linspace(chsv[1],0.25,nsc)
        arhsv[:,2] = np.linspace(chsv[2],1,nsc)
        rgb = matplotlib.colors.hsv_to_rgb(arhsv)
        cols[i*nsc:(i+1)*nsc,:] = rgb       
    cmap = matplotlib.colors.ListedColormap(cols)
    return cmap

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    ''' Extract a subset of a colormap as a new colormap in matplotlib
    From and See example here: https://stackoverflow.com/questions/18926031/how-to-extract-a-subset-of-a-colormap-as-a-new-colormap-in-matplotlib
    cmap = plt.get_cmap('jet')
    new_cmap = truncate_colormap(cmap, 0.2, 0.8)
    '''
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap