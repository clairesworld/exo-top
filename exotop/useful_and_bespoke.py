import numpy as np
def colorize(vector,cmap='plasma', vmin=None, vmax=None):
    """Convert a vector to RGBA colors. @author: jlustigy

    Parameters
    ----------
    vector : array
        Array of values to be represented by relative colors     
    cmap : str (optional)
        Matplotlib Colormap name
    vmin : float (optional)
        Minimum value for color normalization. Defaults to np.min(vector)
    vmax : float (optional)
        Maximum value for color normalization. Defaults to np.max(vector)
        
    Returns
    -------
    vcolors : np.ndarray
        Array of RGBA colors
    scalarmap : matplotlib.cm.ScalarMappable
        ScalerMap to convert values to colors
    cNorm : matplotlib.colors.Normalize
        Color normalization
    """
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import matplotlib.cm as cmx

    if vmin is None: vmin = np.min(vector)
    if vmax is None: vmax = np.max(vector)    
    
    cm = plt.get_cmap(cmap)
    cNorm  = mcolors.Normalize(vmin=vmin, vmax=vmax)
    scalarmap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    vcolors = scalarmap.to_rgba(vector)
    
    return vcolors, scalarmap, cNorm


from matplotlib.colors import LinearSegmentedColormap
def cmap_from_list(clist, n_bin=None, cmap_name=''):
    if n_bin is None:
        n_bin = len(clist)
    cm = LinearSegmentedColormap.from_list(cmap_name, clist, N=n_bin)
    return cm


from collections import Iterable
from six import string_types
def iterable_not_string(obj):
    if isinstance(obj, Iterable) and not isinstance(obj, string_types):
        return True
    else:
        return False


def printe(name, all=False):
    if all:
        print(name, '=', repr(eval(name)))
    print(name, np.shape(eval(name)))
    try:
        print(name, '[0]', np.shape(eval(name[0])))
    except:
        pass
