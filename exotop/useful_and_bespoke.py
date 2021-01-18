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
        return False  # also False if None

def not_iterable(obj):
    if isinstance(obj, string_types):
        return True
    elif not isinstance(obj, Iterable):
        return True
    else:
        return False

def not_string(obj):
    if isinstance(obj, string_types):
        return False
    else:
        return True

from scipy.spatial import distance
def reduced_chisq(O_y, C_y, x=None, n_fitted=1):
    if x is None:
        dist = np.var(O_y)
    else:
        # use mahalanobis distance
        V = np.cov(np.array([O_y, x]).T)
        try:
            IV = np.linalg.inv(V)
        except np.linalg.LinAlgError:
            IV = np.linalg.pinv(V)  # pseudo-inverse
        # print('inv. cov', IV)
        dist = distance.mahalanobis(O_y, x, IV)
        print('D_m', dist)
        print('var(x)', np.var(x))
        print('var(O_y)', np.var(O_y))
        dist - dist**2
    dof = len(O_y) - n_fitted
    chisq = np.sum((np.array(O_y) - np.array(C_y)) / np.array(dist))
    return chisq / dof


def printe(name, obj, showall=False):
    if showall:
        print(name, '=', repr(obj))
    print(name, np.shape(obj))
    try:
        print(name, '[0]', np.shape(obj[0]))
    except:
        pass


from mpl_toolkits.axes_grid1 import make_axes_locatable
def colourbar(mappable, vmin=None, vmax=None, label='', labelsize=16, ticks=None, ticklabels=None, labelpad=17,
              rot=None, discrete=False):
    # from https://joseph-long.com/writing/colorbars/
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    cbar.set_label(label, rotation=270, labelpad=labelpad, fontsize=labelsize)
    if ticks is not None:
        cbar.set_ticks(ticks)
    elif ticks is None and discrete:
        nlabels = len(ticklabels)
        tick_locs = (np.arange(vmin, vmax + 1) + 0.5) * (nlabels - 1) / nlabels
        cbar.set_ticks(tick_locs)
    if ticklabels is not None:
        cbar.ax.set_yticklabels(ticklabels, rotation=rot)
    return cbar


def age_index(times, age, age_scale=1):
    # get index of age in times with optional scaling for age
    return min(enumerate(times), key=lambda x: abs(age - x[1] * age_scale))[0]


def minmaxnorm(x, a=0, b=1):
    # linear normalisation to min, max
    xmin = np.min(x)
    xmax = np.max(x)
    x = (x - xmin) / (xmax - xmin)  # norm to 0, 1
    x = x * (b - a) + a  # scale to a, b
    return x