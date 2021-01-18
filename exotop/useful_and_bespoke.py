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

import pandas as pd

def mahalanobis(x=None, data=None, cov=None):
    """Compute the Mahalanobis Distance between each row of x and the data
    x    : vector or matrix of observed data with, say, p columns.
    data : ndarray of the distribution from which Mahalanobis distance of each observation of x is to be computed.
    cov  : covariance matrix (p x p) of the distribution. If None, will be computed from data but must be df.
    """
    x_minus_mu = x - np.mean(data)
    if not cov:
        cov = np.cov(data.values.T)
    try:
        inv_covmat = np.linalg.inv(cov)
    except np.linalg.LinAlgError:
        inv_covmat = np.linalg.pinv(cov)  # pseudo-inverse
    left_term = np.dot(x_minus_mu, inv_covmat)
    mahal = np.dot(left_term, x_minus_mu.T)
    x['mahala^2'] = mahal.diagonal()**2
    print(x.head(200))
    return mahal.diagonal()


from scipy.spatial import distance
from scipy.stats import chisquare

def reduced_chisq(O_y, C_y, dist=None, n_fitted=2, **kwargs):
    # dist is an array of distance metrics e.g. variance or Mahalanobis for each point in O_y
    print('O_y', O_y)
    print('C_y', C_y)
    print('dist', dist)
    if dist is None:  # default to simple variance if not provided
        dist = np.var(O_y)
    print('var(O_y)', np.var(O_y))
    dof = len(O_y) - n_fitted
    print('dof', dof)
    chisq = np.sum((np.array(O_y) - np.array(C_y))**2 / np.array(dist))
    # chi2 = np.sum(((array(X) - array(X_model)) ** 2 + (array(Y) - array(Y_model)) ** 2) / (s ** 2))  # 2D
    print('chisquare', chisq / dof)
    print('chisquare with y variance', np.sum((np.array(O_y) - np.array(C_y))**2 / np.var(O_y)) / dof)
    print('scipy chisquare', chisquare(O_y, C_y)[0] / dof)
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