""" bunch of custom functions """

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


def cmap_from_list(clist, n_bin=None, cmap_name=''):
    from matplotlib.colors import LinearSegmentedColormap
    if n_bin is None:
        n_bin = len(clist)
    cm = LinearSegmentedColormap.from_list(cmap_name, clist, N=n_bin)
    return cm


def cmap_from_ascii(name, path='', end='.txt'):
    from matplotlib.colors import ListedColormap
    file = path + name + end
    carray = np.genfromtxt(file, comments='#')
    cmap = ListedColormap(carray, name=name)
    return cmap


def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""
    # By Jake VanderPlas
    # License: BSD-style
    import matplotlib.pyplot as plt
    import numpy as np

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)


def iterable_not_string(obj):
    from collections import Iterable
    from six import string_types
    if isinstance(obj, Iterable) and not isinstance(obj, string_types):
        return True
    else:
        return False  # also False if None

def not_iterable(obj):  # but can be a string
    from collections import Iterable
    from six import string_types
    if isinstance(obj, string_types):
        return True
    elif not isinstance(obj, Iterable):
        return True
    else:
        return False

def not_string(obj):
    from six import string_types
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
    if cov is None:
        cov = np.cov(data.values.T)
    try:
        inv_covmat = np.linalg.inv(cov)
    except np.linalg.LinAlgError:
        inv_covmat = np.linalg.pinv(cov)  # pseudo-inverse
    left_term = np.dot(x_minus_mu, inv_covmat)
    mahal = np.dot(left_term, x_minus_mu.T)
    # x['mahala^2'] = mahal.diagonal()**2
    # print(x.head(200))
    return np.sqrt(mahal.diagonal())


def reduced_chisq(O_y, C_y, dist=None, n_fitted=2, **kwargs):
    # from scipy.spatial import distance
    # from scipy.stats import chisquare
    # dist is an array of distance metrics / errors e.g. variance or Mahalanobis for each point in O_y
    dof = len(O_y) - n_fitted
    chisq = np.sum((np.array(O_y) - np.array(C_y))**2 / np.array(dist))
    # chi2 = np.sum(((array(X) - array(X_model)) ** 2 + (array(Y) - array(Y_model)) ** 2) / (s ** 2))  # 2D
    # print('chisquare', chisq / dof)
    return chisq / dof


def printe(name, obj, showall=False):
    if showall:
        print(name, '=', repr(obj))
    print(name, np.shape(obj))
    try:
        print(name, '[0]', np.shape(obj[0]))
    except:
        pass


def colourbar(mappable=None, vector=None, ax=None, vmin=None, vmax=None, label='', labelsize=16, ticksize=14,
              ticks=None, ticklabels=None, labelpad=17,
              rot=None, discrete=False, cmap='rainbow', tickformatter=None, c='k', pad=0.05, log=False):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.colors as colors
    from matplotlib.pyplot import clim
    # from https://joseph-long.com/writing/colorbars/

    if ax is None:
        ax = mappable.axes
    if log:
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = None
    if mappable is None:
        try:
            n = len(vector)
            if vmin is None:
                vmin = np.min(vector)
            if vmax is None:
                vmax = np.max(vector)
        except TypeError as e:
            print(e)
            raise Exception('colourbar: if mappable is None, must provide vector')
        dum = np.linspace(vmin, vmax, n)
        mappable = ax.scatter(dum, dum, c=dum, cmap=cmap, s=0, norm=norm)

    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=pad)
    cbar = fig.colorbar(mappable, cax=cax)
    cbar.set_label(label, rotation=270, labelpad=labelpad, fontsize=labelsize, c=c)

    if ticks is not None:
        cbar.set_ticks(ticks)
    elif ticks is None and discrete:
        nlabels = len(ticklabels)
        tick_locs = (np.arange(vmin, vmax + 1) + 0.5) * (nlabels - 1) / nlabels
        cbar.set_ticks(tick_locs)
    if ticklabels is not None:
        cbar.ax.set_yticklabels(ticklabels, rotation=rot)
    # if tickformatter is not None:
    #     cbar.ax.yaxis.set_major_formatter(tickformatter)

    cbar.ax.yaxis.label.set_color(c)
    cbar.ax.tick_params(axis='y', colors=c, labelsize=ticksize)
    [cbar.ax.spines[s].set_color(c) for s in ['bottom', 'top', 'right', 'left']]
    return cbar


def colourised_legend(ax, clist, cleglabels, lw=0, ls='--', marker='o', markersize=20, alpha=1, legsize=25, ncol=1, title=None, **kwargs):
    import matplotlib.lines as mlines
    handles = []
    for jj, label in enumerate(cleglabels):
        handles.append(mlines.Line2D([], [], color=clist[jj], marker=marker, ls=ls, alpha=alpha,
                                     markersize=markersize, lw=lw, label=str(label)))
    leg = ax.legend(handles=handles, frameon=False, fontsize=legsize, ncol=ncol, bbox_to_anchor=(1.01, 1), loc='upper left', title=title, **kwargs)
    if title is not None:
        leg.get_title().set_fontsize(legsize)  # legend 'Title' fontsize
    ax.add_artist(leg)
    return ax


def age_index(times, age, age_scale=1):
    # get index of age in times with optional scaling for age
    return min(enumerate(times), key=lambda x: abs(age - x[1] * age_scale))[0]


def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return int(idx)

def unique_rows(a):
    # Given a numpy array, return another numpy array with only the unique rows
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

def unique_rows_indices(a):
    # Given a numpy array, return another numpy array with only the unique rows
    a = np.ascontiguousarray(a)
    unique_a, indices = np.unique(a.view([('', a.dtype)]*a.shape[1]), return_index=True)
    return indices
    #return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

def minmaxnorm(x, a=0, b=1):
    # linear normalisation to min, max
    x = np.array(x)
    xmin = np.min(x)
    xmax = np.max(x)
    x = (x - xmin) / (xmax - xmin)  # norm to 0, 1
    x = x * (b - a) + a  # scale to a, b
    return x


# image scatter fn
def imscatter(x, y, image, ax=None, zoom=1):
    # from PIL import Image
    import matplotlib.pyplot as plt
    from matplotlib.offsetbox import OffsetImage, AnnotationBbox
    if ax is None:
        ax = plt.gca()
    try:
        image = plt.imread(image)
    except TypeError:
        # Likely already an array...
        pass
    im = OffsetImage(image, zoom=zoom)
    x, y = np.atleast_1d(x, y)
    artists = []
    for x0, y0 in zip(x, y):
        ab = AnnotationBbox(im, (x0, y0), xycoords='data', frameon=False)
        artists.append(ax.add_artist(ab))
    ax.update_datalim(np.column_stack([x, y]))
    ax.autoscale()
    return artists


def hide_log_ticklabels(ax, axis='both', index='all', hideticks=False, flipped=False):
    # hilariously tricky thing of hiding tick labels for log scale. answer from
    # https://stackoverflow.com/questions/36064477/remove-specific-ticks-on-logarithmic-plot-in-matplotlib
    from matplotlib.ticker import NullFormatter
    if not (axis in ['x', 'y', 'both']):
        raise Exception('axis must be x, y, or both')
    if not (index in ['first', 'last', 'all']):
        raise Exception('index must be first, last, or all')
    axes = []
    if axis in ['x', 'both']:
        axes.append(ax.xaxis)
    if axis == ['y', 'both']:
        axes.append(ax.yaxis)

    for a in axes:
        a.set_minor_formatter(NullFormatter())
        if (index == 'first') or (index == 'last' and flipped):
            tcks = [a.get_major_ticks()[1]]
        elif (index == 'last') or (index == 'first' and flipped):
            tcks = [a.get_major_ticks()[-2]]
        elif index == 'both':
            tcks = a.get_major_ticks()
        for t in tcks:
            t.label1.set_visible(False)
        if hideticks:
            [a.set_tick_params(which=m, size=0) for m in ['minor', 'major']]
            [a.set_tick_params(which=m, width=0) for m in ['minor', 'major']]


def dark_background(fig, ax, fgc='xkcd:off white', bgc='xkcd:black'):
    """ recolour fig and its axes to foreground and background colour - not artists tho """
    import matplotlib.pyplot as plt
    from matplotlib.legend import Legend
    if not_iterable(ax):
        ax = [ax]

    fig.patch.set_facecolor(bgc)

    for a in ax:
        a.set_facecolor(bgc)
        [a.spines[s].set_color(fgc) for s in ['bottom', 'top', 'right', 'left']]
        a.xaxis.label.set_color(fgc)
        a.tick_params(axis='x', colors=fgc)
        a.yaxis.label.set_color(fgc)
        a.tick_params(axis='y', colors=fgc)

        # if there's a legend do that too
        legends = [c for c in a.get_children() if isinstance(c, Legend)]
        for l in legends:
            for text in l.get_texts():
                plt.setp(text, color=fgc)
            # todo: spines? (legend.edgecolor)
        # todo: title color
    print('Remember to add facecolor=fig.get_facecolor() to savefig()')

    return (fig, *ax)


def cornertext(ax, text, pos='top right', size=12, pad=0.05, **kwargs):
    if 'top' in pos:
        y = 1 - pad
        va = 'top'
    elif 'bottom' in pos:
        y = pad
        va = 'bottom'
    if 'left' in pos:
        x = pad
        ha = 'left'
    elif 'right' in pos:
        x = 1 - pad
        ha = 'right'

    # update?
    pass_args = {'x':x, 'y':y, 'va':va, 'ha':ha}
    pass_args.update(kwargs)
    x = pass_args.pop('x')
    y = pass_args.pop('y')

    ax.text(x, y, text, transform=ax.transAxes, fontsize=size, **pass_args)
    return ax

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]


def get_continuous_cmap(hex_list, float_list=None, N=256):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list.

        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.

        Returns
        ----------
        colour map'''
    import matplotlib.colors as mcolors
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list is not None:
        pass
    else:
        float_list = list(np.linspace(0, 1, len(rgb_list)))

    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=N)
    return cmp
