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


def cmap_from_ascii(name, path='', end='.txt', **kwargs):
    from matplotlib.colors import ListedColormap
    file = path + name + end
    carray = np.genfromtxt(file, comments='#', **kwargs)
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
              ticks=None, ticklabels=None, labelpad=17, loc='right',
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
        dum = np.linspace(np.min(vector), np.max(vector), n)
        print('colourmap bounds', np.min(vector), np.max(vector))
        print('colourmax vmin', vmin, 'vmax', vmax)
        mappable = ax.scatter(dum, dum, c=dum, cmap=cmap, s=0, norm=norm, vmin=vmin, vmax=vmax)

    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(loc, size="5%", pad=pad)
    if loc == 'top':
        cbar = fig.colorbar(mappable, cax=cax, orientation='horizontal')
        cax.xaxis.set_ticks_position("top")
        cbar.ax.xaxis.set_label_position('top')
        rotation = 0
    else:
        cbar = fig.colorbar(mappable, cax=cax)
        rotation = 270
    cbar.set_label(label, rotation=rotation, labelpad=labelpad, fontsize=labelsize, c=c)

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


def colourised_legend(ax, clist, cleglabels, lw=0, ls='--', marker='o', markersize=20, alpha=1, legsize=25,
                      titlesize=None, ncol=1, title=None, return_leg=False, **kwargs):
    import matplotlib.lines as mlines
    handles = []
    for jj, label in enumerate(cleglabels):
        handles.append(mlines.Line2D([], [], color=clist[jj], marker=marker, ls=ls, alpha=alpha,
                                     markersize=markersize, lw=lw, label=str(label)))
    leg = ax.legend(handles=handles, frameon=False, fontsize=legsize, ncol=ncol, bbox_to_anchor=(1.01, 1), loc='upper left', title=title, **kwargs)
    if title is not None:
        if titlesize is None:
            titlesize = legsize
        leg.get_title().set_fontsize(titlesize)  # legend 'Title' fontsize
    ax.add_artist(leg)
    if return_leg:
        return ax, leg
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



# from colormath.color_objects import *
# from colormath.color_conversions import convert_color, color_conversion_function
# from colormath import color_diff
#
#
# class MshColor(IlluminantMixin, ColorBase):
#     '''
#     Represents an Msh color as defined by [Moreland2009]. The Msh color space
#     is basically just a polar representation of CIE Lab space.
#     See `Diverging Color Maps for Scientific Visualization
#     <http://www.kennethmoreland.com/color-maps/>` for more information.
#     '''
#
#     VALUES = ['msh_m', 'msh_s', 'msh_h']
#
#     def __init__(self, msh_m, msh_s, msh_h, observer='2', illuminant='d50'):
#         """
#         :param float msh_m: M coordinate.
#         :param float msh_s: s coordinate.
#         :param float msh_h: h coordinate.
#         :keyword str observer: Observer angle. Either ```2``` or ```10``` degrees.
#         :keyword str illuminant: See :doc:`illuminants` for valid values.
#         """
#
#         super(MshColor, self).__init__()
#         #: M coordinate
#         self.msh_m = float(msh_m)
#         #: C coordinate
#         self.msh_s = float(msh_s)
#         #: H coordinate
#         self.msh_h = float(msh_h)
#
#         #: The color's observer angle. Set with :py:meth:`set_observer`.
#         self.observer = None
#         #: The color's illuminant. Set with :py:meth:`set_illuminant`.
#         self.illuminant = None
#
#         self.set_observer(observer)
#         self.set_illuminant(illuminant)
#
#
# @color_conversion_function(LabColor, MshColor)
# def Lab_to_Msh(cobj, *args, **kwargs):
#     """
#     Convert from CIE Lab to Msh.
#     """
#
#     msh_m = math.sqrt(math.pow(float(cobj.lab_l), 2) +
#                       math.pow(float(cobj.lab_a), 2) +
#                       math.pow(float(cobj.lab_b), 2))
#     msh_s = math.acos(float(cobj.lab_l) / msh_m)
#     msh_h = math.atan2(float(cobj.lab_b), float(cobj.lab_a))
#
#     return MshColor(msh_m, msh_s, msh_h,
#                     observer=cobj.observer,
#                     illuminant=cobj.illuminant)
#
#
# @color_conversion_function(MshColor, LabColor)
# def Msh_to_Lab(cobj, *args, **kwargs):
#     """
#     Convert from Msh to Lab.
#     """
#
#     lab_l = cobj.msh_m * math.cos(float(cobj.msh_s))
#     lab_a = cobj.msh_m * math.sin(float(cobj.msh_s)) * math.cos(float(cobj.msh_h))
#     lab_b = cobj.msh_m * math.sin(float(cobj.msh_s)) * math.sin(float(cobj.msh_h))
#     return LabColor(lab_l, lab_a, lab_b,
#                     illuminant=cobj.illuminant,
#                     observer=cobj.observer)
#
#
# class SmoothDivergingColorMap:
#     def __init__(self,
#                  low_color=sRGBColor(0.230, 0.299, 0.754),
#                  high_color=sRGBColor(0.706, 0.016, 0.150),
#                  mid_color=MshColor(88.0, 0.0, 0.0)):
#         """
#         :param color low_color: The color at the low end of the map.
#         :param color high_color: The color at the high end of the map.
#         :param color mid_color: The color at the middle of the map. Should be unsaturated.
#         """
#         self.low_msh = convert_color(low_color, MshColor)
#         self.high_msh = convert_color(high_color, MshColor)
#
#         # If the points are saturated and distinct, then we place a white point
#         # in the middle. Otherwise we ignore it.
#         if self.low_msh.msh_s > 0.05:
#             if self.high_msh.msh_s > 0.05:
#                 if (abs(self.low_msh.msh_h - self.high_msh.msh_h) > math.pi / 3.0) \
#                         and mid_color:
#                     # Both endpoints are saturated and unique and a midpoint was
#                     # given. Interpolate through this midpoint and compute an
#                     # appropriate hue spin.
#                     mid_msh = convert_color(mid_color, MshColor)
#                     self.midpoint_magnitude = mid_msh.msh_m
#                     self.midpoint_low_hue = self.compute_hue_spin(self.low_msh, mid_msh)
#                     self.midpoint_high_hue = self.compute_hue_spin(self.high_msh, mid_msh)
#                 else:
#                     # Both endpoints are distinct colors, but they are either very close
#                     # in hue or no middle point was given. In this case, interpolate
#                     # directly between them.
#                     self.midpoint_magnitude = None
#             else:
#                 # The low color is saturated but the high color is unsaturated.
#                 # Interpolate directly between them, but adjust the hue of the unsaturated
#                 # high color.
#                 self.midpoint_magnitude = None
#                 self.high_msh.msh_h = self.compute_hue_spin(self.low_msh, self.high_msh)
#         else:
#             # The low color is unsaturated. Assume the high color is saturated. (If not,
#             # then this is a boring map no matter what we do.) Interpolate directly
#             # between them, but adjust the hue of the unsaturated low color.
#             self.midpoint_magnitude = None
#             self.low_msh.msh_h = self.compute_hue_spin(self.high_msh, self.low_msh)
#
#     def compute_hue_spin(self, MshSaturated, MshUnsaturated):
#         '''
#         Given a saturated color and unsaturated color, both as MshColor objects,
#         computes a spin component to use during interpolation in Msh space. The spin
#         is considered the target hue to interpolate to.
#         '''
#         if MshSaturated.msh_m >= MshUnsaturated.msh_m:
#             return MshSaturated.msh_h
#         else:
#             hSpin = (MshSaturated.msh_s *
#                      math.sqrt(math.pow(MshUnsaturated.msh_m, 2) -
#                                math.pow(MshSaturated.msh_m, 2)) /
#                      (MshSaturated.msh_m * math.sin(MshSaturated.msh_s)))
#             if MshSaturated.msh_h > -math.pi / 3:
#                 return MshSaturated.msh_h + hSpin
#             else:
#                 return MshSaturated.msh_h - hSpin
#
#     def print_self(self):
#         print('Low Color:')
#         print('\t', self.low_msh)
#         print('\t', convert_color(self.low_msh, LabColor))
#         print('\t', convert_color(self.low_msh, sRGBColor))
#
#         print('Middle Color:')
#         if (self.midpoint_magnitude):
#             print('\t Magnitude', self.midpoint_magnitude)
#             print('\t Low Hue', self.midpoint_low_hue)
#             print('\t High Hue', self.midpoint_high_hue)
#         else:
#             print('\t No Midpoint')
#
#         print('High Color:')
#         print('\t', self.high_msh)
#         print('\t', convert_color(self.high_msh, LabColor))
#         print('\t', convert_color(self.high_msh, sRGBColor))
#
#     def map_scalar(self, scalar, space=MshColor):
#         '''
#         Given a scalar value between 0 and 1, map to a color. The color is
#         returned as a sRGBColor object.
#
#         :param float scalar: The value to map to a color.
#         :param color_object space: The colormath color object to do interpolation in.
#         '''
#         if scalar < 0:
#             return convert_color(self.low_msh, sRGBColor)
#         if scalar > 1:
#             return convert_color(self.high_msh, sRGBColor)
#
#         interp = scalar
#         low_color = convert_color(self.low_msh, space)
#         high_color = convert_color(self.high_msh, space)
#         if self.midpoint_magnitude:
#             # Adjust the interpolation around the midpoint
#             if scalar < 0.5:
#                 interp = 2 * scalar
#                 high_msh = MshColor(self.midpoint_magnitude, 0, self.midpoint_low_hue,
#                                     observer=self.low_msh.observer,
#                                     illuminant=self.low_msh.illuminant)
#                 high_color = convert_color(high_msh, space)
#             else:
#                 interp = 2 * scalar - 1
#                 low_msh = MshColor(self.midpoint_magnitude, 0, self.midpoint_high_hue,
#                                    observer=self.low_msh.observer,
#                                    illuminant=self.low_msh.illuminant)
#                 low_color = convert_color(low_msh, space)
#         low_color = np.array(low_color.get_value_tuple())
#         high_color = np.array(high_color.get_value_tuple())
#
#         mid_color = interp * (high_color - low_color) + low_color
#         rgb = convert_color(space(mid_color[0], mid_color[1], mid_color[2],
#                                   observer=self.low_msh.observer,
#                                   illuminant=self.low_msh.illuminant),
#                             sRGBColor)
#
#         if ((rgb.rgb_r < -0.0019) or (rgb.rgb_r > 1.0019) or
#                 (rgb.rgb_g < -0.0019) or (rgb.rgb_g > 1.0019) or
#                 (rgb.rgb_b < -0.0019) or (rgb.rgb_b > 1.0019)):
#             print('WARNING: Value at scalar %1.4f is out of range' % scalar,
#                   rgb.get_value_tuple())
#
#         # Just in case the color leaves the color gammut, clamp to valid values.
#         return sRGBColor(rgb.clamped_rgb_r,
#                          rgb.clamped_rgb_g,
#                          rgb.clamped_rgb_b)
#
#     def map_scalar_array(self, scalar_array, space=MshColor):
#         '''
#         Given an array of scalar values between 0 and 1, map them to colors.
#         The color is returned as a sRGBColor object.
#
#         :param float scalar_array: Array of values to map to colors.
#         :param color_object space: The colormath color object to do interpolation in.
#         '''
#         f = numpy.vectorize(lambda x: self.map_scalar(x, space))
#         return f(scalar_array)