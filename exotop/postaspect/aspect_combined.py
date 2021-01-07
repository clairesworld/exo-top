""" Functions to plot relevant 1D thermal model results which rely on processing of ASPECT runs """
import numpy as np
import sys

# sys.path.insert(0, '/home/cmg76/Works/exo-top/')
sys.path.insert(0, '/home/claire/Works/exo-top/')
from exotop.postaspect.aspect_scalings import plot_save, dimensionalise_h, nondimensionalise_h, pickleio, \
    data_path_bullard, Ra_i_eff  # noqa: E402
from exotop.model_1D.parameters import M_E  # noqa: E402
from exotop.model_1D.the_results import bulk_planets, plot_output  # noqa: E402
from exotop.useful_and_bespoke import iterable_not_string, not_iterable, colorize
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import NullFormatter  # noqa: E402
import matplotlib.lines as mlines  # noqa: E402
from six import string_types  # noqa: E402
from collections import Iterable  # noqa: E402


def overplot_aspect_data(ax, case, x_param=None, y_param=None, pkl_suffix=None, c='k', markersize=20, marker='*',
                         data_path=data_path_bullard, returndf=False, **kwargs):
    """ plot aspect data on existing ax for comparison, for a single case
    x and y params are either numerical value (in which case plot that) or string of df attribute stored in pkl suff"""

    # correct for differently-named variables
    if y_param == 'heuristic_h':
        y_param = 'h_components'
    elif y_param == 'dyn_top_rms':
        y_param = 'h_rms'
    Ra_str = case[2:5]
    eta_str = case[9:12]

    # first load dataframe if given
    if isinstance(pkl_suffix, Iterable):
        if isinstance(pkl_suffix, string_types):  # single pickle to load
            df = pickleio(case, suffix=pkl_suffix, postprocess_functions=None, load=True, data_path=data_path)
        else:  # multiple pickles to load
            dfs = []
            for ip, ps in enumerate(pkl_suffix):
                df1 = pickleio(case, suffix=ps, postprocess_functions=None, load=True, data_path=data_path)
                dfs.append(df1)
            df = pd.concat(dfs, axis=1)
            df = df.loc[:, ~df.columns.duplicated()]

    # get x value
    if not isinstance(x_param, Iterable):  # i.e. not a string key
        x_data = x_param
    elif isinstance(pkl_suffix, Iterable) and isinstance(x_param, string_types):  # provided a key
        if x_param == 'Ra_i_eff':
            x_data = Ra_i_eff(Ra_1=float(Ra_str), d_eta=float(eta_str), T_i=np.median(df['T_i']),
                              T_l=np.median(df['T_l']), delta_L=np.median(df['delta_L']))
        else:
            try:
                x_data = df[x_param].mean(skipna=True)
            except KeyError:
                raise Exception('x_data key ' + x_param + ' not found in pickle')
    else:
        raise Exception('Must provide numerical x value or dataframe column + pickle file id')

    # get y value
    if not isinstance(y_param, Iterable):
        y_data = y_param
    elif isinstance(pkl_suffix, Iterable) and isinstance(y_param, string_types):
        try:
            y_data = df[y_param].mean(skipna=True)
        except KeyError:
            raise Exception('y_data key ' + y_param + ' not found in pickle')
    else:
        raise Exception('Must provide numerical y value or dataframe column + pickle file id')

    print('(x,y)', x_data, y_data)
    ax.scatter(x_data, y_data, c=c, s=markersize, marker=marker)
    if returndf:
        return ax, df
    else:
        return ax


def plot_h_heuristic_variations(default='Earthbaseline', x_param='M_p', x_min=0.1 * M_E, x_max=10 * M_E, y_params=None,
                                xlabel=r'$M_p$ ($M_E$)', x2_param='Ra_i_eff', x2label=r'Ra$_i$', save=True, age=4.5,
                                fname='heuristic_mass_dependence', overplot_aspect_cases_x_param='Ra_i_eff',
                                ylabels=None, labelpad=10, nplanets=20, markersize=20, marker='*', ini_dict=None,
                                yscales=None, legend=True, xscale=M_E ** -1, x2scale=1, nondimensional=True,
                                fig=None, axes=None, c='xkcd:drab green', c2='xkcd:electric purple', lw=2, labelsize=20,
                                legsize=14, ticksize=14, overplot_aspect_cases=None, logx=True, alpha=0.7,
                                xleglabel=r'$M_p$, 1D model', x2leglabel=r'Ra$_i$, 1D model',
                                aspectleglabel=r'Ra$_i$, 2D model', x_min_planets=None, **kwargs):
    if nondimensional and ylabels is None:
        ylabels = [r'$\Delta h^\prime_{rms}$', r'$\Delta T^\prime_{rh}$', r'$\delta^\prime_u$',
                   r'$\alpha^\prime \Delta T^\prime_{rh} \delta^\prime_u$']
    elif ylabels is None:
        ylabels = [r'$\Delta h_{rms}$ (km)', r'$\Delta T_{rh}$ (K)', r'$\delta_u$ (km)',
                   r'$\alpha \Delta T_{rh} \delta_u$ (km)']
    if y_params is None:
        y_params = ['dyn_top_rms', 'dT_rh', 'delta_rh', 'heuristic_h']
    if overplot_aspect_cases is None:
        overplot_aspect_cases = []
    if x_min_planets is None:
        x_min_planets = x_min

    initial_kwargs = {'tf': age}
    if ini_dict is not None:
        initial_kwargs.update(ini_dict)

    pl_mass = bulk_planets(n=nplanets, name=x_param, mini=x_min, maxi=x_max, like=default, random=False,
                           nondimensional=nondimensional, logscale=logx, initial_kwargs=initial_kwargs,
                           **kwargs)

    if nondimensional and yscales is None:
        yscales = [1, 1, 1, 1]
        for yy in range(len(y_params)):
            y_params[yy] = y_params[yy] + '_prime'

    elif yscales is None:
        yscales = [1e-3, 1, 1e-3, 1e-3]

    if axes is None:
        fig, axes = plt.subplots(len(y_params), 1, figsize=(4, 3 * len(y_params)))

    lines = []
    for row, ax in enumerate(axes):
        ax.set_ylabel(ylabels[row], fontsize=labelsize, labelpad=labelpad)
        x = np.array([vars(pl)[x_param] for pl in pl_mass]) * xscale
        y = np.array([vars(pl)[y_params[row]][-1] for pl in pl_mass]) * yscales[row]
        # only plot planets with x value greater than x_min_planets (this is to extend axis limits for aspect compare)
        i_plot = np.argmax(x >= x_min_planets * xscale)
        p1, = ax.plot(x[i_plot:], y[i_plot:], c=c, lw=lw, label=xleglabel, alpha=alpha)
        # ax.scatter(x[i_plot:], y[i_plot:], c=c, lw=0, alpha=alpha, marker='v')
        ax.set_xlim(x.min(), x.max())
        if logx:
            ax.set_xscale('log')
        if row + 1 == len(axes):  # only label bottom row
            ax.tick_params(axis='both', which='both', labelsize=ticksize)
            ax.set_xlabel(xlabel, fontsize=labelsize, labelpad=labelpad)
            lines.append(p1)
        else:
            ax.tick_params(axis='x', which='both', labelbottom='off')
            ax.tick_params(axis='y', which='both', labelsize=ticksize)
            ax.xaxis.set_major_formatter(NullFormatter())

        if x2_param is not None:
            ax2 = ax.twiny()
            x2 = np.array([vars(pl)[x2_param][-1] for pl in pl_mass]) * x2scale
            p2, = ax2.plot(x2[i_plot:], y[i_plot:], c=c2, lw=lw, label=x2leglabel, alpha=alpha, ls='--')
            # ax2.scatter(x2[i_plot:], y[i_plot:], c=c2, lw=0, alpha=alpha, marker='^')
            ax2.set_xlim(x2.min(), x2.max())
            if logx:
                ax2.set_xscale('log')
            if row == 0:  # only label top row
                ax2.tick_params(axis='both', which='both', labelsize=ticksize)
                ax2.set_xlabel(x2label, fontsize=labelsize, labelpad=labelpad)
                lines.append(p2)
            else:
                ax2.tick_params(axis='x', which='both', labeltop='off')
                ax2.xaxis.set_major_formatter(NullFormatter())

        if not (not overplot_aspect_cases):
            if overplot_aspect_cases_x_param == x_param:
                ax_op = ax
                c_op = c
            elif overplot_aspect_cases_x_param == x2_param:
                ax_op = ax2
                c_op = c2
            else:
                raise Exception('x param for overplotting aspect data must match ax or ax2')
            for case in overplot_aspect_cases:
                ax_op = overplot_aspect_data(ax_op, case, x_param=overplot_aspect_cases_x_param,
                                             y_param=y_params[row][:-6],  # remove prime, assume nondimensional
                                             pkl_suffix=['_T', '_h'], c=c_op, markersize=markersize, marker=marker,
                                             **kwargs)

        # for a in axs:
        # tkw = dict(size=4, width=1.5)
        # a.tick_params(axis='y', colors=p1.get_color(), **tkw)
    if legend:
        if not (not overplot_aspect_cases):
            p3 = mlines.Line2D([], [], lw=0, color=c_op, marker=marker, markersize=markersize / 4, label=aspectleglabel)
            lines.append(p3)
        axes[0].legend(lines, [l.get_label() for l in lines], fontsize=legsize, frameon=True)
    fig.subplots_adjust(hspace=0.05)
    if save:
        plot_save(fig, fname, tight_layout=False, **kwargs)


def plot_1D_evolutions(default, nplanets=45, labelsize=23, ticksize=16, clabelpad=4, save=True,
                       fname='evol_summary', zmin=None, zmax=None, zname=None, zlabel=None, zscale=1, ylabels_dict=None,
                       age=4.5, ini_dict=None, cmap='rainbow',
                       **kwargs):
    if ylabels_dict is None:
        ylabels_dict = {'T_m': ('$T_i$ (K)', 1),
                        'T_scale': ('$\Delta T$ (K)', 1),
                        'eta_m': ('$\eta (T_m)$', 1),
                        'delta_eta': ('$\eta_0/\eta_{\Delta T}$', 1),
                        'delta_rh': ('$\delta_{u}$ (km)', 1e-3),
                        'd_m': ('$d_m$ (km)', 1e-3),
                        'D_l': ('$D_l$ (km)', 1e-3),
                        'Ra_i': ('Ra$_i$', 1),
                        # 'H_rad_m': ('$H_{rad}$ (TW)', 1e-12)
                        }
    initial_kwargs = {'tf': age}
    if ini_dict is not None:
        initial_kwargs.update(ini_dict)
    planets = bulk_planets(nplanets, name=zname, mini=zmin, maxi=zmax, like=default, initial_kwargs=initial_kwargs,
                           nondimensional=True, **kwargs)
    z_vec = np.array([vars(pl)[zname] for pl in planets]) * zscale
    cols = colorize(z_vec, cmap=cmap)[0]

    for ii, pl in enumerate(planets):
        if ii == 0:
            fig = None
            axes = None
        fig, axes = plot_output(pl, ylabels_dict, ncols=2, plots_save=False, fig=fig, axes=axes,
                                legend=False,
                                title='', line_args={'c': cols[ii]}, **kwargs)

    sc = fig.gca().scatter(None, None, vmin=zmin, vmax=zmax, cmap='rainbow')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=zmin * zscale, vmax=zmax * zscale))
    sm._A = []
    fig.subplots_adjust(right=0.85, wspace=0.5)
    cbar_ax = fig.add_axes([0.9, 0.15, 0.04, 0.7])
    fig.colorbar(sm, cax=cbar_ax)
    cbar_ax.set_ylabel(zlabel, fontsize=labelsize, labelpad=clabelpad)
    cbar_ax.tick_params(axis='y', labelsize=ticksize)
    if save:
        plot_save(fig, fname, tight_layout=False, **kwargs)
    return fig, axes
