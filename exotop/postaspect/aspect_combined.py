""" Functions to plot relevant 1D thermal model results which rely on processing of ASPECT runs """
import numpy as np
import sys
# sys.path.insert(0, '/home/cmg76/Works/exo-top/')
sys.path.insert(0, '/home/claire/Works/exo-top/')
from exotop.postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    c_regimes_td, fig_fmt, \
    regime_grid_td, regime_names_td, load_grid, alpha_m, p_Earth  # noqa: E402
from exotop.postaspect.aspect_scalings import plot_save, dimensionalise_h, nondimensionalise_h, pickleio, data_path_bullard, Ra_i_eff  # noqa: E402
from exotop.model_1D.parameters import M_E  # noqa: E402
from exotop.model_1D.the_results import bulk_planets  # noqa: E402
from exotop.useful_and_bespoke import iterable_not_string, not_iterable
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from six import string_types
from collections import Iterable


def overplot_aspect_data(ax, case, x_param=None, y_param=None, pkl_suffix=None, c='k', label=None,
                         data_path=data_path_bullard, legsize=12, **kwargs):
    """ plot aspect data on existing ax for comparison, for a single case
    x and y params are either numerical value (in which case plot that) or string of df attribute stored in pkl suff"""

    # correct for differently-named variables
    if y_param == 'heuristic_h':
        y_param = 'h_components'
    elif y_param == 'dyn_top_rms':
        y_param = 'h_rms'
    Ra_str = case[2:5]
    eta_str = case[9:12]

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

    if not isinstance(x_param, Iterable):
        x_data = x_param
    elif x_param =='Ra_i_eff':
        x_data = Ra_i_eff(Ra_1=float(Ra_str), d_eta=float(eta_str), T_i=np.median(df['T_i']),
                          T_l=np.median(df['T_l']), delta_L=np.median(df['delta_L']))
    elif isinstance(x_param, string_types):
        try:
            x_data = df.x_param.mean
        except KeyError:
            raise Exception('x_data key not found in pickle')
    else:
        raise Exception('Must provide numerical x value or dataframe column')

    if not isinstance(y_param, Iterable):
        y_data = y_param
    elif isinstance(y_param, string_types):
        try:
            y_data = df.y_param.mean
        except KeyError:
            raise Exception('y_data key not found in pickle')
    else:
        raise Exception('Must provide numerical y value or dataframe column')

    print('x_data', x_data)
    print('y_data', y_data)
    ax.scatter(x_data, y_data, c=c, s=20, marker='*', label=label)
    if label is not None:
        ax.legend(frameon=False, fontsize=legsize)

    return ax


def plot_h_heuristic_variations(default='Earthbaseline', x_param='M_p', x_min=0.1*M_E, x_max=5*M_E, y_params=None,
                                xlabel=r'$M_p$ ($M_E$)', x2_param='Ra_i_eff', x2label=r'Ra$_i$', save=True, age=4.5,
                                fname='heuristic_mass_dependence',
                                ylabels=None, labelpad=10,
                                yscales=None, legend=True, xscale=M_E**-1, x2scale=1, nondimensional=True,
                                fig=None, axes=None, c='xkcd:drab green', c2='xkcd:electric purple', lw=2, labelsize=20,
                                legsize=14, ticksize=14, overplot_aspect_cases=None, logx=True, alpha=0.7,
                                xleglabel=r'$M_p$, 1D model', x2leglabel=r'Ra$_i$, 1D model', **kwargs):

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

    pl_mass = bulk_planets(n=20, name=x_param, mini=x_min, maxi=x_max, like=default, random=False,
                           nondimensional=nondimensional,
                           initial_kwargs={'tf': age}, **kwargs)

    if nondimensional and yscales is None:
        yscales = [1, 1, 1, 1]
        for yy in range(len(y_params)):
            y_params[yy] = y_params[yy] + '_prime'

    elif yscales is None:
        yscales = [1e-3, 1, 1e-3, 1e-3]

    if axes is None:
        fig, axes = plt.subplots(len(y_params), 1, figsize=(4, 3 * len(y_params)), sharex=True)

    lines = []
    for iax, ax in enumerate(axes):

        x = np.array([vars(pl)[x_param] for pl in pl_mass])*xscale
        y = np.array([vars(pl)[y_params[iax]][-1] for pl in pl_mass])*yscales[iax]
        p1, = ax.plot(x, y, c=c, lw=lw, label=xleglabel, alpha=alpha)
        ax.set_xlim(x.min(), x.max())
        if iax + 1 == len(axes):
            ax.set_xlabel(xlabel, fontsize=labelsize, labelpad=labelpad)
            lines.append(p1)
        ax.set_ylabel(ylabels[iax], fontsize=labelsize, labelpad=labelpad)

        if x2_param is not None:
            ax2 = ax.twiny()
            x2 = np.array([vars(pl)[x2_param][-1] for pl in pl_mass])*x2scale
            p2, = ax2.plot(x2, y, c=c2, lw=lw, label=x2leglabel, alpha=alpha, ls='--')
            ax2.set_xlim(x2.min(), x2.max())
            if iax == 0:
                ax2.set_xlabel(x2label, fontsize=labelsize, labelpad=labelpad)
                lines.append(p2)
                aspectlabel = '2D model'
            else:
                aspectlabel = None
            for case in overplot_aspect_cases:
                ax2 = overplot_aspect_data(ax2, case, x_param=x2_param, y_param=y_params[iax][:-5], # remove prime
                                     pkl_suffix=['_T', '_h'], c=c2, label=aspectlabel, legsize=legsize, **kwargs)

        for a in [ax, ax2]:
            if logx:
                a.set_xscale('log')
            a.tick_params(axis='both', labelsize=ticksize)
            # tkw = dict(size=4, width=1.5)
            # a.tick_params(axis='y', colors=p1.get_color(), **tkw)
    if legend:
        axes[0].legend(lines, [l.get_label() for l in lines], fontsize=legsize, frameon=False)
    if save:
        plot_save(fig, fname, **kwargs)

