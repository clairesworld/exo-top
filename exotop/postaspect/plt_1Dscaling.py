""" Functions to plot relevant 1D thermal model results which rely on processing of ASPECT runs """
import numpy as np
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import NullFormatter  # noqa: E402
import matplotlib.lines as mlines  # noqa: E402
from six import string_types  # noqa: E402
from collections import Iterable  # noqa: E402
import matplotlib.animation as animation
import matplotlib.ticker as ticker
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
import sys

# sys.path.insert(0, '/home/cmg76/Works/exo-top/')
# sys.path.insert(0, '/home/claire/Works/exo-top/')
from postaspect.plt_aspect import plot_save, plot_getx  # noqa: E402
from postaspect.aspect_post import dimensionalise_h, nondimensionalise_h, pickleio, Ra_i_eff  # noqa: E402
from postaspect.setup_postprocessing import data_path_bullard, fig_path_bullard  # noqa: E402
from model_1D.parameters import M_E, sec2Gyr  # noqa: E402
from model_1D.the_results import bulk_planets, plot_output, build_planet_from_id  # noqa: E402
from useful_and_bespoke import find_nearest_idx, iterable_not_string, not_iterable, colorize, dark_background


def overplot_aspect_data(ax, case, x_param=None, y_param=None, pkl_suffix=['_T', '_h'], c='k', markersize=20,
                         marker='*', zorder=200,
                         data_path=data_path_bullard, return_data=False, visible=True, **kwargs):
    """ plot aspect data on existing ax for comparison, for a single case
    x and y params are either numerical value (in which case plot that) or string of df attribute stored in pkl suff"""
    # TODO avg

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

    # x_data = plot_getx(Ra_str, eta_str, case=case, which_x=x_param, averagescheme='timefirst', data_path=data_path,
    #                    load=True, postprocess_kwargs=None, return_all=False, t1=0, **kwargs)
    # y_data = plot_getx(Ra_str, eta_str, case=case, which_x=y_param, averagescheme='timefirst', data_path=data_path,
    #                    load=True, postprocess_kwargs=None, return_all=False, t1=0, **kwargs)

    ax.scatter(x_data, y_data, c=c, s=markersize, marker=marker, visible=visible, zorder=zorder)
    if return_data:
        return ax, x_data, y_data
    else:
        return ax


def plot_h_heuristic_variations(default='Earthbaseline', x_param='M_p', x_min=0.1 * M_E, x_max=10 * M_E, y_params=None,
                                xlabel=r'$M_p$ ($M_E$)', x2_param='Ra_i_eff', x2label=r'Ra$_i$', save=True, age=4.5,
                                fname='heuristic_mass_dependence', overplot_aspect_cases_x_param='Ra_i_eff',
                                ylabels=None, labelpad=10, nplanets=20, markersize=20, marker='*', ini_dict=None,
                                yscales=None, legend=True, xscale=M_E ** -1, x2scale=1, nondimensional=True,
                                fig=None, axes=None, c='xkcd:drab green', c2='xkcd:electric purple', lw=2, labelsize=20,
                                legsize=14, ticksize=14, overplot_aspect_cases=None, logx=True, alpha=0.7,
                                xleglabel=r'$M_p$, 1D model', x2leglabel=r'Ra$_i$, 1D model', return_artists=False,
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


# # TODO
# def animate_Ra(default='Earthbaseline', fig_path='plat/', figsize=(9, 9), labelsize=16, ylabelpad=10, xlabelpad=10,
#                     ticksize=12, fname='ani_1D',
#                     xticks=[0.01, 0.1, 0.3, 1, 2, 3, 4, 5, 6], yticks=None, aspect_cases=None, markersize=20, marker='o', c_scat='g',
#                     data_path='', x_min=0.03 * M_E, x_max=6 * M_E, y_param=None, x_param='M_p', x2_param='Ra_i_eff',
#                     xscale=M_E ** -1, x2scale=1, yscale=1, fps=15, c='b',
#                     ani_param='Ea', ani_tpref=r'$E_a$ = ', ani_tsuff=r' kJ mol$^{-1}$', ani_scale=1e-3, ani_min=200, ani_max=375, ani_res=30,
#                     tf=9.9, ti=2, ylabel='', x_min_planets=0.1 * M_E, ini_dict=None, logx=True, text=True, dark=True,
#                     lw=3, x_res=8, alpha_lines=1, xlabel=r'$M_p$ ($M_E$)', x2label=r'Ra$_{i, eff}$', **kwargs):
#     # params refer to variable names in 1D model
#     if aspect_cases is None:
#         aspect_cases = ['Ra3e8-eta1e8-wide-ascii', 'Ra3e8-eta1e7-wide-ascii', 'Ra3e8-eta1e6-wide',
#                         'Ra1e8-eta1e8-wide-ascii',
#                         'Ra1e8-eta1e7-wide', 'Ra1e8-eta1e6-wide']
#     ## too much effort to do animated subplots (and massive files)
#
#     # setup plot
#     fig, ax = plt.subplots(figsize=figsize)
#
#
#     # only need to load aspect data once and then re-dimensionalise it
#     cases_x, cases_y = [], []
#     for case in aspect_cases:
#         ax, x_aspect_data, y_aspect_data = overplot_aspect_data(ax, case, x_param=x2_param,
#                                                                 y_param=y_param, data_path=data_path,
#                                                                 pkl_suffix=['_T', '_h'], markersize=markersize,
#                                                                 marker=marker,
#                                                                 return_data=True, visible=False, **kwargs)
#         cases_x.append(x_aspect_data)
#         cases_y.append(y_aspect_data)
#     scat_data = np.vstack((cases_x, cases_y)).T  # shape N, 2
#
#     # setup evolutions
#     ani_vec = np.linspace(ani_min*ani_scale, ani_max*ani_scale, ani_res)
#     if x_min_planets is None:
#         x_min_planets = x_min
#     initial_kwargs = {'tf': tf}
#     if ini_dict is not None:
#         initial_kwargs.update(ini_dict)
#     pl_min = build_planet_from_id(ident=default, nondimensional=True, initial_kwargs=initial_kwargs, verbose=True,
#                                   update_kwargs={x_param: x_min_planets}, **kwargs)
#     i_start_time = np.argmax(pl_min.t >= ti / sec2Gyr)  # don't animate any cases before here
#
#     # run evolutions
#     pl_mass = bulk_planets(n=x_res, name=x_param, mini=x_min, maxi=x_max, like=default, random=False,
#                            nondimensional=True, logscale=logx, initial_kwargs=initial_kwargs, t_eval=pl_min.t,
#                            **kwargs)
#
#     # only plot planets with x value greater than x_min_planets (this is to extend axis limits for aspect compare)
#     x_orig = np.array([vars(pl)[x_param] for pl in pl_mass]) * xscale
#     ax.set_xlim(x_orig.min(),
#                 x_orig.max())  # need to set limits for twin axis before removing small points so u can see Ra
#     i_plot = np.argmax(x_orig >= x_min_planets * xscale)
#
#     pl_pl = pl_mass[i_plot:]
#     x = x_orig[i_plot:]
#     y = np.array([vars(pl)[y_param + '_prime'][i_start_time] for pl in pl_pl]) * yscale
#     p1, = ax.plot(x, y, c=c, lw=lw, alpha=alpha_lines)
#
#     # setup secondary axis
#     ax2 = ax.twiny()
#     x2_orig = np.array([vars(pl)[x2_param][i_start_time] for pl in
#                         pl_mass]) * x2scale  # limits corresponding to tiny planets for extra space
#     ax2.set_xlim(x2_orig.min(), x2_orig.max())
#     x2 = x2_orig[i_plot:]
#     p2, = ax2.plot(x2, y, c=c, lw=lw, alpha=alpha_lines, ls='--')
#
#     # add scatter points for 2D simulations
#     scat = ax2.scatter(cases_x, cases_y, c=c_scat, s=markersize, marker=marker, visible=True)
#
#     # set plot params
#     if xticks is not None:
#         ax.set_xticks(xticks)
#     if yticks is not None:
#         ax.set_yticks(yticks)
#     if logx:
#         ax2.set_xscale('log')
#         ax2.set_yscale('log')
#         ax.set_xscale('log')
#         ax.set_yscale('log')
#     ax2.set_xlabel(x2label, fontsize=labelsize, labelpad=xlabelpad)
#     ax2.tick_params(axis='x', which='major', labelsize=ticksize)
#     ax.set_xlabel(xlabel, fontsize=labelsize, labelpad=xlabelpad)
#     ax.set_ylabel(ylabel, fontsize=labelsize, labelpad=ylabelpad)
#     ax.tick_params(axis='both', which='major', labelsize=ticksize)
#     # ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
#     # ax.xaxis.set_minor_formatter(ticker.ScalarFormatter())
#     ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
#
#     # # set consistent tick labels - do this by tracking a few planets of known mass - will work when lines match
#     # x2_desired_ticks = np.array([1e6, 1e7, 1e8, 1e9, 1e10])
#     # x_at_ticks = np.zeros_like(x2_desired_ticks)
#     # for ii, ttt in enumerate(x2_desired_ticks):
#     #     idx = find_nearest_idx(x2_orig, ttt)
#     #     x_at_ticks[ii] = x_orig[idx]
#     # print('idx', idx)
#     #
#     # print('mass at desired ticks', x_at_ticks)
#     # print('Ra at desired ticks', x2_orig[idx.astype(int)])
#
#     # colour and adjust limits
#     if dark:
#         fig, ax, ax2 = dark_background(fig, [ax, ax2])
#         c_text = 'xkcd:off white'
#     else:
#         c_text = 'k'
#     if text:
#         ani_text = ax.text(0.95, 0.95,
#                            ani_tpref + '{:4.1f}'.format(anim_vec[i_start_time] * ani_scale) + ani_tsuff,
#                            fontsize=ticksize, c=c_text,
#                            horizontalalignment='right',
#                            verticalalignment='top',
#                            transform=ax.transAxes)
#     else:
#         ani_text = ax.text([], [], '')
#     plt.subplots_adjust(bottom=0.2, top=0.8, left=0.15)
#     plt.tight_layout()
#
#     # animation stuff
#     def init():
#         p1.set_ydata(([np.nan] * len(x)))
#         p2.set_ydata(([np.nan] * len(x)))
#         # scat.set_offsets([[], []])
#         ani_text.set_text('')
#         return p1, p2, ani_text,
#
#     def animate(i, cases_x, cases_y, pl_pl, y_param, x2_param, yscale, x2scale, pl_mass, anim_vec, anim_tpref,
#                 anim_scale, anim_tsuff):
#         new_cases_x = cases_x  # Ra_i_eff axis - doesn't move if axis isn't moving
#         new_cases_y = cases_y
#         # new_scat_data = np.vstack((new_cases_x, new_cases_y)).T  # shape N, 2
#         p1y = np.array([vars(pl)[y_param + '_prime'][i] for pl in pl_pl]) * yscale
#         p2x = np.array([vars(pl)[x2_param][i] for pl in pl_pl]) * x2scale
#         p2y = np.array([vars(pl)[y_param + '_prime'][i] for pl in pl_pl]) * yscale
#
#         p1.set_ydata(p1y)  # update the data.
#         p2.set_xdata(p2x)  # note this is attached to Ra axis
#         p2.set_ydata(p2y)
#         # scat.set_offsets(new_scat_data)
#
#         p2xx = np.array([vars(pl)[x2_param][i] for pl in
#                          pl_mass]) * x2scale  # new x limits are masses corresponding to tiny extra planets
#         ax2.set_xlim(min(p2xx), max(p2xx))  # added ax attribute here
#
#         ani_text.set_text(anim_tpref + '{:4.1f}'.format(anim_vec[i] * anim_scale) + anim_tsuff)
#         return p1, p2, ani_text,
#
#     ani = animation.FuncAnimation(fig, animate, frames=range(i_start_time, len(anim_vec)), init_func=init,
#                                   fargs=(cases_x, cases_y, pl_pl, y_param, x2_param, yscale, x2scale, pl_mass, anim_vec,
#                                          ani_tpref, ani_scale, ani_tsuff),
#                                   blit=False, repeat=True,
#                                   interval=200, )  # interval: delay between frames in ms
#     print('finished!')
#     ani.save(fig_path + fname + '.gif', writer='imagemagick', fps=fps,
#              savefig_kwargs={'facecolor': fig.get_facecolor()})
#     print('saved to', fig_path + fname)


def animate_Ra_time(default='Earthbaseline', fig_path='plat/', figsize=(9, 9), labelsize=16, ylabelpad=10, xlabelpad=10,
                    ticksize=12, fname='ani_1D', scalar=False, i_static=None, x_test2=None,
                    xticks=[0.1, 0.3, 1, 2, 3, 4, 5, 6], yticks=None, aspect_cases=None, ms_scat=20, marker_scat='o',
                    c_scat='g',
                    data_path='', x_test=0.815 * M_E, x_min=0.03 * M_E, x_max=6 * M_E, y_param=None, x_param='M_p',
                    x2_param='Ra_i_eff',
                    xscale=M_E ** -1, x2scale=1, yscale=1, fps=15, c='b', y2label=r'$\Delta h_{rms}$ (m)', y2scale=1,
                    anim_param='t', anim_tpref='t = ', anim_tsuff=' Gyr', anim_scale=sec2Gyr, secy=False,
                    tf=10, ti=2, ylabel='', x_min_planets=0.1 * M_E, ini_dict=None, logx=True, text=True, dark=True,
                    lw=3, nplanets=20, alpha_lines=1, xlabel=r'$M_p$ ($M_E$)', x2label=r'Ra$_{i, eff}$', **kwargs):
    # params refer to variable names in 1D model
    if aspect_cases is None:
        aspect_cases = ['Ra3e8-eta1e8-wide-ascii', 'Ra3e8-eta1e7-wide-ascii', 'Ra3e8-eta1e6-wide',
                        'Ra1e8-eta1e8-wide-ascii',
                        'Ra1e8-eta1e7-wide', 'Ra1e8-eta1e6-wide']
    ## too much effort to do animated subplots (and massive files)

    fig, ax = plt.subplots(figsize=figsize)

    # only need to load aspect data once and then re-dimensionalise it
    cases_x, cases_y = [], []
    for case in aspect_cases:
        ax, x_aspect_data, y_aspect_data = overplot_aspect_data(ax, case, x_param=x2_param,
                                                                y_param=y_param, data_path=data_path,
                                                                pkl_suffix=['_T', '_h'], markersize=ms_scat,
                                                                marker=marker_scat,
                                                                return_data=True, visible=False, **kwargs)
        cases_x.append(x_aspect_data)
        cases_y.append(y_aspect_data)
    scat_data = np.vstack((cases_x, cases_y)).T  # shape N, 2

    # run evolutions
    if x_min_planets is None:
        x_min_planets = x_min
    initial_kwargs = {'tf': tf}
    if ini_dict is not None:
        initial_kwargs.update(ini_dict)
    pl_test = build_planet_from_id(ident=default, nondimensional=True, initial_kwargs=initial_kwargs, verbose=True,
                                   update_kwargs={x_param: x_test}, **kwargs)
    if x_test2 is not None:
        pl_test2 = build_planet_from_id(ident=default, nondimensional=True, initial_kwargs=initial_kwargs, verbose=True,
                                        update_kwargs={x_param: x_test2}, **kwargs)
    pl_mass = bulk_planets(n=nplanets, name=x_param, mini=x_min, maxi=x_max, like=default, random=False,
                           nondimensional=True, logscale=logx, initial_kwargs=initial_kwargs, t_eval=pl_test.t,
                           **kwargs)
    i_start_time = np.argmax(pl_test.t >= ti / sec2Gyr)  # don't animate any cases before here
    if i_static is None:
        # index to make static figure
        i_static = i_start_time

    if anim_param == 't':
        anim_vec = np.array(pl_mass[-1].t)

    # only plot planets with x value greater than x_min_planets (this is to extend axis limits for aspect compare)
    x_orig = np.array([vars(pl)[x_param] for pl in pl_mass]) * xscale
    ax.set_xlim(x_orig.min(),
                x_orig.max())  # need to set limits for twin axis before removing small points so u can see Ra
    i_plot = np.argmax(x_orig >= x_min_planets * xscale)
    pl_pl = pl_mass[i_plot:]
    x = x_orig[i_plot:]
    y = np.array([vars(pl)[y_param + '_prime'][i_static] for pl in pl_pl]) * yscale
    p1, = ax.plot(x, y, c=c, lw=lw, alpha=alpha_lines)
    pl_scat = ax.scatter(vars(pl_test)[x_param] * xscale, vars(pl_test)[y_param + '_prime'][i_static] * yscale,
                         marker='$V$',
                         s=ms_scat * 2,
                         c='xkcd:goldenrod', zorder=100)
    if x_test2 is not None:
        pl_scat2 = ax.scatter(vars(pl_test2)[x_param] * xscale, vars(pl_test2)[y_param + '_prime'][i_static] * yscale,
                              marker='$M$',
                              s=ms_scat * 2,
                              c='xkcd:goldenrod', zorder=100)
    else:
        pl_scat2 = ax.scatter([], [], visible=False)

    # secondary x axis

    # attempt to scale x2 axis
    x2_orig = np.array([vars(pl)[x2_param][i_static] for pl in
                        pl_mass]) * x2scale  # limits corresponding to tiny planets for extra space
    x2 = x2_orig[i_plot:]
    global x2_orig_i
    x2_orig_i = x2_orig

    def mass2Ra(u):
        f = InterpolatedUnivariateSpline(x_orig, x2_orig_i, k=3)
        return f(u)

    def Ra2mass(u):
        f = InterpolatedUnivariateSpline(x2_orig_i, x_orig, k=3)
        return f(u)

    ax2 = ax.secondary_xaxis('top', functions=(mass2Ra, Ra2mass))
    p2, = ax.plot(Ra2mass(x2), y, c=c, lw=lw, alpha=alpha_lines, ls='--')
    aspect_scat = ax.scatter(Ra2mass(cases_x), cases_y, c=c_scat, s=ms_scat, marker=marker_scat, zorder=100,
                             visible=True)

    # ax2 = ax.twiny()
    # ax2.set_xlim(x2_orig.min(), x2_orig.max())
    # p2, = ax2.plot(x2, y, c=c, lw=lw, alpha=alpha_lines, ls='--')
    # aspect_scat = ax2.scatter(cases_x, cases_y, c=c_scat, s=markersize, marker=marker, visible=True)

    # secondary y axis
    global alpha_m
    global dT_i
    global dm
    alpha_m = pl_test.alpha_m
    dT_i = pl_test.T_scale[i_static]
    dm = pl_test.L_scale

    def h_todim(u):
        return u * alpha_m * dT_i * dm * y2scale

    def h_tonondim(u):
        return u / (alpha_m * dT_i * dm * y2scale)

    if secy:
        ax3 = ax.secondary_yaxis('right', functions=(h_todim, h_tonondim))
        ax3.set_ylabel(y2label, fontsize=labelsize, labelpad=ylabelpad)
        ax3.tick_params(axis='y', which='major', labelsize=ticksize)

    # setup plot params
    ax2.set_xlabel(x2label, fontsize=labelsize, labelpad=xlabelpad)
    ax2.tick_params(axis='x', which='major', labelsize=ticksize)
    ax.set_xlabel(xlabel, fontsize=labelsize, labelpad=xlabelpad)
    ax.set_ylabel(ylabel, fontsize=labelsize, labelpad=ylabelpad)
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    # ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    # ax.xaxis.set_minor_formatter(ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
    ax2.yaxis.set_minor_formatter(ticker.NullFormatter())
    if logx:
        ax.set_xscale('log')
        ax.set_yscale('log')
        # ax2.set_xscale('log')
        # ax2.set_yscale('log')
    if xticks is not None:
        ax.set_xticks(xticks)
    if yticks is not None:
        ax.set_yticks(yticks)

    # colour and adjust
    if dark:
        if secy:
            fig, ax, ax2, ax3 = dark_background(fig, [ax, ax2, ax3])
        else:
            fig, ax, ax2 = dark_background(fig, [ax, ax2])
        c_text = 'xkcd:off white'
    else:
        c_text = 'k'
    if text:
        ani_text = ax.text(0.95, 0.95,
                           anim_tpref + '{:4.1f}'.format(anim_vec[i_static] * anim_scale) + anim_tsuff,
                           fontsize=ticksize, c=c_text,
                           horizontalalignment='right',
                           verticalalignment='top',
                           transform=ax.transAxes)
    else:
        ani_text = ax.text(0.5, 0.5, '', alpha=0, transform=ax.transAxes)

    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())

    plt.subplots_adjust(bottom=0.2, top=0.8, left=0.15)
    plt.tight_layout()

    # save non-animated snapshot
    plot_save(fig, fname=fname + '_static', fig_path=fig_path, bbox_inches=None)

    # animation stuff
    def init():
        p1.set_ydata(([np.nan] * len(x)))
        p2.set_ydata(([np.nan] * len(x)))
        aspect_scat.set_offsets([[], []])
        pl_scat.set_offsets([[], []])
        if x_test2 is not None:
            pl_scat2.set_offsets([[], []])
        ani_text.set_text('')
        return p1, p2, aspect_scat, pl_scat, ani_text, pl_scat2

    def animate(i, cases_x, cases_y, pl_pl, y_param, x2_param, yscale, x2scale, pl_mass, anim_vec, anim_tpref,
                anim_scale, anim_tsuff, pl_scat2):

        # update secondary x
        global x2_orig_i
        x2_orig_i = np.array([vars(pl)[x2_param][i] for pl in pl_mass]) * x2scale  # new x2 based on original size
        new_cases_x = Ra2mass(cases_x)  # Ra_i_eff axis - doesn't move if axis isn't moving
        new_cases_y = cases_y
        new_scat_data = np.vstack((new_cases_x, new_cases_y)).T  # shape N, 2

        # update secondary y
        global dT_i
        dT_i = pl_test.T_scale[i]

        p1y = np.array([vars(pl)[y_param + '_prime'][i] for pl in pl_pl]) * yscale
        # p2x = np.array([vars(pl)[x2_param][i] for pl in pl_pl]) * x2scale
        # p2y = np.array([vars(pl)[y_param + '_prime'][i] for pl in pl_pl]) * yscale

        p1.set_ydata(p1y)  # update the data.
        # p2.set_xdata(Ra2mass(p2x))  # note this is attached to Ra axis
        # p2.set_ydata(p2y)
        aspect_scat.set_offsets(new_scat_data)

        pl_scat_x = vars(pl_test)[x_param] * xscale
        pl_scat_y = vars(pl_test)[y_param + '_prime'][i] * yscale
        pl_scat_data = np.vstack((pl_scat_x, pl_scat_y)).T  # shape N, 2
        pl_scat.set_offsets(pl_scat_data)

        if pl_scat2 is not None:
            pl_scat2_x = vars(pl_test2)[x_param] * xscale
            pl_scat2_y = vars(pl_test2)[y_param + '_prime'][i] * yscale
            pl_scat2_data = np.vstack((pl_scat2_x, pl_scat2_y)).T  # shape N, 2
            pl_scat2.set_offsets(pl_scat2_data)
        else:
            pl_scat2 = None

        # ax2.set_xlim(min(x2_orig_i), max(x2_orig_i))  # added ax attribute here
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
        ax.xaxis.set_minor_formatter(ticker.NullFormatter())
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
        ax.yaxis.set_minor_formatter(ticker.NullFormatter())

        ani_text.set_text(anim_tpref + '{:4.1f}'.format(anim_vec[i] * anim_scale) + anim_tsuff)
        return p1, p2, aspect_scat, pl_scat, ani_text, pl_scat2

    # print('i frames', i_start_time, len(anim_vec), 'len', len(range(i_start_time, len(anim_vec))))
    # print('len x2', len(vars(pl_mass[-1])[x2_param]))
    # print('len y', len(vars(pl_mass[-1])[y_param + '_prime']))
    # for ii in range(len(pl_mass)):
    #     print('   ii, x2:', len(vars(pl_mass[ii])[x2_param]), pl_mass[ii].M_p/M_E)

    ani = animation.FuncAnimation(fig, animate, frames=range(i_start_time, len(anim_vec)), init_func=init,
                                  fargs=(cases_x, cases_y, pl_pl, y_param, x2_param, yscale, x2scale, pl_mass, anim_vec,
                                         anim_tpref, anim_scale, anim_tsuff, pl_scat2),
                                  blit=False, repeat=True,
                                  interval=200, )  # interval: delay between frames in ms

    print('finished!')
    ani.save(fig_path + fname + '.gif', writer='imagemagick', fps=fps,
             savefig_kwargs={'facecolor': fig.get_facecolor()})
    print('saved to', fig_path + fname)
