""" Run 1D model and plot results """

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from . import inputs
from . import terrestrialplanet as tp
from . import parameters
from . import thermal
from . import topography
from . import oceans
from collections.abc import Iterable
import collections
import six
import pandas as pd
import random as rand
from scipy import interpolate

sys.path.append("..")
import exotop.asharms as harm  # noqa: E402
import matplotlib.animation as animation
from exotop.useful_and_bespoke import age_index, dark_background, not_iterable, colorize, colourbar  # noqa: E402
from exotop.model_1D.parameters import M_E  # noqa: E402
# np.seterr('raise')
import matplotlib.ticker as ticker


def bulk_planets(n=1, name=None, mini=None, maxi=None, like=None, t_eval=None, random=False,
                 initial_kwargs=None, postprocessors=None, update_kwargs=None, logscale=False, **kwargs):
    """varying single parameter 'name' between mini and maxi, use default values otherwise.
    update_kwargs can include any TerrestrialPlanet attribute
    initial_kwargs can include T_m0, T_c0, D_l0, t0, tf"""

    if like is not None:
        pl_kwargs = eval('inputs.' + like + '_in')
        model_kwargs = eval('inputs.' + like + '_run')
    else:
        pl_kwargs = {}  # use defaults given in terrestrialplanet.py
        model_kwargs = {}  # initial conditions defaults given in thermal.py
    if update_kwargs is not None:
        pl_kwargs.update(update_kwargs)
    if initial_kwargs is not None:
        model_kwargs.update(initial_kwargs)
    planets = []
    ii = 0
    if logscale:
        arr = np.logspace(np.log10(mini), np.log10(maxi), num=n)
    else:
        arr = np.linspace(mini, maxi, num=n)
    while ii < n:
        if random:
            val = rand.uniform(mini, maxi)
        else:
            val = arr[ii]
        new_kwargs = pl_kwargs.copy()
        new_kwargs.update({name: val})

        pl = build_planet(new_kwargs, model_kwargs, postprocessors=postprocessors, t_eval=t_eval, **kwargs)
        planets.append(pl)
        ii += 1
    return planets


def build_planet_from_id(ident='Earthbaseline', initial_kwargs=None, update_kwargs=None, postprocessors=None,
                         t_eval=None,
                         **kwargs):
    planet_kwargs = eval('inputs.' + ident + '_in')
    model_kwargs = eval('inputs.' + ident + '_run')
    if initial_kwargs is not None:
        model_kwargs.update(initial_kwargs)
    if update_kwargs is not None:
        planet_kwargs.update(update_kwargs)
    pl = build_planet(planet_kwargs, model_kwargs, postprocessors, t_eval, **kwargs)
    return pl


def build_planet(planet_kwargs=None, initial_kwargs=None, postprocessors=None, t_eval=None, nondimensional=False,
                 **kwargs):
    if postprocessors is None:
        postprocessors = ['topography']
    if initial_kwargs is None:
        initial_kwargs = {}
    if planet_kwargs is None:
        planet_kwargs = {}
    pl = tp.TerrestrialPlanet(**planet_kwargs)
    pl = thermal.solve(pl, t_eval=t_eval, **initial_kwargs, **kwargs)  # T_m, T_c, D_l

    if 'topography' in postprocessors:
        pl = topography.topography(pl, **kwargs)
    if 'ocean_capacity' in postprocessors:
        pl = oceans.max_ocean(pl, **kwargs)

    if nondimensional:
        pl.nondimensionalise()
    return pl


def build_solarsystem(run_args=None, ident_list=None, dicts=False):
    if ident_list is None:
        ident_list = ['Moon1', 'Mercury1', 'Mars1', 'Venus', 'Earth']
    planets = []
    for ii, ident in enumerate(ident_list):
        pl = build_planet_from_id(ident, run_args)
        planets.append(pl)
    if dicts:
        # Create a zip object from two lists
        z = zip(ident_list, planets)
        # Create a dictionary from zip object
        return dict(z)
    return planets


"                      PLOTTING                           "


def plot_save(fig, fname, fig_path='plat/', fig_fmt='.png', bbox_inches='tight', tight_layout=True, **kwargs):
    path = fig_path + fname + fig_fmt
    directory = os.path.dirname(path)
    os.makedirs(directory, exist_ok=True)
    if tight_layout:
        fig.tight_layout()
    fig.savefig(path, bbox_inches=bbox_inches, **kwargs)
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  saved to ', path, '!')


def plot_output(pl, names, ncols=6, tspan=None, title=None, plots_save=False, write_out=False,
                compare_dir=None, fig_path='plat/', labelpad=None, labelsize=15, legsize=10, fname=None,
                line_args=None, cmp_line_args=None, annotate_colour='xkcd:bright purple',
                print_tf=False, colorbar=False, legend=True, hidex=False, fformat='.png',
                ident=None, fig=None, axes=None, label=None, cmp_label=None, ticksize=12,
                fontname=None, suptitlepad=1.04, legax=0, **kwargs):
    # names: y param
    if ident is None:
        ident = pl.ident
    if fname is None:
        fname = ident
    if line_args is None:
        line_args = {'lw': 2, 'ls': '-', 'c': 'k', 'marker': None, 'ms': 5}
    if cmp_line_args is None:
        cmp_line_args = {'lw': 1, 'ls': '-', 'c': 'r', 'marker': None, 'ms': 5}

    t = pl.t  # time points of ode solutions in s

    nrows = int(np.ceil(len(names) / ncols))
    if (fig is None) and (axes is None):
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3))
    if tspan is None:
        tspan = (0, t[-1] * 1e-9 / parameters.years2sec)
    out_vars = list(names.keys())
    ylabels = list(names.values())  # tuple (ylabel, yscale)
    if label is None:
        label = 'this work'
    for n, par in enumerate(out_vars):
        # loop across axes: each y axis variable
        ax = axes.flatten()[n]
        y = eval('pl.' + par)
        if np.size(y) == 1:
            y = [y] * len(t)
        try:
            #             print('y', y)
            # if the name of the parameter you want to plot exists
            yl = str(ylabels[n][0])
            if (par == 'eta_m') or (par == 'Ra_i') or (par == 'Ra_i_eff') or (
                    par == 'delta_eta'):  # always log scale for viscosity, Ra
                ax.set_yscale('log')  # y = np.log10(y)
            plot_one(ax, t * 1e-9 / parameters.years2sec, y * ylabels[n][1], xlabel='', ylabel=yl, ticksize=ticksize,
                     labelpad=labelpad,
                     label=label, fontname=fontname, labelsize=labelsize, legsize=legsize, line_args=line_args)
            if compare_dir is not None:
                # if data exists to benchmark this param
                try:
                    if (isinstance(compare_dir, collections.Iterable)) and (
                            not isinstance(compare_dir, six.string_types)):
                        for cc, cdir in enumerate(compare_dir):
                            df = pd.read_csv(cdir + '/' + par + '.csv', header=None, names=['time', 'value'],
                                             index_col=False)
                            if cmp_label is None:
                                cmp_label = cdir
                            plot_one(ax, df['time'], df['value'],
                                     '', yl, labelsize=labelsize, legsize=legsize, ticksize=ticksize,
                                     label=cmp_label[cc], fontname=fontname, line_args=cmp_line_args[cc])
                    else:
                        # not iterable
                        df = pd.read_csv(compare_dir + '/' + par + '.csv', header=None, names=['time', 'value'],
                                         index_col=False)
                        if cmp_label is None:
                            cmp_label = compare_dir
                        plot_one(ax, df['time'], df['value'],
                                 '', yl, labelsize=labelsize, legsize=legsize, ticksize=ticksize,
                                 label=cmp_label, fontname=fontname, line_args=cmp_line_args)
                except IOError:
                    print('file', str(par + '.csv'), 'not found')
                    pass
            if par == 'urey' and print_tf:  # print final value of urey ratio
                ii = np.where(t * 1e-9 / parameters.years2sec <= tspan[-1])
                ax.annotate('%.2f' % (y[ii][-1]), xy=(tspan[-1], y[ii][-1]), fontsize=legsize,
                            color=annotate_colour,
                            textcoords="axes fraction", xytext=(0.95, 0.2),
                            ha='right', va='bottom',
                            arrowprops=dict(arrowstyle='->', connectionstyle="arc3,rad=-0.1",
                                            ec=annotate_colour))
            ax.set_xlim(tspan)
            if legend and (n == legax):
                ax.legend(frameon=False, fontsize=legsize)
        except ValueError as e:
            print('could\'t plot', par)
            print(e)

    # hide unused axes
    while n + 1 < ncols * nrows:
        fig.delaxes(axes.flatten()[n + 1])
        n += 1

    plot_setxlabel(axes, 'Age (Gyr)', 'bottom', fontname=fontname, labelpad=labelpad, labelsize=labelsize)
    if title is None:
        title = pl.ident
    fig.suptitle(title, fontsize=labelsize, y=suptitlepad, fontname=fontname)

    # if colorbar:
    #     fig.subplots_adjust(right=0.8)
    #     cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    #     fig.colorbar(plat, cax=cbar_ax)
    #     plat.set_visible = False

    plt.tight_layout()
    if plots_save:
        plt.savefig(fig_path + fname + fformat, bbox_inches='tight')
    if write_out:
        print('\n        n timesteps =', len(t))
        print('$t_f$ =', t[-1] * parameters.sec2Gyr, 'Gyr')
        print(r'$R_p$ =', '%.2f' % (pl.R_p / parameters.R_E), 'R_E =', '%.2f' % (pl.R_p * 1e-3), 'km')
        #         print(r'$R_c$ =', '%.2f'%(kwargs['CRF']*kwargs['R_p']*1e-3), 'km')
        print('M_c', '%.2e' % pl.M_c, 'kg')
        print(r'$M_{m+lid}$ =', '%.2e' % (pl.M_m), 'kg')
        print(r'$g_{sfc}$ =', '%.2f' % pl.g_sfc, 'm s^-2')
        print(r'$g_{cmb}$ =', '%.2f' % pl.g_cmb, 'm s^-2')
        print(r'$\kappa_m$', '%.6f' % pl.kappa_m, 'm^2 s^-1')
        print(r'CRF =', '%.2f' % pl.CRF)
        print(r'$h_0$ =', '%.2f' % (pl.h_rad_m[0] * 1e12), 'pW kg^-1')
        print(r'$h_{4.5}$ =', '%.2f' % (pl.h_rad_m[-1] * 1e12), 'pW kg^-1')
        #         print(r'$H_0$ =', '%.2f'%(H_rad_m[0] + H_rad_lid[0]), 'TW')
        #         print(r'$H_{4.5}$ =', '%.2f'%(H_rad_m[-1] + H_rad_lid[-1]), 'TW')
        print(r'Urey ratio @ $t_f$ =', '%.2f' % pl.urey[-1])
        print('q_sfc(t=0)', '%.2f' % (pl.q_sfc[0] * 1e3), 'mW m^-3')
    return fig, axes


def snaps(pl, plot_snapshots=None, fig_path=None, plots_save=False, fformat='.png', ident=None, **kwargs):
    if ident is None:
        ident = pl.ident

    t = pl.t  # time points of ode solutions in s
    try:
        n_col = len(plot_snapshots)
    except:
        n_col = 1
    fig2, axes2 = plt.subplots(1, n_col, figsize=(3 * n_col, 5))
    for iax, tx in enumerate(plot_snapshots):  # tx is the time value u want nearest
        ii = min(enumerate(t), key=lambda x: abs(tx - x[1] * parameters.sec2Gyr))[0]
        plot_structure(ax=axes2[iax], t=t[ii], T_m=pl.T_m[ii], T_c=pl.T_c[ii], T_s=pl.T_s,
                       T_l=pl.T_l[ii], R_l=pl.R_l[ii], R_p=pl.R_p, R_c=pl.R_c, h_rad_m=pl.h_rad_m[ii],
                       d_lbl=pl.TBL_c[ii], d_ubl=pl.delta_rh[ii], q_ubl=pl.q_ubl[ii], a0=pl.a0[ii],
                       k_m=pl.k_m, legsize=10, **kwargs)
    plt.tight_layout()
    if plots_save:
        fig2.plot_save(fig_path + pl.ident + '_profiles' + fformat, bbox_inches='tight')
    return fig2, axes2


def plot_structure(ax=None, t=None, T_m=None, T_c=None, R_p=None, R_l=None, R_c=None, T_l=None,
                   T_s=None, h_rad_m=None, d_lbl=None, d_ubl=None, q_ubl=None, a0=None, k_m=None,
                   labelsize=16, legsize=14, Tlid_ini=None, **kwargs):
    """ plot temp structure (for a given time) """
    r_c = np.linspace(0, R_c * 1e-3)
    r_lbl = np.linspace(R_c * 1e-3, (R_c + d_lbl) * (1e-3))
    r_m = np.linspace((R_c + d_lbl) * 1e-3, (R_l - d_ubl) * 1e-3)  # radius for mantle in km
    r_ubl = np.linspace((R_l - d_ubl) * 1e-3, (R_l) * 1e-3)
    r_l = np.linspace(R_l * 1e-3, R_p * 1e-3)  # radius for lid
    T_cond = thermal.sph_conduction(r_l * 1e3, a0=a0, T_l=T_l, R_p=R_p, R_l=R_l, T_s=T_s, k_m=k_m, **kwargs)
    q = thermal.sph_flux(r_l * 1e3, a0=a0, T_l=T_l, T_s=T_s, R_p=R_p, R_l=R_l, k_m=k_m, **kwargs)
    if Tlid_ini == 'linear':
        T_cond = thermal.sph_conduction(r_l * 1e3, a0=0, T_l=T_l, R_p=R_p, R_l=R_l, T_s=T_s, k_m=k_m, **kwargs)
        q = thermal.sph_flux(r_l * 1e3, a0=0, T_l=T_l, T_s=T_s, R_p=R_p, R_l=R_l, k_m=k_m, **kwargs)

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(3, 5))
    ax.axhline(y=R_l * 1e-3, ls='--', lw=1, c='xkcd:bluish purple')
    ax.axhline(y=R_c * 1e-3, ls='--', lw=1, c='xkcd:clay')

    ax.plot(T_cond, r_l, c='xkcd:bluish purple')
    ax.plot(T_l + q_ubl / k_m * (R_l - r_ubl * 1e3), r_ubl, c='xkcd:greenish')
    ax.plot([T_m] * len(r_m), r_m, c='xkcd:greenish')
    ax.plot([T_m, T_c], [(R_c + d_lbl) * 1e-3, R_c * 1e-3], c='xkcd:greenish', marker=None)
    ax.plot([T_c] * len(r_c), r_c, c='xkcd:clay')

    ax.set_ylabel('Radius (km)', fontsize=labelsize)
    ax.set_xlabel('Temperature (K)', fontsize=labelsize)
    ax.set_ylim([0, R_p * 1e-3])
    #     #ax.fill_between(x=x, y0=[0]*len(x), y1=[R_cmb*1e-3]*len(x), color='xkcd:gold') # core
    ax.text(T_cond[-1], 0, 'Core', ha='left', va='bottom', fontsize=legsize, c='xkcd:clay')
    #     #ax.fill_between(x=x, y0=[R_cmb*1e-3]*len(x), y1=[R_p*1e-3]*len(x), color='xkcd:tomato') # mantle
    ax.text(T_cond[-1], R_c * 1e-3, 'Convecting region', ha='left', va='bottom', fontsize=legsize, c='xkcd:greenish')
    ax.text(T_cond[-1], R_l * 1e-3, 'Lid', ha='left', va='bottom', fontsize=legsize, c='xkcd:bluish purple')

    ax2 = ax.twiny()
    ax2.set_xlabel('Flux, steady-state (mW m$^{-2}$)', color='xkcd:grey')
    ax2.plot(q * 1e3, r_l, color='xkcd:grey')
    ax2.plot(q_ubl * 1e3, r_ubl[0], marker='o', color='xkcd:grey')
    ax2.annotate('$q_{ubl}$', (q_ubl * 1e3, r_ubl[-1]), color='xkcd:grey', fontsize=12, ha="left", va="top")
    ax2.tick_params(axis='x', labelcolor='xkcd:grey')

    ax.set_title(('%.1f' % (t * 1e-9 / parameters.years2sec)) + ' Gyr', fontsize=labelsize)

    return ax


def interp_benchmark(path, yscale=1):
    df = pd.read_csv(path, header=None, names=['time', 'value'], index_col=False)
    f = interpolate.interp1d(np.array(df['time']), np.array(df['value']) * yscale, kind='linear')
    times = df['time']  # in Gyr
    return np.array(times), f


def plot_qsfc_error(pl, ax3=None, compare_dir=None, fig_path=None, plots_save=False, ident=None, fformat='.png',
                    **kwargs):
    """ sanity check on q_sfc """
    if ident is None:
        ident = pl.ident

    t = pl.t  # time points of ode solutions in s
    if ax3 is None:
        fig3, ax3 = plt.subplots(1, 1, figsize=(5, 5))

    t_D_l, f_D_l_interp = interp_benchmark(path=compare_dir + '/D_l.csv', yscale=1e3)  # in Gyr, m
    temp = t * 1e-9 / parameters.years2sec  # in Gyr

    try:
        t_T_l, f_T_l_interp = interp_benchmark(path=compare_dir + '/T_l.csv')  # in Gyr, K
        iii = np.where((temp >= t_T_l.min()) & (temp <= t_T_l.max()))
        times0 = temp[iii]  # time points of ODE solver subset to interpolation range
        T_l_interp = f_T_l_interp(times0)  # published plot interpolated to model times, in K
    except FileNotFoundError as e:
        t_T_l, f_T_avg_interp = interp_benchmark(path=compare_dir + '/T_avg.csv')  # in Gyr, K
        iii = np.where((temp >= t_T_l.min()) & (temp <= t_T_l.max()))
        times0 = temp[iii]  # time points of ODE solver subset to interpolation range
        T_avg_interp = f_T_avg_interp(times0)  # published plot interpolated to model times, in K
        T_l_interp = Tl_from_Tmean(R_l=pl.R_l[iii], T_avg=T_avg_interp, a0=pl.a0[iii], **kwargs)

    D_l_interp = f_D_l_interp(times0)  # published plot interpolated to model times, in m
    R_l_interp = pl.R_p - D_l_interp
    q_sfc_interp = thermal.sph_flux(pl.R_p, a0=pl.a0[iii], T_l=T_l_interp, T_s=pl.T_s, R_l=R_l_interp,
                                    R_p=pl.R_p, k_m=pl.k_m, **kwargs)  # sfc flux in W m^-2

    ax3.plot(times0, pl.q_sfc[iii] * 1e3, c='xkcd:black', label='this work')
    ax3.plot(times0, q_sfc_interp * 1e3, c='xkcd:blue', label='Thiriet interp')
    df = pd.read_csv(compare_dir + '/q_sfc.csv', header=None, names=['time', 'value'],
                     index_col=False)  # in Gyr, mW m^-3
    ax3.plot(df['time'], df['value'], c='xkcd:red', label='Thiriet digitised')
    ax3.legend(frameon=False, fontsize=14)
    ax3.set_xlabel('Time (Gyr)', fontsize=16)
    ax3.set_ylabel('$q_{sfc}$ (mW m$^{-2}$)', fontsize=16)

    #     plt.tight_layout()

    fig0, ax0 = plt.subplots(1, 1, figsize=(4, 4))
    ax0.plot(times0, pl.q_sfc[iii] * 1e3 - q_sfc_interp * 1e3, c='xkcd:grey')
    ax0.set_xlabel('Time (Gyr)', fontsize=14)
    ax0.set_ylabel('$\Delta q_{sfc}$ (mW m$^{-2}$)', fontsize=14)
    ax0.set_title(
        'Mean error: $\pm$' + '%.2f' % np.mean(np.absolute(pl.q_sfc[iii] * 1e3 - q_sfc_interp * 1e3)) + ' mW m$^{-2}$',
        fontsize=14)
    plt.tight_layout()
    if plots_save:
        fig3.plot_save(fig_path + pl.ident + '_test_qsfc' + fformat)
        fig0.plot_save(fig_path + pl.ident + '_q_error' + fformat)


def plot_Tavg(pl, ax3=None, compare_dir=None, fig_path=None, plots_save=False, ident=None, fformat='.png', **kwargs):
    """ sanity check on T_avg """
    if ident is None:
        ident = pl.ident

    t = pl.t  # time points of ode solutions in s

    if ax3 is None:
        fig3, ax3 = plt.subplots(1, 1, figsize=(5, 5))

    # plot your T_avg calculation but using interpolated published D_l    
    t_D_l, f_D_l_interp = interp_benchmark(path=compare_dir + '/D_l.csv', yscale=1e3)  # in Gyr, m

    # select model time points in interpolation range
    temp = t * 1e-9 / parameters.years2sec  # in Gyr
    iii = np.where((temp >= t_D_l.min()) & (temp <= t_D_l.max()))
    times0 = temp[iii]  # time points of ODE solver subset to interpolation range
    D_l_interp = f_D_l_interp(times0)  # published D_l at model time points in m
    R_l_interp = pl.R_p - D_l_interp  # m

    T_avg_interp = thermal.T_mean(T_m=pl.T_m[iii], T_l=pl.T_l[iii], R_p=pl.R_p, R_l=R_l_interp,
                                  R_c=pl.R_c, a0=pl.a0[iii], T_s=pl.T_s, k_m=pl.k_m, **kwargs)
    ax3.plot(times0, pl.T_avg[iii], c='xkcd:black', label='this work')
    ax3.plot(times0, T_avg_interp, c='xkcd:blue', label='this work with Thiriet D_l')
    df = pd.read_csv(compare_dir + '/T_avg.csv', header=None, names=['time', 'value'],
                     index_col=False)  # in Gyr, mW m^-3
    ax3.plot(df['time'], df['value'], c='xkcd:red', label='Thiriet digitised')
    ax3.legend(frameon=False, fontsize=14)
    ax3.set_xlabel('Time (Gyr)', fontsize=16)
    ax3.set_ylabel('$T_{avg}$ (K)', fontsize=16)
    ax3.set_title('blue should match red')

    plt.tight_layout()
    if plots_save:
        fig3.plot_save(fig_path + pl.ident + '_test_Tavg' + fformat)


def plot_one(ax, x, y, xlabel, ylabel, labelsize=15, legsize=16, ticksize=12, line_args=None,
             text=None, xticks=True, ylim=None, label=None, labelpad=None, fontname=None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(3, 3))
    ax.set_xlabel(xlabel, fontsize=labelsize, fontname=fontname, labelpad=labelpad)
    ax.set_ylabel(ylabel, fontsize=labelsize, fontname=fontname, labelpad=labelpad)
    if not xticks:
        ax.set_xticks([])
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.plot(x, y, label=label, **line_args)
    ax.set_xlim(x.min(), x.max())
    ax.tick_params(axis='x', labelsize=ticksize)
    ax.tick_params(axis='y', labelsize=ticksize)
    if text is not None:
        ax.text(0.05, 0.95, text, ha='left', va='top', transform=ax.transAxes, fontsize=legsize)
    return ax


def plot_setxlabel(axes, label, style, labelsize=15, fontname=None, labelpad=None, **kwargs):
    try:
        ny, nx = np.shape(axes)  # assume 2D
        for ii, ax in enumerate(axes):
            for jj, a in enumerate(ax):
                if (style == 'all') or ((style == 'bottom') and (ii == ny - 1)):
                    a.set_xlabel(label, fontsize=labelsize, fontname=fontname)
                if (style == 'bottom') and (ii < ny - 1):
                    a.set_xticks([])
    except ValueError:  # 1D
        ny = 1
        try:
            for jj, a in enumerate(axes):
                a.set_xlabel(label, fontsize=labelsize, fontname=fontname, labelpad=labelpad)
        except TypeError:  # single
            axes.set_xlabel(label, fontsize=labelsize, fontname=fontname, labelpad=labelpad)


def Tl_from_Tmean(R_c=None, R_l=None, R_p=None, T_avg=None, T_s=None, a0=None, k_m=None, Ea=None, R_b=None, a_rh=None,
                  **kwargs):
    """solved for T_l using sympy"""
    return (k_m * (4 * R_l ** 4 - 4 * R_l ** 3 * R_p - 3 * R_l ** 2 * R_p ** 2 + 2 * R_l * R_p ** 3 + R_p ** 4) * (
            -60 * Ea * R_c ** 6 * k_m + 60 * Ea * R_c ** 3 * R_l ** 3 * k_m + 30 * Ea * R_c ** 3 * R_l ** 2 * R_p * k_m + 30 * Ea * R_c ** 3 * R_l * R_p ** 2 * k_m - 30 * Ea * R_l ** 5 * R_p * k_m - 30 * Ea * R_l ** 4 * R_p ** 2 * k_m + 120 * R_b * R_c ** 3 * R_l ** 3 * T_avg * a_rh * k_m - 60 * R_b * R_c ** 3 * R_l ** 2 * R_p * T_avg * a_rh * k_m - 60 * R_b * R_c ** 3 * R_l * R_p ** 2 * T_avg * a_rh * k_m - 8 * R_b * R_l ** 8 * a0 * a_rh + 14 * R_b * R_l ** 7 * R_p * a0 * a_rh + 9 * R_b * R_l ** 6 * R_p ** 2 * a0 * a_rh - 20 * R_b * R_l ** 5 * R_p ** 3 * a0 * a_rh - 60 * R_b * R_l ** 5 * R_p * T_s * a_rh * k_m - 10 * R_b * R_l ** 4 * R_p ** 4 * a0 * a_rh - 30 * R_b * R_l ** 4 * R_p ** 2 * T_s * a_rh * k_m + 18 * R_b * R_l ** 3 * R_p ** 5 * a0 * a_rh - 120 * R_b * R_l ** 3 * R_p ** 3 * T_avg * a_rh * k_m + 180 * R_b * R_l ** 3 * R_p ** 3 * T_s * a_rh * k_m + R_b * R_l ** 2 * R_p ** 6 * a0 * a_rh + 60 * R_b * R_l ** 2 * R_p ** 4 * T_avg * a_rh * k_m - 30 * R_b * R_l ** 2 * R_p ** 4 * T_s * a_rh * k_m - 4 * R_b * R_l * R_p ** 7 * a0 * a_rh + 60 * R_b * R_l * R_p ** 5 * T_avg * a_rh * k_m - 60 * R_b * R_l * R_p ** 5 * T_s * a_rh * k_m) - 2 * np.sqrt(
        15) * np.sqrt(Ea * k_m ** 3 * (
            60 * Ea * R_c ** 6 * k_m - 60 * Ea * R_c ** 3 * R_l ** 2 * R_p * k_m - 60 * Ea * R_c ** 3 * R_l * R_p ** 2 * k_m + 15 * Ea * R_l ** 4 * R_p ** 2 * k_m + 30 * Ea * R_l ** 3 * R_p ** 3 * k_m + 15 * Ea * R_l ** 2 * R_p ** 4 * k_m - 240 * R_b * R_c ** 3 * R_l ** 3 * T_avg * a_rh * k_m + 120 * R_b * R_c ** 3 * R_l ** 2 * R_p * T_avg * a_rh * k_m + 120 * R_b * R_c ** 3 * R_l * R_p ** 2 * T_avg * a_rh * k_m + 16 * R_b * R_l ** 8 * a0 * a_rh - 28 * R_b * R_l ** 7 * R_p * a0 * a_rh - 18 * R_b * R_l ** 6 * R_p ** 2 * a0 * a_rh + 40 * R_b * R_l ** 5 * R_p ** 3 * a0 * a_rh + 120 * R_b * R_l ** 5 * R_p * T_s * a_rh * k_m + 20 * R_b * R_l ** 4 * R_p ** 4 * a0 * a_rh + 60 * R_b * R_l ** 4 * R_p ** 2 * T_s * a_rh * k_m - 36 * R_b * R_l ** 3 * R_p ** 5 * a0 * a_rh + 240 * R_b * R_l ** 3 * R_p ** 3 * T_avg * a_rh * k_m - 360 * R_b * R_l ** 3 * R_p ** 3 * T_s * a_rh * k_m - 2 * R_b * R_l ** 2 * R_p ** 6 * a0 * a_rh - 120 * R_b * R_l ** 2 * R_p ** 4 * T_avg * a_rh * k_m + 60 * R_b * R_l ** 2 * R_p ** 4 * T_s * a_rh * k_m + 8 * R_b * R_l * R_p ** 7 * a0 * a_rh - 120 * R_b * R_l * R_p ** 5 * T_avg * a_rh * k_m + 120 * R_b * R_l * R_p ** 5 * T_s * a_rh * k_m)) * (
                    R_c - R_l) * (R_l - R_p) ** 2 * (2 * R_l + R_p) ** 2 * (R_c ** 2 + R_c * R_l + R_l ** 2)) / (
                   30 * R_b * R_l ** 2 * a_rh * k_m ** 2 * (R_l - R_p) ** 2 * (2 * R_l + R_p) ** 2 * (
                   4 * R_l ** 4 - 4 * R_l ** 3 * R_p - 3 * R_l ** 2 * R_p ** 2 + 2 * R_l * R_p ** 3 + R_p ** 4))


def eta_from_Ra(rho=None, g=None, alpha=None, dT=None, d=None, kappa=None, Ra=None):
    return rho * g * alpha * dT * d ** 3 / (kappa * Ra)


def Ra_from_RaB(F=None, dT=1000, l=None, k=None, Ra_B=None):
    # convert between basal heating Ra and internal Ra, given F flux into base
    return Ra_B * k * dT / (F * l)  # basal heating Ra_B


def Ra_from_RaF_2(F=74e-3, dT=1000, l=750e3, kappa=8e-7, rho=3300, c_p=1200, Ra_F=2.4e6):  # moresi and parsons
    return Ra_F * (rho * c_p * kappa * dT) / (l * F)


def Ra_from_RaF(F=None, dT_m=None, k=None, l=None, Ra_F=None, **kwargs):  # F is surface flux
    return Ra_F / (l * F / (k * dT_m))


def benchmark_thermal_plots(ident, show_qsfc_error=False, show_Tavg=False, names=None, pl_update_args=None,
                            model_update_args=None, verbose=False, **kwargs):
    if names is None:
        names = {'T_avg': ('$T_{avg}$ (K)', 1),
                 'q_sfc': ('$q_{sfc}$ (mW m$^{-2}$)', 1e3),
                 'D_l': ('$D_l$ (km)', 1e-3),
                 'q_core': ('$q_{B}$ (mW m$^{-2}$)', 1e3),
                 'urey': ('Ur', 1),
                 'T_l': ('$T_l$ (K)', 1),
                 }
    planet_kwargs = eval('inputs.' + ident + '_in')
    model_kwargs = eval('inputs.' + ident + '_run')
    if pl_update_args is not None:
        planet_kwargs.update(pl_update_args)
    if model_update_args is not None:
        model_kwargs.update(model_update_args)
    pl = tp.TerrestrialPlanet(**planet_kwargs)
    pl = thermal.solve(pl, **model_kwargs)  # T_m, T_c, D_l
    pl = topography.topography(pl, C=1)

    if verbose:
        print('T_c', pl.T_c[-1])
        print('d_m', pl.d_m[-1])
        print('dT_m', abs(pl.dT_m[-1]))
        print('g_sfc', pl.g_sfc)
        print('q_ubl', pl.q_ubl[-1])
        print('eta', pl.eta_m[-1])
        print('Ra_i_eff', pl.Ra_i_eff[-1])
        print('h_rms', pl.dyn_top_rms[-1])
        print('h_iso', pl.dyn_top_rms_isoviscous[-1])

    fig, axes = plot_output(pl, names, verbose=False, **kwargs)
    if show_qsfc_error:
        plot_qsfc_error(pl, ident=ident, **kwargs)
    if show_Tavg:
        plot_Tavg(pl, **kwargs)
    return fig, axes


def plot_vs_x(scplanets=None, lplanets=None, xname=None, ynames=None, planets2=None, fig=None, axes=None,
              labels=False, labelsize=15, legsize=12, alpha=1, legend=False, snap=4.5, labelpad=None,
              plots_save=False, s=30, ls='-', lw=1, cmap='rainbow', marker='o', legtitle=None, legendtop=False,
              colorbar=False, c='k', set_ylabel=True, ymin=None, ymax=None, set_ylim=True, set_xlim=False,
              fformat='.png',
              zorder_l=None, zorder_sc=None, label_l=None, fname=None, ticksize=12, xmin=None, xmax=None,
              fig_path='plat/', printrange=False, log=False, relative=False, **kwargs):
    # for a list of planets, plot some parameter on the y axis vs. parameter x
    if (c is None) and (scplanets is not None):
        c = np.arange(len(scplanets))
    #         colour = cm.get_cmap(cmap)
    #         norm = colors.Normalize(vmin=0, vmax=len(planets))
    nax = len(ynames)
    if xmin is not None:
        set_xlim = True
    if axes is None:
        fig, axes = plt.subplots(1, nax, figsize=(5 * nax, 4))

    xparam = list(xname.keys())[0]
    xlabels = list(xname.values())[0]
    yparam = list(ynames.keys())
    ylabels = list(ynames.values())  # tuple (ylabel, yscale)

    ii = 0
    while ii < nax:
        try:
            ax = axes[ii]
        except TypeError:  # single ax
            ax = axes
        if scplanets is not None:
            x = []
            y = []
            for ip, pl in enumerate(scplanets):  # planets to plot as scatter
                data_x = eval('pl.' + xparam) * xlabels[1]
                if isinstance(data_x, Iterable):
                    data_x = data_x[-1]  # get final value
                x.append(data_x)
                data_y = eval('pl.' + yparam[ii]) * ylabels[ii][1]
                if isinstance(data_y, Iterable):
                    data_y = data_y[-1]  # get final value
                y.append(data_y)
                if labels:
                    ax.annotate(xy=(data_x, data_y), s=pl.ident[0:2], fontsize=legsize)
            sc = ax.scatter(x, y, c=c, s=s, marker=marker, cmap=cmap, zorder=zorder_sc)
            if colorbar and (ii == 0):
                plt.colorbar(sc)
        if lplanets is not None:
            x = []
            y = []
            xscale = xlabels[1]
            yscale = ylabels[ii][1]
            try:
                for ip, pl in enumerate(lplanets):  # planets to plot as line
                    t = pl.t
                    it = age_index(t, snap, parameters.sec2Gyr)
                    data_x = eval('pl.' + xparam) * xscale
                    if isinstance(data_x, Iterable):
                        data_x = data_x[it]  # if an evolution model then take certain snap
                    x.append(data_x)

                    if log and relative:
                        # data_y = 10**( np.log10(eval('pl.' + yparam[ii])) * np.log10(yscale) )
                        # print('data_y', np.min(data_y), np.max(data_y))
                        data_y = eval('pl.' + yparam[ii]) * yscale
                        # print('data_y', np.min(data_y), np.max(data_y))
                    else:
                        data_y = eval('pl.' + yparam[ii]) * yscale
                    if isinstance(data_y, Iterable):
                        data_y = data_y[it]  # if an evolution model then take certain snap
                    y.append(data_y)
            except TypeError:  # if given a single planet (not iterable) - get values across evol
                x = eval('lplanets.' + xparam) * xscale
                y = eval('lplanets.' + yparam[ii]) * yscale
            # sort
            if isinstance(y, Iterable) and isinstance(x, Iterable):
                x, y = zip(*sorted(zip(x, y)))
            elif isinstance(y, Iterable) and (not isinstance(x, Iterable)):
                x, y = zip(*sorted(zip([x] * np.ones_like(y), y)))
            elif (not isinstance(y, Iterable)) and isinstance(x, Iterable):
                x, y = zip(*sorted(zip(x, [y] * np.ones_like(x))))

            ax.plot(x, y, ls=ls, c=c, lw=lw, alpha=alpha, zorder=zorder_l, label=label_l)
            if printrange:
                print('range:', np.min(y), np.max(y))
        if set_ylabel:
            ax.set_ylabel(ylabels[ii][0], fontsize=labelsize, labelpad=labelpad)
            ax.tick_params(axis='y', labelsize=ticksize)
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
            # ax.yaxis.set_ticklabels([])
        ax.tick_params(axis='x', labelsize=ticksize)

        # log scale for viscosity
        if (yparam[ii] is 'eta_m') or (yparam[ii] is 'nu_m') or (yparam is 'Ra_i_eff'):
            ax.set_yscale('log')
        if (xparam is 'eta_m') or (xparam is 'nu_m') or (xparam is 'eta_0') or (xparam is 'Ra_i_eff'):
            ax.set_xscale('log')
        if log:
            ax.set_xscale('log')
            ax.set_yscale('log')

        if set_ylim:
            ax.set_ylim(ymin, ymax)
        if set_xlim:
            if xmin is None:
                xmin = np.min(x)
            if xmax is None:
                xmax = np.min(y)
            ax.set_xlim(xmin, xmax)

        if legend:
            if legendtop:
                legend = ax.legend(frameon=False, fontsize=legsize,
                                   borderaxespad=0, title=legtitle,  # mode="expand",
                                   loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=2, )
            else:
                legend = ax.legend(frameon=False, fontsize=legsize, loc='upper left',
                                   bbox_to_anchor=(1.05, 0.9),
                                   borderaxespad=0, ncol=1, title=legtitle)
            if legtitle is not None:
                plt.setp(legend.get_title(), fontsize=legsize)
                legend._legend_box.align = "left"
        ii += 1
    plot_setxlabel(axes, xlabels[0], 'all', labelsize=labelsize, labelpad=labelpad, **kwargs)

    plt.tight_layout()
    if plots_save:
        if fname is None:
            fname = 'scatter_' + xparam
        plt.savefig(fig_path + fname + fformat, bbox_inches='tight')
    return fig, axes


#
# def plot_change_with_observeables(defaults='Earthbaseline', wspace=0.1, tickwidth=1, relative=True, textc='k',
#                                   age=4.5, x_vars=None, ylabel='$\Delta h$ / $\Delta h_0$  ',
#                                   xlabels=None, nplanets=20, log=False,
#                                   fig=None, axes=None, model_param='dyn_top_rms', legend=False, legsize=12, yscale=1,
#                                   pl_baseline=None, update_kwargs={}, initial_kwargs={}, verbose=False, **kwargs):
#     if x_vars is None:
#         x_vars = ['age', 'M_p', 'CMF', 'H0', 'Ea']
#     if xlabels is None:
#         xlabels = x_vars
#     if axes is None:
#         fig, axes = plt.subplots(1, len(x_vars), figsize=(4 * len(x_vars), 4), sharey=True)
#
#     model_baseline = eval('pl_baseline.' + model_param)
#
#     it = age_index(pl_baseline.t, age, parameters.sec2Gyr)
#     model_baseline = model_baseline[it]
#     i_ax = 0
#     if not_iterable(axes):
#         axes = [axes]
#
#     if relative:
#         yscale = model_baseline ** -1
#     print('yscale', yscale)
#
#     ylabel=True
#     legendd=legend
#
#     if 'age' in x_vars:
#         # time/age variation
#         fig, ax = plot_vs_x(legend=legendd, legsize=legsize, log=log,
#                             lplanets=pl_baseline, xname={'t': ('Age (Gyr)', parameters.sec2Gyr)}, set_xlim=True,
#                             ynames={model_param: (ylabel, yscale)}, ylabel=ylabel, relative=relative,
#                             plots_save=False, fig=fig, axes=axes[i_ax], xmin=1.5, xmax=4.5, **kwargs)
#         if relative:
#             ax.axhline(y=1, lw=1, alpha=0.7, zorder=0)
#         if legend and relative:
#             ax.text(0.95, 0.95,
#                     '{:1.0f}'.format(pl_baseline.M_p / parameters.M_E) + ' $M_E$ \n 300 kJ mol$^{-1}$ \n 0.3 CMF \n 4.6 pW kg$^{-1}$',
#                     fontsize=legsize, c=textc,
#                     horizontalalignment='right',
#                     verticalalignment='top',
#                     transform=ax.transAxes)
#         i_ax += 1
#         ylabel = False
#         legendd = False
#
#     if 'M_p' in x_vars:
#         # mass variation
#         print('generating planets...')
#         planets_mass = bulk_planets(n=nplanets, name='M_p', mini=0.1 * parameters.M_E, maxi=6 * parameters.M_E, like=defaults,
#                                     t_eval=pl_baseline.t, random=False,
#                                     initial_kwargs=initial_kwargs, update_kwargs=update_kwargs, **kwargs)
#         print('finished!')
#         fig, ax = plot_vs_x(legend=legendd, set_xlim=True,  legsize=legsize, log=log,
#             lplanets=planets_mass, xname={'M_p': ('$M_p$ ($M_E$)', parameters.M_E ** -1)},
#             ynames={model_param: ('', yscale)}, snap=age, relative=relative,
#             plots_save=False, fig=fig, axes=axes[i_ax], ylabel=ylabel, **kwargs)
#         if relative:
#             ax.axhline(y=1, lw=1, alpha=0.7, zorder=0)
#         if legend and relative:
#             ax.text(0.05, 0.95, '4.5 Ga\n300 kJ mol$^{-1}$\n0.3 CMF\n4.6 pW kg$^{-1}$', fontsize=legsize,
#                     horizontalalignment='left', c=textc,
#                     verticalalignment='top',
#                     transform=ax.transAxes)
#         i_ax += 1
#         ylabel = False
#         legendd = False
#
#     if 'CMF' in x_vars:
#         # CMF variation
#         planets_CMF = bulk_planets(n=nplanets, name='CMF', mini=0.1, maxi=0.7, like=defaults, t_eval=pl_baseline.t,
#                                    random=False, initial_kwargs=initial_kwargs, update_kwargs=update_kwargs, **kwargs)
#         # where did mini=0.07829, maxi=0.544, come from?
#         fig, ax = plot_vs_x(legend=legendd, set_xlim=True,  legsize=legsize, log=log,
#             lplanets=planets_CMF, xname={'CMF': ('Core Mass Fraction', 1)},
#             ynames={model_param: ('', yscale)}, snap=age, relative=relative,
#             plots_save=False, fig=fig, axes=axes[i_ax], ylabel=ylabel, **kwargs)
#         if relative:
#             ax.axhline(y=1, lw=1, alpha=0.7, zorder=0)
#         if legend and relative:
#             ax.text(0.95, 0.95,
#                     '4.5 Ga \n' + '{:1.0f}'.format(pl_baseline.M_p / parameters.M_E) + ' $M_E$ \n' + '300 kJ mol$^{-1}$  \n 4.6 pW kg$^{-1}$',
#                     fontsize=legsize, c=textc,
#                     horizontalalignment='right',
#                     verticalalignment='top',
#                     transform=ax.transAxes)
#         i_ax += 1
#         ylabel = False
#         legendd = False
#
#     if 'H0' in x_vars:
#         # H0 variation
#         planets_H0 = bulk_planets(n=nplanets, name='H_0', mini=10e-12, maxi=40e-12, like=defaults, t_eval=pl_baseline.t,
#                                   random=False, initial_kwargs=initial_kwargs, update_kwargs=update_kwargs, **kwargs)
#         fig, ax = plot_vs_x(legend=legendd, xmin=10, xmax=40, set_xlim=True,  legsize=legsize,
#                             lplanets=planets_H0, xname={'H_0': ('$H_0$ (pW kg$^{-1}$)', 1e12)},
#                             ynames={model_param: ('', yscale)}, snap=age, log=log, relative=relative,
#                             plots_save=False, fig=fig, axes=axes[i_ax], ylabel=ylabel, **kwargs)
#         if relative:
#             ax.axhline(y=1, lw=1, alpha=0.7, zorder=0)
#         if legend and relative:
#             ax.text(0.95, 0.95, '4.5 Ga \n ' +
#                     '{:1.0f}'.format(pl_baseline.M_p / parameters.M_E) + ' $M_E$ \n 300 kJ mol$^{-1}$  \n 0.3 CMF',
#                     fontsize=legsize,
#                     horizontalalignment='right',
#                     verticalalignment='top', c=textc,
#                     transform=ax.transAxes)
#         i_ax += 1
#         ylabel = False
#         legendd = False
#
#     if 'Ea' in x_vars:
#         # Ea variation
#         planets_Ea = bulk_planets(n=nplanets, name='Ea', mini=250e3, maxi=350e3, like=defaults, t_eval=pl_baseline.t,
#                                   random=False, initial_kwargs=initial_kwargs, update_kwargs=update_kwargs, **kwargs)
#         fig, ax = plot_vs_x(legend=legendd, xmin=250, xmax=350, set_xlim=True,  legsize=legsize,
#                             lplanets=planets_Ea, xname={'Ea': ('$E_a$ (kJ mol$^{-1}$)', 1e-3)},
#                             ynames={model_param: ('', yscale)}, snap=age, log=log, relative=relative,
#                             plots_save=False, fig=fig, axes=axes[i_ax], ylabel=ylabel, **kwargs)
#         if relative:
#             ax.axhline(y=1, lw=1, alpha=0.7, zorder=0)
#         if legend and relative:
#             ax.text(0.95, 0.95, '4.5 Ga \n ' +
#                     '{:1.0f}'.format(pl_baseline.M_p / parameters.M_E) + ' $M_E$ \n 0.3 CMF \n 4.6 pW kg$^{-1}$',
#                     fontsize=legsize, c=textc,
#                     horizontalalignment='right',
#                     verticalalignment='top',
#                     transform=ax.transAxes)
#
#     for ax in axes:
#         ax.xaxis.set_tick_params(width=tickwidth)
#         ax.yaxis.set_tick_params(width=tickwidth)
#     plt.subplots_adjust(wspace=wspace)
#     return fig, axes


def plot_change_with_observeables(defaults='Earthbaseline', wspace=0.1, tickwidth=1, relative=True, textc='k',
                                  age=4.5, x_vars=None, ylabel='$\Delta h$ / $\Delta h_0$  ', fig_height=4,
                                  xlabels=None, nplanets=20, log=False, x_range=None, xscales=None, units=None,
                                  fig=None, axes=None, model_param='dyn_top_rms', legend=False, legsize=12, yscale=1,
                                  pl_baseline=None, update_kwargs={}, initial_kwargs={}, verbose=False, **kwargs):
    if x_vars is None:
        x_vars = ['t', 'M_p', 'CMF', 'H_0', 'Ea']
    if units is None:
        units = ['Gyr', '$M_E$', 'CMF', 'pW kg$^{-1}$', 'kJ mol$^{-1}$']
    if xlabels is None:
        xlabels = x_vars
    if x_range is None:
        x_range = [(1.5, age), (0.1 * parameters.M_E, 6 * parameters.M_E), (0.1, 0.7), (10e-12, 40e-12), (250e3, 350e3)]
    if xscales is None:
        xscales = [parameters.sec2Gyr, parameters.M_E ** -1, 1, 1e12, 1e-3]
    if axes is None:
        fig, axes = plt.subplots(1, len(x_vars), figsize=(4 * len(x_vars), fig_height), sharey=True)
    if not_iterable(axes):
        axes = [axes]

    it = age_index(pl_baseline.t, age, parameters.sec2Gyr)
    model_baseline = eval('pl_baseline.' + model_param)[it]
    if relative:
        yscale = model_baseline ** -1
    set_ylabel = True
    legendd = legend
    for i_ax, x_var in enumerate(x_vars):
        xmin, xmax = x_range[i_ax]
        if x_var == 't':
            # time/age variation - plot single planet evol
            fig, ax = plot_vs_x(lplanets=pl_baseline, xname={'t': (xlabels[i_ax], xscales[i_ax])},
                                ynames={model_param: (ylabel, yscale)}, fig=fig, axes=axes[i_ax], legsize=legsize,
                                legend=legendd, plots_save=False, set_ylabel=set_ylabel, set_xlim=True, xmin=xmin,
                                xmax=xmax, log=log, relative=relative, **kwargs)
        else:
            if verbose:
                print('generating planets across', x_var, '...')
            planets = bulk_planets(n=nplanets, name=x_var, mini=xmin, maxi=xmax,
                                   like=defaults,
                                   t_eval=pl_baseline.t, random=False, verbose=verbose,
                                   initial_kwargs=initial_kwargs, update_kwargs=update_kwargs, **kwargs)
            fig, ax = plot_vs_x(lplanets=planets, xname={x_var: (xlabels[i_ax], xscales[i_ax])},
                                ynames={model_param: ('', yscale)}, fig=fig, axes=axes[i_ax], legsize=legsize,
                                legend=legendd, snap=age, plots_save=False, set_ylabel=ylabel, set_xlim=True, log=log,
                                relative=relative, **kwargs)

        if relative:
            ax.axhline(y=1, lw=1, alpha=0.7, zorder=0)
        if legend and relative:
            string = ''
            for jj, u in enumerate(units):
                if jj != i_ax:
                    if x_vars[jj] == 't':
                        string = string + '{:.3g} '.format(age * xscales[jj]) + u + '\n'
                    else:
                        string = string + '{:.3g} '.format(eval('pl_baseline.' + x_vars[jj]) * xscales[jj]) + u + '\n'
            ax.text(0.95, 0.95, string,
                    fontsize=legsize, c=textc,
                    horizontalalignment='right',
                    verticalalignment='top',
                    transform=ax.transAxes)

        # turn off
        set_ylabel = False
        legendd = False

    for ax in axes:
        ax.xaxis.set_tick_params(width=tickwidth)
        ax.yaxis.set_tick_params(width=tickwidth)
    plt.subplots_adjust(wspace=wspace)
    return fig, axes


def plot_h_relative_multi(defaults='Earthbaseline', save=False, fname='relative_h',
                          models=None, labels=None, c=None, fig=None, axes=None, age=4.5,
                          initial_kwargs={}, update_kwargs={}, legend=True,
                          ylabel='$\Delta h$ / $\Delta h_0$  ', **kwargs):
    initial_kwargs.update({'tf': age})
    pl_baseline = build_planet_from_id(ident=defaults,
                               initial_kwargs=initial_kwargs, update_kwargs=update_kwargs,
                               postprocessors=['topography'], t_eval=None)
    legendd = False
    for ii, h_param in enumerate(models):
        if ii == len(models) - 1:
            legendd = True
        fig, axes = plot_change_with_observeables(defaults=defaults, fig=fig, axes=axes, model_param=models[ii],
                                                  legend=legendd, pl_baseline=pl_baseline, label_l=labels[ii], c=c[ii],
                                                  age=age, ylabel=ylabel, initial_kwargs=initial_kwargs,
                                                  update_kwargs=update_kwargs, **kwargs)

    if save:
        plot_save(fig, fname, **kwargs)

    return fig, axes


def plot_ocean_capacity_relative(age=4.5, legsize=16, fname='ocean_vol', mass_frac_sfcwater=None, textc='k', M0=0.815,
                                 titlesize=24, save=False, spectrum_fname='', spectrum_fpath='', c='#81f79f', title='',
                                 ticksize=14, labelsize=16, clabel='Surface water mass fraction', clabelpad=20,
                                 relative=False,
                                 mass_iax=0, leg_bbox=(1.7, 1.01), log=False, figsize=(10, 10), ytitle=1.1,
                                 cmap='terrain_r',
                                 defaults='Venusbaseline', ylabel=r'$V_{\mathrm{max}}/V_{\mathrm{max, Ve}}$', **kwargs):
    phi0, degree = harm.load_spectrum(fpath=spectrum_fpath, fname=spectrum_fname)
    h_rms0 = harm.powerspectrum_RMS(power_lm=phi0, degree=degree)
    pl0 = bulk_planets(n=1, name='M_p', mini=M0 * parameters.M_E, maxi=M0 * parameters.M_E, like=defaults, t_eval=None,
                       random=False, phi0=phi0, h_rms0=h_rms0, postprocessors=['topography', 'ocean_capacity'],
                       **kwargs)[0]
    fig, axes = plt.subplots(figsize=figsize)
    fig, axes = plot_change_with_observeables(defaults=defaults, model_param='max_ocean', legend=True, pl_baseline=pl0,
                                              textc=textc,
                                              label_l=None, c=c, ylabel=ylabel, age=age, h_rms0=h_rms0, legsize=legsize,
                                              postprocessors=['topography', 'ocean_capacity'], phi0=phi0, log=log,
                                              fig=fig,
                                              axes=axes, ticksize=ticksize, labelsize=labelsize, relative=relative,
                                              **kwargs)

    if mass_frac_sfcwater is not None:
        # how does actual vol scale assuming constant mass fraction of surface water (bad assumption)?
        ax = axes[mass_iax]
        rho_w = 1000
        vol_0 = pl0.max_ocean[-1]
        print('vol_0', vol_0)
        masses = np.logspace(np.log10(0.1), np.log10(6))  # mass in M_E
        colours = colorize([np.log10(m) for m in mass_frac_sfcwater], cmap=cmap)[0]

        for ii, X in enumerate(mass_frac_sfcwater):
            M_w = masses * parameters.M_E * X  # mass of sfc water in kg
            vol_w = M_w / rho_w  # volume of surface water if same mass fraction
            if relative:
                vol_w = vol_w / vol_0
            ax.plot(masses, vol_w, alpha=0.4, lw=0, zorder=0, c=colours[ii], marker='o', markersize=15,
                    label='Maximum water budget')

        colourbar(mappable=None, ax=ax, vmin=np.min(mass_frac_sfcwater), vmax=np.max(mass_frac_sfcwater), label=clabel,
                  labelsize=labelsize * 0.8, ticksize=ticksize, labelpad=clabelpad, ticks=mass_frac_sfcwater,
                  cmap=cmap, c=textc, log=True, pad=0.2)
        #
        # # title and legend
        # legend = ax.legend(frameon=False, fontsize=legsize,
        #                         borderaxespad=0,  # mode="expand",
        #                         loc='lower left', bbox_to_anchor=leg_bbox, ncol=1)

    fig.suptitle(title, fontsize=titlesize, y=ytitle, c=textc)  # x=0.365,

    if save:
        plot_save(fig, fname, **kwargs)
    return fig, axes


def read_JFR(fname='', path='benchmarks/JFR/'):
    df = pd.read_csv(path + fname, header=0, index_col=False)
    Ra = np.array(df.Ra)
    h_peak = np.array(df.peak_topo)
    h_rms = np.array(df.RMS_topo)
    Nu = np.array(df.Nu)
    return Ra, h_peak, h_rms, Nu
