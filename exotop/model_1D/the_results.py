""" Run 1D model and plot results """

import numpy as np
import matplotlib.pyplot as plt
import os
from model_1D import evolve as evol
from model_1D import inputs
from model_1D import terrestrialplanet as tp
from model_1D import parameters
from model_1D import thermal
from model_1D import topography
from model_1D import oceans
from collections.abc import Iterable
import collections
import six
import pandas as pd
from scipy import interpolate
import sh_things as sh
import matplotlib.animation as animation
from useful_and_bespoke import age_index, dark_background, not_iterable, colorize, colourbar
from model_1D.parameters import M_E
import matplotlib.ticker as ticker
import matplotlib.lines as mlines

"                      PLOTTING                           "


def plot_save(fig, fname, fig_path='plat/', fig_fmt='.png', bbox_inches='tight', tight_layout=True, **kwargs):
    path = fig_path + fname + fig_fmt
    if fig_path != '':
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
                     label=label, fontname=fontname, labelsize=labelsize, legsize=legsize, line_args=line_args,
                     **kwargs)
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
                                     label=cmp_label[cc], fontname=fontname, line_args=cmp_line_args[cc], **kwargs)
                    else:
                        # not iterable
                        df = pd.read_csv(compare_dir + '/' + par + '.csv', header=None, names=['time', 'value'],
                                         index_col=False)
                        if cmp_label is None:
                            cmp_label = compare_dir
                        plot_one(ax, df['time'], df['value'],
                                 '', yl, labelsize=labelsize, legsize=legsize, ticksize=ticksize,
                                 label=cmp_label, fontname=fontname, line_args=cmp_line_args, **kwargs)
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

            if legend and (n == legax):
                ax.legend(frameon=False, fontsize=legsize)
        except ValueError as e:
            print('could\'t plot', par)
            print(e)

    # hide unused axes
    try:
        while n + 1 < ncols * nrows:
            fig.delaxes(axes.flatten()[n + 1])
            n += 1
    except:
        pass

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
             text=None, xticks=True, xlim=None, ylim=None, label=None, labelpad=None, fontname=None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(3, 3))
    if xlim is None:
        xlim = (x.min(), x.max())
    ax.set_xlabel(xlabel, fontsize=labelsize, fontname=fontname, labelpad=labelpad)
    ax.set_ylabel(ylabel, fontsize=labelsize, fontname=fontname, labelpad=labelpad)
    if not xticks:
        ax.set_xticks([])
    mask = (x >= xlim[0]) & (x <= xlim[1])
    try:
        ax.plot(x[mask], y[mask], label=label, **line_args)
    except TypeError:
        ax.plot(x[mask], [y] * len(x[mask]), label=label, **line_args)
    ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        ax.set_ylim(None, None, auto=True)
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
    pl = evol.solve(pl, verbose=verbose, run_kwargs=model_kwargs)  # T_m, T_c, D_l
    pl = topography.topography(pl, C=1)
    pl.nondimensionalise()

    if verbose:
        print('\nT_c', pl.T_c[-1])
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


def ensemble_marginal_distribution(yvar, xvar, default='baseline', dist_res=100, update_kwargs=None, run_kwargs=None,
                                   yscale=1, age=4.5, x_res=8, minx=None, maxx=None, verbose=False,
                                   names=None, mini=None, maxi=None, t_eval=None,
                                   n_sigma=1, log=False, **kwargs):
    """generate ensemble of planets depending on x independent variable over some random variations of other
     parameters -- do this for """
    if maxi is None:
        maxi = [300e3, 2.5e12, 2000, 200e3]
    if mini is None:
        mini = [240e3, 1.5e10, 1000, 50e3]
    if names is None:
        names = ['Ea', 'eta_pre', 'T_m0', 'D_l0']

    if default is not None:
        pl_kwargs_base = eval('inputs.' + default + '_in').copy()
        model_kwargs_base = eval('inputs.' + default + '_run').copy()
    else:
        pl_kwargs_base = {}  # use defaults given in terrestrialplanet.py
        model_kwargs_base = {}  # initial conditions defaults given in thermal.py
    if update_kwargs is not None:
        pl_kwargs_base.update(update_kwargs)
    if run_kwargs is not None:
        model_kwargs_base.update(run_kwargs)

    grid = np.zeros((x_res, dist_res))
    if log:
        x_vec = np.logspace(np.log10(minx), np.log10(maxx), x_res)
    else:
        x_vec = np.linspace(minx, maxx, x_res)
    for ii, x in enumerate(x_vec):
        # reset
        pl_kwargs = pl_kwargs_base.copy()
        pl_kwargs.update({xvar: x})
        model_kwargs = model_kwargs_base.copy()
        if ii == 1 and t_eval is None:
            t_eval = pl_ensemble[0].t

        pl_ensemble = evol.bulk_planets_mc(n=dist_res, names=names, mini=mini, maxi=maxi, pl_kwargs=pl_kwargs,
                                           model_kwargs=model_kwargs, t_eval=t_eval, log=False, **kwargs)

        it = age_index(pl_ensemble[0].t, age, age_scale=parameters.sec2Gyr)  # should all be evaluated at same time
        col = np.array([vars(pl)[yvar][it] for pl in pl_ensemble]) * yscale
        grid[ii, :] = col
        if verbose:
            print('    generated ensemble for', x)

    # mean, std
    y_av = np.mean(grid, axis=1)
    y_std = np.std(grid, axis=1)
    y_upper = y_av + y_std * n_sigma  # todo for log scape
    y_lower = y_av - y_std * n_sigma

    return x_vec, y_av, y_upper, y_lower


def ensemble_time_distribution(yvar, default='baseline',
                               dist_res=100, update_kwargs=None, run_kwargs=None, yscale=1,
                               names=['Ea', 'eta_pre', 'T_m0', 'D_l0'],
                               mini=[240e3, 1.5e10, 1000, 100e3],
                               maxi=[300e3, 2.5e12, 2000, 300e3], min_t=0, max_t=4.5,
                               n_sigma=1, log=False, **kwargs):
    # generate ensemble of planets depending on x independent variable over some random variations of other parameters and plot evol
    if default is not None:
        pl_kwargs = eval('inputs.' + default + '_in').copy()
        model_kwargs = eval('inputs.' + default + '_run').copy()
    else:
        pl_kwargs = {}  # use defaults given in terrestrialplanet.py
        model_kwargs = {}  # initial conditions defaults given in thermal.py
    if update_kwargs is not None:
        pl_kwargs.update(update_kwargs)
    if run_kwargs is not None:
        model_kwargs.update(run_kwargs)
    model_kwargs.update({'tf': max_t})

    pl_ensemble = evol.bulk_planets_mc(n=dist_res, names=names, mini=mini, maxi=maxi, pl_kwargs=pl_kwargs,
                                       model_kwargs=model_kwargs, **kwargs)

    # mean, std
    y_all = []
    for pl in pl_ensemble:
        t = pl.t
        y = eval('pl.' + yvar) * yscale
        y_all.append(y)

    y_all = np.array(y_all)
    y_av = np.mean(y_all, axis=0)
    y_std = np.std(y_all, axis=0)
    y_upper = y_av + y_std * n_sigma  # todo for log scape
    y_lower = y_av - y_std * n_sigma

    return t, y_av, y_upper, y_lower


def plot_distribution(yvars, default='baseline',
                      # xmin=0.1*M_E, xmax=6*M_E, logx=True, xres=100,
                      num=100, update_kwargs=None, run_kwargs=None,
                      names=['Ea', 'eta_pre', 'T_m0', 'T_c0', 'D_l0'],
                      mini=[240e3, 1.5e10, 1000, 2000, 100e3],
                      maxi=[300e3, 2.5e12, 2000, 2500, 300e3],
                      xlabelpad=None, ylabelpad=None, n_sigma=1, ylims=None, tickpad=10,
                      fig=None, axes=None, c='k', lw=0.5, alpha=0.7, c_mean='k', log=None, xticks=None, yticks=None,
                      xlabel='Time (Gyr)', ylabels=None, yscales=None, labelsize=16, ticksize=12, save=False,
                      fname='evol_dist', fig_path='', legtext=None, legsize=16, **kwargs):
    # generate ensemble of planets depending on x independent variable over some random variations of other parameters and plot evol
    if yscales is None:
        yscales = [1] * len(yvars)
    if ylabels is None:
        ylabels = [''] * len(yvars)
    if log is None:
        log = [False] * len(yvars)
    if ylims is None:
        ylims = [None] * len(yvars)
    if default is not None:
        pl_kwargs = eval('inputs.' + default + '_in')
        model_kwargs = eval('inputs.' + default + '_run')
    else:
        pl_kwargs = {}  # use defaults given in terrestrialplanet.py
        model_kwargs = {}  # initial conditions defaults given in thermal.py
    if update_kwargs is not None:
        pl_kwargs.update(update_kwargs)
    if run_kwargs is not None:
        model_kwargs.update(run_kwargs)

    # if logx:
    #     xvec = np.logspace(np.log10(xmin), np.log10(xmax), xres)
    # else:
    #     xvec = np.linspace(xmin, xmax, xres)
    # for x in xvec:
    #     pl_kwargs.update({xvar: x})

    pl_ensemble = evol.bulk_planets_mc(n=num, names=names, mini=mini, maxi=maxi, pl_kwargs=pl_kwargs,
                                       model_kwargs=model_kwargs, **kwargs)

    # plot ensemble
    if fig is None and axes is None:
        fig, axes = plt.subplots(len(yvars), 1, figsize=(2.5, 2 * len(yvars)))

    for ii, yvar in enumerate(yvars):
        try:
            ax = axes[ii]
        except TypeError:
            ax = axes
        y_all = []
        for pl in pl_ensemble:
            t = pl.t * parameters.sec2Gyr
            y = eval('pl.' + yvar) * yscales[ii]
            y_all.append(y)
            ax.plot(t, y, c=c, lw=lw, alpha=alpha)

        # mean, std
        y_all = np.array(y_all)
        y_av = np.mean(y_all, axis=0)
        ax.plot(t, y_av, c=c_mean, lw=4)
        y_std = np.std(y_all, axis=0)
        y_upper = y_av + y_std * n_sigma  # todo for log scape
        y_lower = y_av - y_std * n_sigma
        ax.plot(y, y_lower, c=c_mean, lw=1, ls='--')
        ax.plot(y, y_upper, c=c_mean, lw=1, ls='--')

        # format
        if ii == len(yvars) - 1:
            ax.set_xlabel(xlabel, fontsize=labelsize, labelpad=xlabelpad)
            ax.tick_params(axis='x', labelsize=ticksize, pad=tickpad)
            if xticks is not None:
                ax.set_xticks(xticks)
        else:
            ax.set_xticks([])

        ax.set_ylabel(ylabels[ii], fontsize=labelsize, labelpad=ylabelpad)
        if log[ii]:
            ax.set_yscale('log')
            ax.tick_params(axis='both', which='minor', labelsize=0)
        if yticks is not None:
            ax.set_yticks(yticks[ii])
        ax.tick_params(axis='y', labelsize=ticksize, pad=tickpad)
        ax.set_xlim(t[0], t[-1])
        if ylims[ii] is not None:
            ax.set_ylim(ylims[ii])

        if legtext is not None and ii == 0:
            ax.text(0.95, 0.95, legtext, fontsize=legsize, transform=ax.transAxes, ha='right', va='top')

    plt.tight_layout()
    if save:
        plot_save(fig, fname, fig_path=fig_path)
    # else:
    #     plt.show()
    return fig, axes


def plot_h_v_obvs(default='baseline', age=4.5, labelsize=28, legsize=16, ticksize=20, xlabelpad=20, fig_path='',
                  show_ss=False, leg=True,
                  ylabel='$\Delta h_{rms}$ (m)', log=True, nplanets=20, save=False,
                  x_vars=['t', 'M_p', 'H_0', 'CMF', 'Ea'],
                  units=['Gyr', '$M_E$', 'pW kg$^{-1}$', 'CMF', 'kJ mol$^{-1}$'],
                  x_range=[(1.5, 4.5), (0.1 * parameters.M_E, 6 * parameters.M_E), (10e-12, 40e-12),
                           (0.1, 0.7), (250e3, 350e3)],
                  xticks=[[1, 2, 4], [0.1, 1, 6], [10, 20, 40], [0.1, 0.2, 0.5], [300, 350]],
                  xscales=[parameters.sec2Gyr, parameters.M_E ** -1, 1e12, 1, 1e-3], ylim=(500, 1000),
                  yticks=[500, 600, 700, 800, 900, 1000],
                  xlabels=['Age\n(Gyr)', 'Planet mass\n($M_E$)',
                           'Rad. heating $t_0$\n(pW kg$^{-1}$)',
                           'Core mass fraction', 'Activation energy\n(kJ mol$^{-1}$)'], dark=False, show_Huang=False,
                  models=['dyn_top_rms'],
                  **kwargs):
    # how h varies across key input parameters
    textc = 'k'
    if dark:
        textc = 'w'

    fig, axes = plot_h_relative_multi(default=default, age=age, alpha=1, wspace=0.15, legsize=legsize,
                                      ticksize=ticksize, labelsize=labelsize, fig_height=6,
                                      yset_ylim=False, legend=True, labels=x_vars,
                                      lw=4,  # legtitle=r'\textbf{\textit{Model}}',
                                      ylabel=ylabel,
                                      nplanets=nplanets, log=False,
                                      relative=False, x_vars=x_vars, units=units,
                                      labelpad=20, legendtop=True, tickwidth=2,
                                      initial_kwargs={'T_m0': 1000, 'T_c0': 3000},
                                      models=models,
                                      x_range=x_range, xscales=xscales, xlabels=xlabels,
                                      c=['#d88868', '#749af3'], textc=textc, **kwargs)

    for ax in axes:
        if log:
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
            ax.xaxis.set_minor_formatter(ticker.NullFormatter())
        ax.set_ylim(ylim)
    axes[0].yaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
    axes[0].yaxis.set_minor_formatter(ticker.NullFormatter())
    # axes[0].set_ylabel('$\Delta h_{rms}$ (m)', fontsize=labelsize, c='xkcd:off white')
    axes[0].set_yticks(yticks)

    if leg:
        handles = [mlines.Line2D([], [], color='#d88868', ls='-',
                                 markersize=0, lw=4, label='$\Delta h^\prime = 0.094$ Ra$_{i, eff}^{-0.151}$')]

        if show_ss:
            handles.append(
                mlines.Line2D([], [], color='xkcd:goldenrod', marker='$V$',
                              markersize=15, lw=0, label=r'Venus'))
        if show_Huang:
            handles.append(mlines.Line2D([], [], color='xkcd:orchid', marker='^',
                                         markersize=15, lw=0, label=r'Huang+ (2013) 3D model'))

        legend = axes[0].legend(handles=handles, frameon=False, fontsize=legsize,
                                borderaxespad=0,
                                loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=3, )

    if dark:
        fig, *axes = dark_background(fig, axes)

    for i, ax in enumerate(axes):
        if x_vars[i] == 't':
            ax.set_xlim(x_range[i])
        else:
            ax.set_xlim((x_range[i][0] * xscales[i], x_range[i][1] * xscales[i]))
        if xticks is not None:
            ax.set_xticks(xticks[i])

    if show_ss:
        # VENUS: 850 m
        h_Venus = 865.4906656355711
        M_Venus = 0.815

        h_Mars = 6688.627942023225
        M_Mars = 0.107

        ax = axes[1]
        # imscatter(M_Venus, h_Venus, planet_icon_path + 'Venus.png', zoom=0.04, ax=ax)
        # imscatter(M_Mars, h_Mars, planet_icon_path + 'Mars.png', zoom=0.08, ax=ax)

        ax.scatter(M_Venus, h_Venus, marker='$V$', c='xkcd:goldenrod', s=200, zorder=100)

    if show_Huang:
        # Huang cases 1-13, 15
        h_Huang = np.array(
            [200.15279436132423, 688.2014927583677, 673.7880493468331, 402.07565967751117, 695.2136989391211,
             672.4561163950626, 214.12066607342535, 488.4601789919337, 878.5607285545191, 292.43829959982384,
             311.3352436867767, 339.3664129742059, 640.1361418805931, 430.1894190342128])
        for h in h_Huang:
            ax.scatter(M_Venus, h, marker='^', s=70, alpha=0.5, c='xkcd:orchid', label=r'Huang+ (2013)', zorder=1)

    if save:
        plot_save(fig, fname='h_obvs', fig_path=fig_path)
    else:
        plt.show()
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
            # print('1/yscale', 1/yscale)
            try:
                for ip, pl in enumerate(lplanets):  # planets to plot as line
                    data_x = eval('pl.' + xparam) * xscale
                    if isinstance(data_x, Iterable):
                        t = pl.t
                        it = age_index(t, snap, parameters.sec2Gyr)
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
                        t = pl.t
                        it = age_index(t, snap, parameters.sec2Gyr)
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
        if yparam[ii] in ['eta_m', 'nu_m', 'Ra_i_eff']:
            ax.set_yscale('log')
        if xparam in ['eta_m', 'nu_m', 'eta_0', 'Ra_i_eff']:
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


def plot_change_with_observeables_ensemble(defaults='Earthbaseline', wspace=0.1, tickwidth=1, textc='k',
                                           age=4.5, x_vars=None, ylabel='$\Delta h$ / $\Delta h_0$  ', fig_height=4,
                                           dist_res=10, ylim=None, leg_loc='upper left',
                                           xlabels=None, log=None, x_range=None, xscales=None, units=None, x_res=8,
                                           fig=None, axes=None, model_param='dyn_top_rms', legend=False, legsize=12,
                                           yscale=1, alpha=0.2,
                                           linec='k', labelsize=16, lw=3, ticksize=12,
                                           update_kwargs={}, initial_kwargs={}, verbose=False, **kwargs):
    if x_vars is None:
        x_vars = ['t', 'M_p', 'CMF', 'H_0']
    if units is None:
        units = ['Gyr', '$M_E$', 'CMF', 'pW kg$^{-1}$']
    if xlabels is None:
        xlabels = x_vars
    if x_range is None:
        x_range = [(1.5, age), (0.1 * parameters.M_E, 6 * parameters.M_E), (0.1, 0.7), (10e-12, 40e-12)]
    if xscales is None:
        xscales = [parameters.sec2Gyr, parameters.M_E ** -1, 1, 1e12]
    if axes is None:
        fig, axes = plt.subplots(1, len(x_vars), figsize=(6 * len(x_vars), fig_height), sharey=True)
    if not_iterable(axes):
        axes = [axes]
    if log is None:
        log = [False] * len(x_vars)

    # pl_baseline = evol.build_planet_from_id(ident=defaults,
    #                                         initial_kwargs=initial_kwargs, update_kwargs=update_kwargs,
    #                                         postprocessors=['topography'], t_eval=None)

    for i_ax, x_var in enumerate(x_vars):
        print('axis', i_ax + 1, '/', len(axes))
        xmin, xmax = x_range[i_ax]
        if x_var == 't':
            # time/age variation - plot single planet evol
            x_vec, y_av, y_upper, y_lower = ensemble_time_distribution(yvar=model_param, xvar=x_var,
                                                                       default=defaults, dist_res=dist_res,
                                                                       update_kwargs=update_kwargs,
                                                                       run_kwargs=initial_kwargs,
                                                                       yscale=yscale,
                                                                       min_t=xmin, max_t=xmax,
                                                                       **kwargs)
        else:
            if verbose:
                print('generating planets across', x_var, '...')
            x_vec, y_av, y_upper, y_lower = ensemble_marginal_distribution(yvar=model_param, xvar=x_var,
                                                                           default=defaults, dist_res=dist_res,
                                                                           update_kwargs=update_kwargs,
                                                                           run_kwargs=initial_kwargs,
                                                                           yscale=yscale, age=age,
                                                                           x_res=x_res, minx=xmin, maxx=xmax,
                                                                           log=log[i_ax],
                                                                           **kwargs)

        print('      range:', y_av[0], '-', y_av[-1], '| % diff:',  abs(y_av[-1] - y_av[0])/y_av[0])
        x_vec = x_vec * xscales[i_ax]
        axes[i_ax].plot(x_vec, y_av, c=linec, lw=lw)
        # print(x_var, ': x', x_vec, 'y', y_av)
        axes[i_ax].fill_between(x_vec, y_lower, y_upper, color=linec, alpha=alpha)
        axes[i_ax].set_xlabel(xlabels[i_ax], fontsize=labelsize)
        axes[i_ax].tick_params(axis='both', labelsize=ticksize)
        if log[i_ax]:
            axes[i_ax].set_xscale('log')
        if ylim is not None:
            axes[i_ax].set_ylim(ylim)
        if legend:
            string = ''
            if 't' not in x_vars:
                string = string + '{:.3g} '.format(age) + 'Gyr' + '\n'
            for jj, u in enumerate(units):
                if jj != i_ax:
                    if x_vars[jj] == 't':
                        string = string + '{:.3g} '.format(age) + u + '\n'
                    else:
                        string = string + '{:.3g} '.format(
                            eval('inputs.' + defaults + '_in')[x_vars[jj]] * xscales[jj]) + u + '\n'
            string = string[:-1]  # remove last \n
            if leg_loc == 'upper left':
                axes[i_ax].text(0.06, 0.96, string,
                                fontsize=legsize, c=textc,
                                horizontalalignment='left',
                                verticalalignment='top',
                                transform=axes[i_ax].transAxes)
            elif leg_loc == 'upper right':
                axes[i_ax].text(0.96, 0.96, string,
                                fontsize=legsize, c=textc,
                                horizontalalignment='right',
                                verticalalignment='top',
                                transform=axes[i_ax].transAxes)
    axes[0].set_ylabel(ylabel, fontsize=labelsize)

    for ax in axes:
        ax.xaxis.set_tick_params(width=tickwidth)
        ax.yaxis.set_tick_params(width=tickwidth)
    # plt.subplots_adjust(wspace=wspace)
    return fig, axes


def plot_change_with_observeables(defaults='Earthbaseline', wspace=0.1, tickwidth=1, relative=True, textc='k', c='k',
                                  age=4.5, x_vars=None, ylabel='$\Delta h$ / $\Delta h_0$  ', fig_height=4,
                                  xlabels=None, nplanets=20, log=False, x_range=None, xscales=None, units=None,
                                  fig=None, axes=None, model_param='dyn_top_rms', legend=False, legsize=12, yscale=1,
                                  pl_baseline=None, update_kwargs={}, initial_kwargs={}, relval=None, verbose=False, **kwargs):
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

    if type(model_param) == str:
        model_param = [model_param]
        c = [c]

    if relative and (relval is None):
        try:
            it = age_index(pl_baseline.t, age, parameters.sec2Gyr)
            model_baseline = eval('pl_baseline.' + model_param[0])[it]
        except IndexError:
            # scalar
            model_baseline = eval('pl_baseline.' + model_param[0])
        except AttributeError:
            # steady state?
            model_baseline = eval('pl_baseline.' + model_param[0])
        yscale = model_baseline ** -1
    elif relative:
        yscale = relval ** -1

    set_ylabel = True
    legendd = legend
    for i_ax, x_var in enumerate(x_vars):
        xmin, xmax = x_range[i_ax]
        if x_var == 't':
            # time/age variation - plot single planet evol
            for ip, param in enumerate(model_param):
                print('yscale', yscale)
                fig, ax = plot_vs_x(lplanets=pl_baseline, xname={'t': (xlabels[i_ax], xscales[i_ax])},
                                    ynames={param: (ylabel, yscale)}, fig=fig, axes=axes[i_ax], legsize=legsize,
                                    legend=legendd, plots_save=False, set_ylabel=set_ylabel, set_xlim=True, xmin=xmin,
                                    xmax=xmax, log=log, relative=relative, c=c[ip], **kwargs)
        else:
            if verbose:
                print('generating planets across', x_var, '...')
            if hasattr(pl_baseline, 't'):
                t_eval = pl_baseline.t
            else:
                t_eval = None
            planets = evol.bulk_planets(n=nplanets, name=x_var, mini=xmin, maxi=xmax,
                                        like=defaults,
                                        t_eval=t_eval, random=False, verbose=verbose,
                                        initial_kwargs=initial_kwargs, update_kwargs=update_kwargs, **kwargs)

            # for pl in planets:
            #     print('mass:', pl.M_p/parameters.M_E, '| rel. vol:', pl.max_ocean*yscale, '| mass ocn:',
            #           pl.max_ocean*1000/parameters.M_E, 'M_E | mass frac:',
            #           pl.max_ocean*1000/pl.M_p)

            for ip, param in enumerate(model_param):
                print(param, 'pl[-1]', eval('planets[-1].' + param))
                fig, ax = plot_vs_x(lplanets=planets, xname={x_var: (xlabels[i_ax], xscales[i_ax])},
                                    ynames={param: ('', yscale)}, fig=fig, axes=axes[i_ax], legsize=legsize,
                                    legend=legendd, snap=age, plots_save=False, set_ylabel=ylabel, set_xlim=True, log=log,
                                    relative=relative, c=c[ip], **kwargs)

        # if relative:
        #     ax.axhline(y=1, lw=1, alpha=0.5, zorder=0)
        if legend:
            string = ''
            for jj, u in enumerate(units):
                if jj != i_ax:
                    if x_vars[jj] == 't':
                        string = string + '{:.3g} '.format(age) + u + '\n'
                    else:
                        string = string + '{:.3g} '.format(eval('pl_baseline.' + x_vars[jj]) * xscales[jj]) + u + '\n'
            string = string[:-1]  # remove last \n
            ax.text(0.02, 0.96, string,
                    fontsize=legsize, c=textc,
                    horizontalalignment='left',
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


def plot_h_ensemble(x_vars=None, x_range=None, defaults='Earthbaseline', save=False, fname='relative_h',
                    models=None, labels=None, c=None, fig=None, axes=None, age=4.5,
                    initial_kwargs={}, update_kwargs={}, legend=True,
                    ylabel='$\Delta h$ / $\Delta h_0$  ', **kwargs):
    if x_vars is None:
        x_vars = ['t', 'M_p', 'H_0', 'CMF']
    if x_range is None:
        x_range = [(1.5, 4.5), (0.1 * parameters.M_E, 6 * parameters.M_E), (10e-12, 40e-12), (0.1, 0.7)]

    initial_kwargs.update({'tf': age})
    pl_baseline = evol.build_planet_from_id(ident=defaults,
                                            initial_kwargs=initial_kwargs, update_kwargs=update_kwargs,
                                            postprocessors=['topography'], t_eval=None)  # not used
    legendd = False
    for ii, h_param in enumerate(models):
        if ii == len(models) - 1:
            legendd = True

        fig, axes = plot_change_with_observeables(defaults=defaults, relative=False, fig=fig, axes=axes,
                                                  model_param=h_param,
                                                  legend=legendd, pl_baseline=pl_baseline, label_l=labels[ii], c=c[ii],
                                                  age=age, ylabel=ylabel, initial_kwargs=initial_kwargs,
                                                  update_kwargs=update_kwargs, **kwargs)

    if save:
        plot_save(fig, fname, **kwargs)

    return fig, axes


def plot_h_relative_multi(defaults='Earthbaseline', save=False, fname='relative_h',
                          models=None, labels=None, c=None, fig=None, axes=None, age=4.5,
                          initial_kwargs={}, update_kwargs={}, legend=True,
                          ylabel='$\Delta h$ / $\Delta h_0$  ', **kwargs):
    initial_kwargs.update({'tf': age})
    pl_baseline = evol.build_planet_from_id(ident=defaults,
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
                                 relative=True, fig=None, axes=None, vol_0=None, simple_scaling=False, c2='g',
                                 mass_iax=0, leg_bbox=(1.7, 1.01), log=False, figsize=(10, 10), ytitle=1.1, x_range=None,
                                 cmap='terrain_r', vmin=None, vmax=None, mass_EO=1.4e21, rho_w=1000, version=0, alpha_w=0.5,
                                 defaults='Venusbaseline', ylabel=r'$V_{\mathrm{max}}/V_{\mathrm{max, Ve}}$', **kwargs):
    # phi0, degree = harm.load_spectrum(fpath=spectrum_fpath, fname=spectrum_fname)
    # h_rms0 = harm.powerspectrum_RMS(power_lm=phi0, degree=degree)
    degree, phi0 = sh.load_model_spectrum_pkl(fname=spectrum_fname, path=spectrum_fpath)

    print('\n...........\nfname', spectrum_fname, 'phi0', phi0[:5])
    pl0 = evol.bulk_planets(n=1, name='M_p', mini=M0 * parameters.M_E, maxi=M0 * parameters.M_E, like=defaults,
                            # verbose=True,
                            t_eval=None, random=False, phi0=phi0, postprocessors=['topography'],
                            **kwargs)[0]
    pl0 = oceans.max_ocean(pl0, at_age=age, #name_rms='dyn_top_rms',
                           phi0=phi0, **kwargs)
    pl0 = oceans.simple_vol_scaling(pl0, at_age=age, #name_rms='dyn_top_rms',
                              **kwargs)
    # print('max ocean', pl0.max_ocean, 'm3')
    # print('peak scaling ocean', pl0.simple_ocean, 'm3')

    if vol_0 is None:
        vol_0 = pl0.max_ocean
    elif vol_0 == 'Earth':
        vol_0 = mass_EO / rho_w
    print('vol_0', vol_0)

    if axes is None:
        fig, axes = plt.subplots(figsize=figsize)
    if simple_scaling:
        model_param = ['max_ocean', 'simple_ocean']
        c = [c, c2]
    else:
        model_param = 'max_ocean'

    fig, axes = plot_change_with_observeables(defaults=defaults, model_param=model_param, legend=True, pl_baseline=pl0,
                                              textc=textc, at_age=age, x_range=x_range,
                                              label_l=None, c=c, ylabel=ylabel, age=age, legsize=legsize,
                                              postprocessors=['topography', 'ocean_capacity'], phi0=phi0, log=log,
                                              fig=fig, relval=vol_0,
                                              axes=axes, ticksize=ticksize, labelsize=labelsize, relative=relative,
                                               **kwargs)


    if mass_frac_sfcwater is not None:
        # how does actual vol scale assuming constant mass fraction of surface water (bad assumption)?
        ax = axes[mass_iax]
        masses = np.logspace(np.log10(0.1), np.log10(6), num=60)  # mass in M_E

        if version == 0:
            f_water = np.logspace(np.log10(mass_frac_sfcwater[0]), np.log10(mass_frac_sfcwater[-1]), num=25)
            colours = colorize([np.log10(m) for m in f_water], cmap=cmap)[0]
            for ii, X in enumerate(f_water):
                M_w = masses * parameters.M_E * X  # mass of sfc water in kg
                vol_w = M_w / rho_w  # corresponding volume
                if relative:
                    vol_w = vol_w / vol_0
                # print('relative water budget change', vol_w)
                ax.plot(masses, vol_w, alpha=alpha_w, lw=0, zorder=0, c=colours[ii], marker='o', markersize=15)
        elif version == 1:
            # contourfill with some individual lines
            f_water = np.logspace(np.log10(mass_frac_sfcwater[0]), np.log10(mass_frac_sfcwater[-1]), num=100)
            vol_ws = np.zeros((len(masses), len(f_water)))
            colours = colorize([np.log10(m) for m in f_water], cmap=cmap)[0]
            for ii, X in enumerate(f_water):
                M_w = masses * parameters.M_E * X  # mass of sfc water in kg
                vol_w = M_w / rho_w  # corresponding volume
                if relative:
                    vol_w = vol_w / vol_0
                vol_ws[:, ii] = vol_w
                if np.mod(ii, 14) == 0:
                    ax.plot(masses, vol_w, alpha=1, lw=0.5, zorder=1, c=colours[ii])
                ax.plot(masses[::14], vol_w[::14], alpha=alpha_w, lw=0, markersize=10, zorder=0, c=colours[ii])
            # ax.pcolormesh(masses, vol_ws, vol_ws, alpha=0.2, zorder=0, cmap=cmap)

        if vmin is None:
            vmin = np.min(mass_frac_sfcwater)
        if vmax is None:
            vmax = np.max(mass_frac_sfcwater)
        colourbar(mappable=None, vector=mass_frac_sfcwater, ax=ax, vmin=vmin,
                  vmax=vmax, label=clabel,
                  labelsize=labelsize, ticksize=ticksize, labelpad=clabelpad,
                  #ticks=mass_frac_sfcwater,
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


def read_JFR(fname='', path='/home/claire/Works/exo-top/benchmarks/JFR/'):
    df = pd.read_csv(path + fname, header=0, index_col=False)
    Ra = np.array(df.Ra)
    h_peak = np.array(df.peak_topo)
    h_rms = np.array(df.RMS_topo)
    Nu = np.array(df.Nu)
    return Ra, h_peak, h_rms, Nu
