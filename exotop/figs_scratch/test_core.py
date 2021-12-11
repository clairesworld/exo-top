import numpy as np
from model_1D.core import TerrestrialCore, show_inner_core, plot_flux_evol, core_params_default, print_core
from model_1D import evolve as evol
from model_1D.thermal import potential_temperature
import matplotlib.pyplot as plt
from model_1D import parameters as p
from model_1D import rheology as rh
from matplotlib import rc
import pandas as pd
from scipy.interpolate import interp1d

rc('text', usetex=False)


def bench_static(**kwargs):
    c = TerrestrialCore(test='Nimmo_static')  # , **core_params_default, **test_params_default)
    c.core_gradient(t=4.5 / p.sec2Gyr, **kwargs)  # Q_r=2.9e12,
    print_core(c)


def plot_bench(ax, fname, compare_dir='', c='r', ls='-', lw=1, alpha=0.6, label=None):
    df = pd.read_csv(compare_dir + fname + '.csv', header=None, names=['x', 'y'],
                     index_col=False)
    ax.plot(df['x'], df['y'], c=c, ls=ls, lw=lw, label=label, alpha=alpha)
    # if fname == 'T_mp':
    #     eta_m = rh.eta_Nim(df['y'].to_numpy(), eta_0=1e21, zeta=1e-2, T_ref=1573)
    #     print('eta m using Nimmo Tm prime =', eta_m[-1])
    return ax


def bench_evolution(compare_dir='/home/claire/Works/exo-top/benchmarks/nimmo_Earth/',
                    planet_kwargs={'CMF': 0.32264755186184685, 'R_c0': 3480e3, 'R_p0': 6400e3, 'alpha_m': 2.2e-5,
                                   'rho_m': 4800, 'c_m': 1200, 'T_s': 293, 'g_sfc': 9.8,
                                   'kappa_m': 6e-7, 'Ra_crit_u': 600, 'Ra_crit_c0': 600,
                                   'kappa_lm': 10e-7, 'k_m': 3.456, 'k_lm': 5.76,
                                   'c_i': np.array([30.4e-9, 22.7e-9, 0.26e-9, 120e-9])},  # tweak so you match
                    core_kwargs={'test': 'Nimmo_evol'}, visc_type='Nim', **kwargs):
    pl = evol.build_planet_from_id(ident='baseline', use_core=True, t0_buffer=0,
                                   # t_eval=np.arange(0, 4.5, 0.004)/p.sec2Gyr,
                                   run_kwargs={'use_lid': False, 'visc_type': visc_type},
                                   update_kwargs=planet_kwargs, core_params=core_kwargs,
                                   )
    c = pl.core
    print_core(c)
    print('pl.eta_m', pl.eta_m[-1], 'Pa s')
    print('pl.eta_cmb', pl.eta_cmb[-1], 'Pa s')

    R_i = np.squeeze(pl.core.R_ic)
    print('dR_i/dt', np.gradient(R_i*1e-3, pl.t*p.sec2Gyr), 'km/Gyr')

    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot(pl.t*p.sec2Gyr, R_i*1e-3 - np.mean(R_i*1e-3), label='R_i')
    ax1.plot(pl.t*p.sec2Gyr, pl.core.T_cmb - np.mean(pl.core.T_cmb), label='T_c')
    ax1.legend()
    ax2.plot(pl.t * p.sec2Gyr, np.gradient(R_i*1e-3, pl.t*p.sec2Gyr), label='dR_i/dt')
    ax2.plot(pl.t * p.sec2Gyr,  np.gradient(np.squeeze(pl.core.T_cmb), pl.t*p.sec2Gyr), ls='--', label='dT_c/dt')
    ratio = np.gradient(R_i*1e-3, pl.t*p.sec2Gyr) / np.gradient(np.squeeze(pl.core.T_cmb), pl.t * p.sec2Gyr)
    ax2.plot(pl.t * p.sec2Gyr, ratio, label='C_r')
    ax2.set_ylim((-20, 20))
    ax2.legend()
    # print('C_r inferred', ratio, 'km/K')


    # plot_flux_evol(c=c, t=pl.t, pl=None, **kwargs)
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(pl.t * p.sec2Gyr, pl.core.T_ic, 'k:', lw=1, label='$T_{icb}$')
    ax1.plot(pl.t * p.sec2Gyr, pl.core.T_cmb, 'k-', lw=2, label='$T_c$')
    ax1.plot(pl.t * p.sec2Gyr, pl.T_m, 'k--', lw=2, label='$T_m$')
    ax1.plot(pl.t * p.sec2Gyr,
             potential_temperature(T=pl.core.T_cmb, z=pl.R_p - pl.R_c, g=pl.g_sfc, Cp=pl.c_m, alpha=pl.alpha_m),
             'k-', lw=0.5, label='$T_c^\prime$')
    ax1.plot(pl.t * p.sec2Gyr, pl.Tp, 'k--', lw=0.5, label='$T_m^\prime$')
    ax1.legend()

    ax2.plot(pl.t * p.sec2Gyr, pl.Q_ubl * 1e-12, 'k--', lw=2, label='$Q_M$')
    ax2.plot(pl.t * p.sec2Gyr, pl.H_rad_m * 1e-12, 'k:', lw=1, label='$H_M M_M$')
    ax2.plot(pl.t * p.sec2Gyr, pl.Q_core * 1e-12, 'k-', lw=2, label='$Q_C$')
    ax2.legend()

    # inner core age
    tau_ic = c.inner_core_age(pl.t)
    ax1.axvline(tau_ic * p.sec2Gyr, c='k', lw=0.5, label='inner core age')

    ax1.set_title('core and mantle evolution')
    ax2.set_xlabel('time (Gyr)')
    ax1.set_ylabel('T (K)')
    ax2.set_ylabel('Power (TW)')
    ax1.set_ylim((1000, 8000))
    ax2.set_ylim((-20, 200))

    # show_inner_core(c=c, t_plot=tau_ic, times=pl.t)
    # show_inner_core(c=c, t_plot=tau_ic - 0.1 / p.sec2Gyr, times=pl.t)

    # benchmarking data
    if compare_dir is not None:
        ax1.axvline(3.5, c='r', lw=0.5)  # Nimmo
        ls = [':', '-', '--', '-', '--']
        lw = [1, 2, 2, 0.5, 0.5]
        for ii, par in enumerate(['T_ic', 'T_c', 'T_m', 'T_cp', 'T_mp']):
            ax1 = plot_bench(ax1, par, compare_dir=compare_dir, c='r', ls=ls[ii], lw=lw[ii], label=None)

        ls = ['--', ':', '-']
        lw = [2, 1, 2]
        for ii, par in enumerate(['Q_ubl', 'H_rad_m', 'Q_core']):
            ax2 = plot_bench(ax2, par, compare_dir=compare_dir, c='r', ls=ls[ii], lw=lw[ii], label=None)

        # test viscosity
        T_mp = pd.read_csv(compare_dir + 'T_mp' + '.csv', header=None, names=['x', 'y'],
                           index_col=False)['y'].to_numpy()
        eta_m = rh.eta_Nim(T=T_mp, T_ref=1573)
        # print('Nimmo present day T_mp is', T_mp[-1], 'K')
        print('eta m using Nimmo temperatures =', eta_m[-1], 'should be 6.7e20 Pa s |', rh.eta_Nim(T=1613, T_ref=1573))

        T_m = pd.read_csv(compare_dir + 'T_m' + '.csv', header=None, names=['x', 'y'],
                          index_col=False)['y'].to_numpy()[-1]
        T_c = pd.read_csv(compare_dir + 'T_c' + '.csv', header=None, names=['x', 'y'],
                          index_col=False)['y'].to_numpy()[-1]
        eta_cmb = 10 * rh.eta_Nim(T=(T_c + T_m) / 2, T_ref=3400)
        # show_inner_core(c=c, title=str(3.5) + ' Gyr', t_plot=3.5 / p.sec2Gyr, times=pl.t)

        # print('eta cmb using Nimmo temperatures =', eta_cmb, 'should be 6.7e21 Pa s')
        # print('z from mantle temperatures', np.log(T_mp[-1]/T_m) * pl.c_m/(-pl.alpha_m * pl.g_sfc) * 1e-3, 'km')
        # print('mantle depth', pl.d*1e-3, 'km', 'half depth', pl.d/2*1e-3, 'km')

    plot_flux_evol(c=pl.core, t=pl.t, pl=pl)

    print('M_p', pl.M_p / p.M_E, 'M_E | M_m', pl.M_m, 'M_p - M_c', pl.M_p - pl.M_c)
    plt.show()


bench_evolution(compare_dir=None,
                planet_kwargs={'kappa_m': 6e-7, 'kappa_lm': 10e-7, 'k_m': 3.456, 'k_lm': 5.76, 'x_Eu': 1},
                core_kwargs={'CMF': 0.32264755186184685},
                visc_type='KW')

# # L_test = np.sqrt((3 * core_params_default['K_0'] * np.log(core_params_default['rho_cen'] / core_params_default['rho_0'] + 1)) / (
# #                     2 * np.pi * p.G * core_params_default['rho_0'] * core_params_default['rho_cen']))
# # M_c = mass_profile(3480e3, rmin=0, rho_cen=12500, L=L_test)
# # print('predicted cmf:', M_c/p.M_E)


""" test coupling """


def plot_T_core_planet(M_p=1 * p.M_E, x_Eu=1, Ea=300e3, CMF=0.33, t_IC=None, **kwargs):
    pl = evol.build_planet_from_id(baseline='baseline', use_core=True, t0_buffer=None,
                                   update_kwargs={
                                       'M_p': M_p, 'Ea': Ea, 'x_Eu': x_Eu, 'CMF': CMF,
                                   }, postprocessors=None)

    plt.plot(pl.t * p.sec2Gyr, pl.core.T_cmb, 'k-', label='$T_c$ in core')
    plt.plot(pl.t * p.sec2Gyr, pl.T_c, 'g--', label='$T_c$ in planet')
    plt.plot(pl.t * p.sec2Gyr, pl.T_m, 'r--', label='$T_m$')
    plt.legend()
    plt.title('core and mantle evolution')
    plt.xlabel('time (Gyr)')
    plt.ylabel('T (K)')

    if t_IC is not None:
        show_inner_core(c=pl.core, t_plot=t_Gyr / p.sec2Gyr, times=pl.t)
    plot_flux_evol(c=pl.core, t=pl.t, pl=pl)
    plt.show()


def plot_fluxes_var(x_vec=None, x_name='CMF', core_params=None, pl_id='baseline',
                flux_colours=['xkcd:pink', 'xkcd:chartreuse', 'xkcd:peach', 'xkcd:cerulean', 'xkcd:scarlet', 'k'],
                **kwargs):
    fig, ax = plt.subplots(1, 1)
    if core_params is None:
        core_params = core_params_default
    if x_vec is None:
        x_vec = np.linspace(0.1, 0.5, num=10)  # default for CMF test

    for ii, x in enumerate(x_vec):
        print('\n\n---------------------------->', ii, x_name, '=', x)
        core_params.update({x_name: x})

        if pl_id is not None:
            pl = evol.build_planet_from_id(pl_id, use_core=True, t0_buffer=None, core_params=core_params,
                                           update_kwargs={x_name: x}, postprocessors=None)

            c = pl.core
        else:
            print('non coupled version - fixed cooling rate')
            c = TerrestrialCore(**core_params)
            c.solve(plot=False, **kwargs)

        Q_s = c.Q_s[-1]
        Q_l = c.Q_l[-1]
        Q_g = c.Q_g[-1]
        Q_r = c.Q_r[-1]
        Q_cmb = c.Q_cmb[-1]
        total_flux = Q_s + Q_l + Q_g + Q_r + Q_cmb

        fluxes = np.array([Q_s, Q_l, Q_g, Q_r, Q_cmb, total_flux]) * 1e-12
        labels = ['specific heat', 'latent heat', 'grav. heat', 'rad. heat', 'CMB loss', 'total']
        for jj in range(len(fluxes)):
            if ii == 0:
                label = labels[jj]
            else:
                label = None
            ax.plot(x, fluxes[jj], marker='x', c=flux_colours[jj], label=label)
        # print('R_ic', c.R_ic)

    ax.set_xlabel(x_name)
    ax.set_ylabel('flux (TW)')
    ax.legend()
    plt.show()

# x_vec, x_name = np.linspace(0.1, 0.5, num=9), 'CMF'
# x_vec, x_name = np.linspace(0.3, 3, num=9), 'x_Eu'
# # x_vec, x_name = np.linspace(0.5, 2, num=9)*p.M_E, 'M_p'
# plot_fluxes(x_vec=x_vec, x_name=x_name, core_params=core_params_default, pl_id='baseline')
