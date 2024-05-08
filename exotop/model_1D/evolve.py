import numpy as np
from . import terrestrialplanet as tp
from . import inputs
from . import parameters as p
from . import thermal as th
from . import rheology as rh
from . import topography
from . import oceans
from . import core as tcore
from scipy import integrate
import random as rand
from scipy.stats import loguniform
from pprint import pprint
import sys


def bulk_planets(n=1, name=None, mini=None, maxi=None, like=None, random=False,
                 run_kwargs=None, update_kwargs=None, logscale=False, **kwargs):
    """varying single parameter 'name' between mini and maxi, use default values otherwise.
    update_kwargs can include any TerrestrialPlanet attribute
    initial_kwargs can include T_m0, T_c0, D_l0, t0, tf"""

    if like is not None:
        pl_kwargs = eval('inputs.' + like + '_in').copy()
        model_kwargs = eval('inputs.' + like + '_run').copy()
    else:
        pl_kwargs = {}  # use defaults given in terrestrialplanet.py
        model_kwargs = {}  # initial conditions defaults given in thermal.py
    if update_kwargs is not None:
        pl_kwargs.update(update_kwargs)
    if run_kwargs is not None:
        model_kwargs.update(run_kwargs)
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

        pl = build_planet(new_kwargs, model_kwargs, **kwargs)
        planets.append(pl)
        ii += 1
    return planets


def bulk_planets_mc(n=100, names=None, mini=None, maxi=None, pl_kwargs={}, model_kwargs={}, t_eval=None, log=False,
                    T_m0_options=None, propagate_fit_err=False, verbose=False, beta_h=p.beta_h, cov_beta_h=p.cov_beta_h,
                    **kwargs):
    """varying multiple parameters in 'names' between mini and maxi, use default values otherwise.
    update_kwargs can include any TerrestrialPlanet attribute
    initial_kwargs can include T_m0, T_c0, D_l0, t0, tf. names, mini, and maxi are in order and must have same lenghts"""
    # pl_kwargs and model_kwargs should already be taken from inputs file
    if T_m0_options is None:
        T_m0_options = np.array([1000, 1250, 1500, 1750, 2000]) + 273  # [ 1273, 1523, 1773, 2023, 2273]

    planets = []
    ii = 0
    while ii < n:
        new_kwargs_pl = pl_kwargs.copy()
        new_kwargs_model = model_kwargs.copy()

        for iii, name in enumerate(names):
            if name == 'T_m0':
                # for initial temps do discrete values after johnny seales
                val = np.random.choice(T_m0_options, size=1)[0]
            elif log:
                val = loguniform.rvs(mini[iii], maxi[iii], size=1)
            else:
                val = rand.uniform(mini[iii], maxi[iii])

            if name in ['T_m0', 'T_c0', 'D_l0']:
                new_kwargs_model.update({name: val})
            else:
                new_kwargs_pl.update({name: val})
            if verbose:
                print(ii, '/', n - 1, ': drew random', name, val)
        if propagate_fit_err:
            beta_h_new = np.random.multivariate_normal(beta_h, cov_beta_h)
            new_kwargs_model.update({'beta_h': beta_h_new})
            if verbose:
                print('       & drew random errors on h scaling', beta_h_new)

        if t_eval is None and ii > 0:
            t_eval = planets[0].t
        pl = None
        while pl is None:  # returns None if numerical issue e.g. D_l0 negative
            pl = build_planet(planet_kwargs=new_kwargs_pl, run_kwargs=new_kwargs_model, t_eval=t_eval, verbose=verbose,
                              **kwargs)
        planets.append(pl)
        ii += 1
    return planets


def build_planet_from_id(ident='Earthbaseline', run_kwargs=None, update_kwargs=None, **kwargs):
    planet_kwargs = eval('inputs.' + ident + '_in').copy()
    model_kwargs = eval('inputs.' + ident + '_run').copy()
    if run_kwargs is not None:
        model_kwargs.update(run_kwargs)
    if update_kwargs is not None:
        planet_kwargs.update(update_kwargs)
    pl = build_planet(planet_kwargs=planet_kwargs, run_kwargs=model_kwargs, **kwargs)
    return pl


def build_planet(planet_kwargs, run_kwargs, solve_ODE=True, steady=False, use_core=False, core_params={}, **kwargs):
    if run_kwargs is None:
        run_kwargs = {}
    if planet_kwargs is None:
        planet_kwargs = {}
    pl = tp.TerrestrialPlanet(**planet_kwargs)
    if use_core:
        core_params.update({'M_p': pl.M_p, 'R_p': pl.R_p})
        pl.core = tcore.TerrestrialCore(**core_params)
        run_kwargs['T_c0'] = pl.core.T_cmb  # update
        run_kwargs['T_m0'] = pl.core.T_cmb - 500  # start mantle at core temp as in Nimmo
        # print('using initial T_cmb', run_kwargs['T_c0'], 'K', 'T_m', run_kwargs['T_m0'], 'K')

        # set planet structure to accord with (better) core structure
        pl.reinitialise(R_c=pl.core.R_c, CMF=pl.core.CMF, R_p=pl.core.R_p)

        # check
        for k, v in pl.core.__dict__.items():
            if hasattr(pl, k):
                if eval('pl.' + k) != v:
                    print('!value mismatch between', k, ': planet has', eval('pl.' + k), 'core has',
                          eval('pl.core.' + k))

    if steady:
        age = run_kwargs['tf'] / p.sec2Gyr
        pl = steadystate(pl, age=age, **run_kwargs, **planet_kwargs, **kwargs)
        pl = postprocess_planet(pl, **kwargs)
    elif solve_ODE:
        pl = solve(pl, run_kwargs=run_kwargs, use_core=use_core, **kwargs)  # T_m, T_c, D_l
        pl = postprocess_planet(pl, **kwargs, **run_kwargs)
    return pl


def postprocess_planet(pl, postprocessors=None, nondimensional=False, **kwargs):
    if postprocessors is None:
        postprocessors = ['topography']

    if 'topography' in postprocessors:
        pl = topography.topography(pl, **kwargs)
    if 'ocean_capacity' in postprocessors:
        pl = oceans.max_ocean_fast(pl, **kwargs)
        # pl = oceans.max_ocean(pl, **kwargs)
        pl = oceans.simple_vol_scaling(pl, **kwargs)
    if nondimensional:
        pl.nondimensionalise()
    return pl


def solve(pl, run_kwargs={}, t0=0, t_eval=None, t0_buffer=5, verbose=False, use_core=False, use_lid=True, atol=(1e-6, 1e-6, 1e-6),
          **kwargs):
    """ call ODE solver """

    tf = run_kwargs['tf']
    T_m0 = run_kwargs['T_m0']
    T_c0 = run_kwargs['T_c0']
    D_l0 = run_kwargs['D_l0']
    pl.age = tf

    if t0_buffer is not None and t0_buffer != 0:
        t0i = t0 - t0_buffer
        if verbose:
            print('getting initial conditions with', t0_buffer, 'Gyr buffer')
        pre_f = integrate.solve_ivp(
            fun=lambda t, y: LHS(t, y, **dict(pl=pl, use_core=use_core, use_lid=use_lid, **run_kwargs, **kwargs)),
            t_span=(t0i * 1e9 * p.years2sec, t0 * 1e9 * p.years2sec), y0=[T_m0, T_c0, D_l0],
            t_eval=None,
            max_step=100e6 * p.years2sec,
            method='BDF',  # 'RK45',
            # atol=atol,
        )
        T_m0, T_c0, D_l0 = pre_f.y[:, -1]

    if T_m0 > 10000 or D_l0 < 0:
        return None

    if verbose:
        print('Solving 1D thermal history with T_m0 =', T_m0, 'K, T_c0 =', T_c0, 'K, D_l0 =', D_l0, 'm, tf =', tf,
              'Gyr')

    f = integrate.solve_ivp(
        fun=lambda t, y: LHS(t, y, **dict(pl=pl, use_core=use_core, verbose=verbose, **run_kwargs, **kwargs)),
        t_span=(t0 * 1e9 * p.years2sec, tf * 1e9 * p.years2sec), y0=[T_m0, T_c0, D_l0],
        t_eval=t_eval,
        max_step=100e6 * p.years2sec,
        method='RK45',
        # atol=atol
    )

    # initial postprocessing to get values of derived parameters across time
    pl.T_m = f.y[0]
    pl.T_c = f.y[1]
    pl.D_l = f.y[2]
    pl.t = f.t
    pl = recalculate(pl.t, pl, verbose=verbose, use_core=use_core, **run_kwargs, **kwargs)
    if use_core:
        pl.core.T_cmb = pl.T_c
        pl.core.Q_cmb = pl.Q_core
        pl.core.evolve_profiles(t=pl.t)
        pl.core.core_gradient(t=pl.t)

    return pl


def LHS(t, y, pl=None, use_core=False, use_lid=True, **kwargs):
    """ ODE equation to solve, LHS = 0 """

    # print('t', t*p.sec2Gyr, 'Gyr, T_m', y[0])

    pl.T_m = y[0]
    pl.T_c = y[1]
    if use_lid:
        pl.D_l = y[2]
    else:
        pl.D_l = 0
    # pl.D_l = 0
    # print('lid = 0')
    # print('  LHS in planet: t =', t*p.sec2Gyr, 'Gyr', 'T_cmb', pl.T_c, 'K', 'T_m', pl.T_m, 'K')
    if np.isnan(pl.T_c):
        print('pl.T_c', pl.T_c, 'K')
        # sys.exit()
    pl = recalculate(t, pl, use_core=use_core, use_lid=use_lid, **kwargs)

    if not use_core:
        # print('not using core dTdc')
        dTdt_c = th.dTdt(-pl.Q_core, pl.M_c, pl.c_c)
    else:
        pl.core.Q_cmb = pl.Q_core
        # dTdt_c, dR_idt = tcore.temperature_LHS(t, y=[pl.T_c, pl.core.R_ic], c=pl.core, **kwargs)
        dTdt_c = tcore.temperature_LHS(t, y=pl.T_c, c=pl.core, **kwargs)

    # print('    dTdt_c in planet', dTdt_c)
    dTdt_m = th.dTdt(-pl.Q_ubl + pl.H_rad_m + pl.Q_core - pl.Q_melt, pl.M_conv, pl.c_m)
    # print('    dTdt_m in planet', dTdt_m)
    if use_lid:
        dDdt = th.lid_growth(T_m=pl.T_m, q_ubl=pl.q_ubl, h0=pl.h_rad_m, R_p=pl.R_p, R_l=pl.R_l, T_l=pl.T_l,
                             rho_m=pl.rho_m,
                             T_s=pl.T_s,
                             c_m=pl.c_m, k_m=pl.k_m, **kwargs)
    else:
        dDdt = 0
    # dDdt = 0
    # print('t =', t*p.sec2Gyr, 'Gyr, T_m =', pl.T_m, 'dD/dt =', dDdt*1e-3/p.sec2Gyr, 'km/Gyr', 'h =', pl.h_rad_m*1e12, 'pW/kg')

    return [dTdt_m, dTdt_c, dDdt]


def recalculate(t, pl, verbose=False, rad_type=None, q_type='Ra_crit', lid_heating=True, use_core=False, use_lid=True, **kwargs):
    """ calcualte all the thermal convection parameters (single time step or postprocess all timesteps"""

    # if verbose:
    #     print('\nrecalculate kwargs:\n', kwargs)
    #     print('\nplanet:')
    #     for attr in dir(pl):
    #         print("pl.%s = %r" % (attr, getattr(pl, attr))

    # check if there are any weird issues
    try:
        lid_too_thick = any([list(pl.D_l)]) > pl.R_p - pl.R_c
    except TypeError:
        lid_too_thick = pl.D_l > pl.R_p - pl.R_c
    if lid_too_thick:
        if verbose:
            print('{:.2f} M_E: convection stopped at {:.2f} Gyr because lid grew to core'.format(pl.M_p / p.M_E,
                                                                                                 t * p.sec2Gyr))
        pl.D_l = pl.R_p - pl.R_c

    # geometry update
    pl.R_l = pl.R_p - pl.D_l
    pl.d_m = pl.R_l - pl.R_c
    V_lid = 4 / 3 * np.pi * (pl.R_p ** 3 - pl.R_l ** 3)
    pl.M_lid = V_lid * pl.rho_m  # not considering density change in lid due to melting/differentiation
    pl.M_conv = pl.M_m - pl.M_lid  # mass of convecting region

    Tp = th.potential_temperature(T=pl.T_m, z=pl.R_p - pl.R_c, g=pl.g_sfc, Cp=pl.c_m,
                                  alpha=pl.alpha_m)  # halfway depth or CMB?
    pl.Tp = Tp
    # print('think lid thickness is', pl.D_l*1e-3, 'km', 'mantle depth', pl.d_m*1e-3, 'km', 'depth available', (pl.R_p - pl.R_c)*1e-3, 'km')


    # radiogenic heating
    # pl.h_rad_m = th.h_rad(t, H_0=pl.H_0, H_f=pl.H_f, c_n=pl.c_n, p_n=pl.p_n, lambda_n=pl.lambda_n, verbose=verbose, **kwargs)  # W kg^-1
    if rad_type == 'H0':
        # pl.h_rad_m = th.h_rad_old_old(t, age=pl.age)
        pl.h_rad_m = th.h_rad_old(t, H_0=pl.H_0, H_f=pl.H_f, c_n=pl.c_n, p_n=p.p_n, lambda_n=p.lambda_n,
                                  verbose=verbose, **kwargs)
        if np.size(t) > 1:  # and verbose:
            print('h rad using H0', pl.h_rad_m[0] * 1e12, pl.h_rad_m[-1] * 1e12)
    else:
        pl.h_rad_m = th.h_rad(t, pl.c_i, pl.h_i, pl.tau_i, pl.age, pl.x_Eu)  # W kg^-1
        # if np.size(t) > 1:
        #     print('h rad using xEu', pl.h_rad_m * 1e12)
    pl.a0 = pl.h_rad_m * pl.rho_m  # volumetric heating in W m^-3
    pl.H_rad_m = th.H_rad(h=pl.h_rad_m, M=pl.M_conv, **kwargs)  # mantle radiogenic heating in W
    pl.H_rad_lid = th.H_rad(h=pl.h_rad_m, M=pl.M_lid, **kwargs)  # lid radiogenic heating in W

    # temperature update

    if use_core:
        pl.dT_visc = p.R_b * Tp ** 2 / pl.Ea
    else:
        pl.dT_visc = p.R_b * pl.T_m ** 2 / pl.Ea
    pl.dT_rh = pl.a_rh * pl.dT_visc
    if not use_lid:
        pl.T_l = pl.T_s
    elif use_core:
        pl.T_l = th.T_lid(T_m=Tp, a_rh=pl.a_rh, dT_visc=pl.dT_visc)
    else:
        pl.T_l = th.T_lid(T_m=pl.T_m, a_rh=pl.a_rh, dT_visc=pl.dT_visc)
    pl.dT_m = pl.T_c - pl.T_l
    pl.delta_T = pl.T_c - pl.T_s
    pl.b = rh.visc_factor(pl, **kwargs)

    try:
        bad_T = any([list(pl.dT_m)]) < 0
    except TypeError:
        bad_T = pl.dT_m < 0
    if bad_T and verbose:
        print('{:.2f} M_E: SERIOUS WARNING: -ve T_m at {:.2f} Gyr'.format(pl.M_p / p.M_E, t * p.sec2Gyr))

    # viscosity update
    # pl.eta_s = rh.dynamic_viscosity(T=pl.T_s, pl=pl, **kwargs)  # viscosity at surface temperature
    # pl.eta_l = rh.dynamic_viscosity(T=pl.T_l, pl=pl, **kwargs)  # viscosity at lid base temperature
    if use_core:
        pl.eta_m = rh.dynamic_viscosity(T=pl.Tp, T_ref=1573,
                                        **kwargs)  # viscosity at interior mantle temperature
        pl.eta_cmb = 10 * rh.dynamic_viscosity(T=(pl.T_c + pl.T_m) / 2, T_ref=3400,
                                               **kwargs)  # viscosity at cmb - T_ref may not be used
    else:
        pl.eta_m = rh.dynamic_viscosity(T=pl.T_m, pl=pl, T_ref=None,
                                        **kwargs)  # viscosity at interior mantle temperature
        pl.eta_cmb = rh.dynamic_viscosity(T=(pl.T_c + pl.T_m) / 2, pl=pl, T_ref=None, **kwargs)  # viscosity at cmb

    # pl.delta_eta = pl.eta_s / rh.dynamic_viscosity(T=(pl.T_c - pl.T_s), pl=pl, **kwargs)  # total viscosity contrast

    # equivalent to effective Ra_i in 2D models with d and dT decreased by lid thickness and temp drop
    pl.Ra_i_eff = th.Ra(eta=pl.eta_m, rho=pl.rho_m, g=pl.g_sfc, deltaT=abs(pl.dT_m), l=pl.d_m, kappa=pl.kappa_m,
                        alpha=pl.alpha_m)

    # for calculating the upper tbl thickness
    if use_core:
        deltaT = pl.Tp - pl.T_l
    else:
        deltaT = pl.T_m - pl.T_l
    pl.Ra_rh_u = th.Ra(eta=pl.eta_m, rho=pl.rho_m, g=pl.g_sfc, deltaT=deltaT, l=pl.d_m, kappa=pl.kappa_m,
                       alpha=pl.alpha_m)

    # the 'overall' interior Ra covering the whole domain from cmb to surface
    pl.Ra_i = th.Ra(eta=pl.eta_m, rho=pl.rho_m, g=pl.g_sfc, kappa=pl.kappa_m, alpha=pl.alpha_m,
                    deltaT=pl.delta_T,
                    l=pl.d)
    pl.theta_FK = rh.theta_FK(Ea=pl.Ea, T_i=pl.T_m, T_s=pl.T_s)

    pl.delta_rh = th.bdy_thickness(d_m=pl.d_m, Ra_crit=pl.Ra_crit_u, beta=pl.beta_u, Ra_rh=pl.Ra_rh_u, **kwargs)
    # if use_core:
    #     delta_rh = (600 * pl.kappa_m * pl.eta_m / (pl.rho_m * pl.g_sfc * pl.alpha_m * (pl.Tp - pl.T_s)))**(1/3)
        # print('calculated delta_rh',pl.delta_rh*1e-3, 'should be', delta_rh*1e-3)


    if q_type == 'theta':  # i.e. like Foley & Smye
        pl.q_ubl = th.q_theta_scaling(theta_FK=pl.theta_FK, Ra_i=pl.Ra_i, T_m=pl.T_m, T_s=pl.T_s, d=pl.d, k=pl.k_m,
                                      c1=0.5)
    elif use_core:
        pl.q_ubl = th.q_bl(deltaT=Tp - pl.T_l, k=pl.k_m, d_bl=pl.delta_rh, beta=pl.beta_u, **kwargs)
    else:
        pl.q_ubl = th.q_bl(deltaT=pl.T_m - pl.T_l, k=pl.k_m, d_bl=pl.delta_rh, beta=pl.beta_u, **kwargs)
    pl.Q_ubl = th.Q_bl(q=pl.q_ubl, R=pl.R_l)

    # core
    if not hasattr(pl, 'Ra_crit_c0'):
        pl.Ra_crit_c = 0.28 * pl.Ra_i ** 0.21  # Thiriet eqn 16, Deschamps & Sotin (2001)
    else:
        pl.Ra_crit_c = pl.Ra_crit_c0
    pl.Ra_rh_c = th.Ra(eta=pl.eta_cmb, rho=pl.rho_m, g=pl.g_cmb, deltaT=abs(pl.T_c - pl.T_m), l=pl.d_m,
                       kappa=pl.kappa_m, alpha=pl.alpha_m)
    pl.TBL_c = th.bdy_thickness(d_m=pl.d_m, Ra_crit=pl.Ra_crit_c, Ra_rh=pl.Ra_rh_c, beta=pl.beta_c, **kwargs)

    # if use_core:
    #     delta_c = (600 * pl.kappa_lm * pl.eta_cmb / (pl.rho_m * pl.g_sfc * pl.alpha_m * pl.T_c - pl.T_m))**(1/3)
    #     print('calculated delta_b', pl.TBL_c * 1e-3, 'should be', delta_c * 1e-3)

    pl.dT_cmb = pl.T_c - pl.T_m
    pl.q_core = th.q_bl(deltaT=pl.dT_cmb, k=pl.k_lm, d_bl=pl.TBL_c, beta=pl.beta_c, **kwargs)
    pl.Q_core = th.Q_bl(q=pl.q_core, R=pl.R_c)

    pl.Ra_F = th.Ra_F(eta=pl.eta_m, kappa=pl.kappa_m, H=pl.h_rad_m, alpha=pl.alpha_m, k=pl.k_m, rho=pl.rho_m,
                      g=pl.g_sfc, l=pl.R_l - pl.R_c, F_b=pl.q_core)

    if lid_heating:
        a0_lid = pl.a0
    else:
        a0_lid = 0
    pl.q_sfc = th.sph_flux(pl.R_p, a0=a0_lid, T_l=pl.T_l, R_l=pl.R_l, k_m=pl.k_m, T_s=pl.T_s, R_p=pl.R_p,
                           **kwargs)  # W m^-2
    # pl.q_sfc = th.sfc_flux_Nu(eta=pl.eta_m, kappa=pl.kappa_m, alpha=pl.alpha_m, k=pl.k_m, rho=pl.rho_m, g_sfc=pl.g_sfc,
    #                           d=pl.d, deltaT_visc=pl.dT_visc)

    pl.Q_sfc = th.Q_bl(q=pl.q_sfc, R=pl.R_p)
    pl.urey = (pl.H_rad_m + pl.H_rad_lid) / pl.Q_sfc

    pl.T_avg = th.T_mean_old(T_m=pl.T_m, T_l=pl.T_l, R_p=pl.R_p, R_l=pl.R_l, R_c=pl.R_c, T_s=pl.T_s, k_m=pl.k_m,
                             a0=a0_lid, **kwargs)
    # pl.T_avg = th.T_mean(R_c=pl.R_c, R_l=pl.R_l, R_p=pl.R_p, k_m=pl.k_m, T_m=pl.T_m, T_l=pl.T_l, T_s=pl.T_s, a0=a0_lid,
    #                      delta_rh=pl.delta_rh)

    # melting
    pl = th.melting(pl, **kwargs)

    return pl


def steadystate(pl, age=4.5 / p.sec2Gyr, T_m0=2000, tol=1e-1, verbose=False, **kwargs):
    # ADAPTED FROM KITE+ 2009 - only works with reference viscosity model
    T_s = pl.T_s
    kappa_m = pl.kappa_m
    alpha_m = pl.alpha_m
    rho_m = pl.rho_m
    Ra_crit = pl.Ra_crit_u
    g = pl.g_sfc
    k_m = pl.k_m
    beta = pl.beta_u
    a_rh = pl.a_rh
    Ea = pl.Ea
    eta_0 = pl.eta_0
    T_ref = pl.T_ref  # for viscosity law
    H0 = pl.H_0
    R_p = pl.R_p
    R_c = pl.R_c
    SA_p = pl.SA_p
    d_m0 = R_p - pl.R_c  # mantle depth ignoring lid
    M_m0 = pl.M_m  # mantle mass ignoring lid
    kwargs.update({'backwards_cooling': False, 'visc_type': 'Thi'})
    h_t = th.h_rad(t=age, H_0=pl.H_0, c_n=pl.c_n, p_n=pl.p_n, lambda_n=pl.lambda_n,
                   **kwargs)  # rad heating per unit mass
    # print('h(t)', h_t*1e12, 'pW/kg')
    T_m1 = T_m0
    # d_m = d_m0
    diff = 10

    while diff > tol:
        # print('\n')
        # loop starts
        T_m = T_m1
        d_lid = th.d_lid_ss(T_m, a_rh, k_m, Ea, H0, Ra_crit, eta_0, T_ref, kappa_m, alpha_m, g, rho_m, T_s)
        # print('d_lid', d_lid*1e-3, 'km')
        R_l = R_p - d_lid  # steady state
        d_m = R_l - R_c
        V_lid = 4 / 3 * np.pi * (R_p ** 3 - R_l ** 3)
        M_lid = V_lid * pl.rho_m  # not considering density change in lid due to melting/differentiation
        M_m = M_m0 - M_lid  # mass of convecting region
        # print('M_m/M_m0', M_m/M_m0)

        eta_m = rh.dynamic_viscosity(T_m, visc_type='Thi')
        # print('eta_m', eta_m)
        T_l = th.T_lid(T_m, a_rh, Ea)
        Nu = (rho_m * g * alpha_m * (T_m - T_l) * d_m ** 3 / (kappa_m * eta_m * Ra_crit)) ** beta
        T_m1 = M_m / SA_p * h_t * d_m / (Nu * k_m) + T_l
        diff = abs(T_m - T_m1)
        # print('T_m1', T_m1)
    if verbose:
        print('T_m', T_m1, 'K | d_lid', d_lid * 1e-3, 'km')

    # save
    pl.T_m = T_m1
    pl.T_c = pl.T_m  # no core thermally
    pl.D_l = d_lid
    pl.eta_m = eta_m
    pl.d_m = d_m
    pl.R_l = R_l
    pl.T_l = T_l
    pl.Nu = Nu
    pl.dT_rh = pl.a_rh * (p.R_b * pl.T_m ** 2 / pl.Ea)
    pl.dT_m = pl.T_c - pl.T_l
    pl.Ra_i_eff = th.Ra(eta=pl.eta_m, rho=pl.rho_m, g=pl.g_sfc, deltaT=pl.dT_m, l=pl.d_m, kappa=pl.kappa_m,
                        alpha=pl.alpha_m)
    pl.Ra_i = th.Ra(eta=pl.eta_m, rho=pl.rho_m, g=pl.g_sfc, kappa=pl.kappa_m, alpha=pl.alpha_m,
                    deltaT=abs(pl.dT_m) + (pl.T_l - pl.T_s),
                    l=pl.R_p - pl.R_c)
    pl.delta_rh = th.bdy_thickness(d_m=pl.d_m, Ra_crit=pl.Ra_crit_u, beta=pl.beta_u, Ra_rh=pl.Ra_i_eff, **kwargs)
    pl.q_ubl = th.q_bl(deltaT=pl.T_m - pl.T_l, k=pl.k_m, d_bl=pl.delta_rh, beta=pl.beta_u, **kwargs)
    pl.eta_s = rh.dynamic_viscosity(T=pl.T_s, pl=pl, **kwargs)  # viscosity at surface temperature
    pl.delta_eta = pl.eta_s / rh.dynamic_viscosity(T=(pl.T_c - pl.T_s), pl=pl, **kwargs)  # total viscosity contrast

    return pl

# def steady_state(t, pl, tol=1e-1, T_m0=1700, D_l0=100e3, **kwargs):
#     """ steady state temps for Urey ratio 1 planet """
#     difference = 1
#     q_rad = th.h_rad(t, H_0=pl.H_0, H_f=pl.H_f, c_n=pl.c_n, p_n=pl.p_n, lambda_n=pl.lambda_n, **kwargs)
#     a0 = q_rad * pl.rho_m  # volumetric heating in W m^-3
#     D_l = D_l0
#     T_m = T_m0
#     while difference > tol:
#         T_c = T_m
#         R_l = pl.R_p - D_l
#         d_m = pl.R_l - pl.R_c
#
#         # heating
#         V_lid = 4 / 3 * np.pi * (pl.R_p ** 3 - pl.R_l ** 3)
#         M_lid = V_lid * pl.rho_m  # not considering density change in lid due to melting/differentiation
#         M_conv = pl.M_m - pl.M_lid  # mass of convecting region
#         H_rad = th.H_rad(q_rad, M_conv)  # mantle rad heating in W (equal to surface flux here)
#         H_rad_lid = th.H_rad(q_rad, M_lid)  # lid radiogenic heating in W
#
#         # temperature update
#         T_l = th.T_lid(T_m=T_m, a_rh=pl.a_rh, Ea=pl.Ea)
#         dT_rh = pl.a_rh * (p.R_b * T_m ** 2 / pl.Ea)
#         dT_m = T_c - T_l
#
#         eta_m = rh.dynamic_viscosity(T=T_m, pl=pl, **kwargs)  # viscosity at interior mantle temperature
#         Ra_i_eff = th.Ra(eta=eta_m, rho=pl.rho_m, g=pl.g_sfc, deltaT=dT_m, l=d_m, kappa=pl.kappa_m,
#                          alpha=pl.alpha_m)
#
#         delta_rh = th.bdy_thickness(d_m=d_m, Ra_crit=pl.Ra_crit_u, beta=1/3, Ra_rh=Ra_i_eff, **kwargs)
#         q_ubl = th.q_bl(deltaT=dT_m, k=pl.k_m, d_bl=delta_rh, beta=1/3, **kwargs)
#         Q_ubl = th.Q_bl(q=q_ubl, R_b=R_l)
#         q_sfc = th.sph_flux(pl.R_p, a0=pl.a0, T_l=T_l, R_l=R_l, k_m=pl.k_m, T_s=pl.T_s, R_p=pl.R_p, **kwargs)
#
#         Q_sfc = th.Q_bl(q=q_sfc, R_b=pl.R_p)  # W
#
#     # TODO: sympy?
