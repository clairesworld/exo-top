import numpy as np
from . import terrestrialplanet as tp
from . import inputs
from . import parameters as p
from . import thermal as th
from . import rheology as rh
from . import topography
from . import oceans
from scipy import integrate
import random as rand
from scipy.stats import loguniform
from pprint import pprint


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
                    T_m0_options=[1273, 1523, 1773, 2023, 2273], **kwargs):
    """varying multiple parameters in 'names' between mini and maxi, use default values otherwise.
    update_kwargs can include any TerrestrialPlanet attribute
    initial_kwargs can include T_m0, T_c0, D_l0, t0, tf. names, mini, and maxi are in order and must have same lenghts"""
    # pl_kwargs and model_kwargs should already be taken from inputs file

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

        if t_eval is None and ii > 0:
            t_eval = planets[0].t
        pl = build_planet(new_kwargs_pl, new_kwargs_model, t_eval=t_eval, **kwargs)
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
    pl = build_planet(planet_kwargs, model_kwargs, **kwargs)
    return pl


def build_planet(planet_kwargs, run_kwargs, solve_ODE=True, steady=False, **kwargs):
    if run_kwargs is None:
        run_kwargs = {}
    if planet_kwargs is None:
        planet_kwargs = {}
    pl = tp.TerrestrialPlanet(**planet_kwargs)
    if steady:
        age = run_kwargs['tf']/p.sec2Gyr
        pl = steadystate(pl, age=age, **run_kwargs, **planet_kwargs, **kwargs)
    elif solve_ODE:
        pl = solve(pl, run_kwargs=run_kwargs, **kwargs)  # T_m, T_c, D_l
    pl = postprocess_planet(pl, **kwargs)

    return pl


def postprocess_planet(pl, postprocessors=None, nondimensional=True, **kwargs):
    if postprocessors is None:
        postprocessors = ['topography']

    if 'topography' in postprocessors:
        pl = topography.topography(pl, **kwargs)
    if 'ocean_capacity' in postprocessors:
        pl = oceans.max_ocean(pl, **kwargs)
        pl = oceans.simple_vol_scaling(pl, **kwargs)
    if nondimensional:
        pl.nondimensionalise()
    return pl


def solve(pl, run_kwargs={}, t0=0, t_eval=None, verbose=False, **kwargs):
    """ call ODE solver """

    tf = run_kwargs['tf']
    T_m0 = run_kwargs['T_m0']
    T_c0 = run_kwargs['T_c0']
    D_l0 = run_kwargs['D_l0']

    if verbose:
        print('Solving 1D thermal history with T_m0 =', T_m0, 'K, T_c0 =', T_c0, 'K, D_l0 =', D_l0, 'm, tf =', tf,
              'Gyr')
        back = run_kwargs['backwards_cooling']
        if back:
            print('radiogenic calculated backwards =', back, 'H_f =', pl.H_f*1e12, 'pW/kg')
        else:
            print('radiogenic calculated backwards =', back, 'H_0 =', pl.H_0 * 1e12, 'pW/kg')
        # print('\nplanet object', pprint(vars(pl)))

    t0 = t0 * 1e9 * p.years2sec  # convert input times to Gyr
    tf = tf * 1e9 * p.years2sec

    f = integrate.solve_ivp(fun=lambda t, y: LHS(t, y, **dict(pl=pl, **run_kwargs, **kwargs)),
                            t_span=(t0, tf), y0=[T_m0, T_c0, D_l0], t_eval=t_eval, max_step=100e6 * p.years2sec,
                            method='RK45', dense_output=False)

    # initial postprocessing to get values of derived parameters across time
    pl.T_m = f.y[0]
    pl.T_c = f.y[1]
    pl.D_l = f.y[2]
    pl.t = f.t
    pl = recalculate(pl.t, pl, verbose=True, **run_kwargs, **kwargs)

    return pl


def LHS(t, y, pl=None, **kwargs):
    """ ODE equation to solve, LHS = 0 """

    # print('t', t*p.sec2Gyr, 'Gyr, T_m', y[0])

    pl.T_m = y[0]
    pl.T_c = y[1]
    pl.D_l = y[2]
    pl = recalculate(t, pl, **kwargs)

    dTdt_c = th.dTdt(-pl.Q_core, pl.M_c, pl.c_c)
    dTdt_m = th.dTdt(-pl.Q_ubl + pl.H_rad_m + pl.Q_core, pl.M_conv, pl.c_m)
    dDdt = th.lid_growth(T_m=pl.T_m, q_ubl=pl.q_ubl, h0=pl.h_rad_m, R_p=pl.R_p, R_l=pl.R_l, T_l=pl.T_l, rho_m=pl.rho_m,
                         T_s=pl.T_s,
                         c_m=pl.c_m, k_m=pl.k_m, **kwargs)

    return [dTdt_m, dTdt_c, dDdt]


def recalculate(t, pl, verbose=False, **kwargs):
    """ calcualte all the thermal convection parameters (single time step or postprocess all timesteps"""

    # if verbose:
    #     print('\nrecalculate kwargs:\n', kwargs)
    #     print('\nplanet:')
    #     for attr in dir(pl):
    #         print("pl.%s = %r" % (attr, getattr(pl, attr)))

    # check if there are any weird issues
    try:
        lid_too_thick = any([list(pl.D_l)]) > pl.R_p - pl.R_c
    except TypeError:
        lid_too_thick = pl.D_l > pl.R_p - pl.R_c
    if lid_too_thick:
        print('{:.2f} M_E: convection stopped at {:.2f} Gyr because lid grew to core'.format(pl.M_p / p.M_E,
                                                                                             t * p.sec2Gyr))
        pl.D_l = pl.R_p - pl.R_c

    # geometry update
    pl.R_l = pl.R_p - pl.D_l
    pl.d_m = pl.R_l - pl.R_c
    V_lid = 4 / 3 * np.pi * (pl.R_p ** 3 - pl.R_l ** 3)
    pl.M_lid = V_lid * pl.rho_m  # not considering density change in lid due to melting/differentiation
    pl.M_conv = pl.M_m - pl.M_lid  # mass of convecting region

    # radiogenic heating
    pl.h_rad_m = th.h_rad(t, H_0=pl.H_0, H_f=pl.H_f, c_n=pl.c_n, p_n=pl.p_n, lambda_n=pl.lambda_n, verbose=verbose, **kwargs)  # W kg^-1
    pl.a0 = pl.h_rad_m * pl.rho_m  # volumetric heating in W m^-3
    pl.H_rad_m = th.H_rad(h=pl.h_rad_m, M=pl.M_conv, **kwargs)  # mantle radiogenic heating in W
    pl.H_rad_lid = th.H_rad(h=pl.h_rad_m, M=pl.M_lid, **kwargs)  # lid radiogenic heating in W

    # temperature update
    pl.T_l = th.T_lid(T_m=pl.T_m, a_rh=pl.a_rh, Ea=pl.Ea)
    pl.dT_rh = pl.a_rh * (p.R_b * pl.T_m ** 2 / pl.Ea)
    pl.dT_m = pl.T_c - pl.T_l
    try:
        bad_T = any([list(pl.dT_m)]) < 0
    except TypeError:
        bad_T = pl.dT_m < 0
    if bad_T:
        print('{:.2f} M_E: SERIOUS WARNING: -ve T_m at {:.2f} Gyr'.format(pl.M_p / p.M_E, t * p.sec2Gyr))

    # viscosity update
    pl.eta_s = rh.dynamic_viscosity(T=pl.T_s, pl=pl, **kwargs)  # viscosity at surface temperature
    pl.eta_l = rh.dynamic_viscosity(T=pl.T_l, pl=pl, **kwargs)  # viscosity at lid base temperature
    pl.eta_m = rh.dynamic_viscosity(T=pl.T_m, pl=pl, **kwargs)  # viscosity at interior mantle temperature
    pl.eta_cmb = rh.dynamic_viscosity(T=(pl.T_c + pl.T_m) / 2, pl=pl, **kwargs)  # viscosity at cmb
    pl.delta_eta = pl.eta_s / rh.dynamic_viscosity(T=(pl.T_c - pl.T_s), pl=pl, **kwargs)  # total viscosity contrast

    # equivalent to effective Ra_i in 2D models with d and dT decreased by lid thickness and temp drop
    pl.Ra_i_eff = th.Ra(eta=pl.eta_m, rho=pl.rho_m, g=pl.g_sfc, deltaT=abs(pl.dT_m), l=pl.d_m, kappa=pl.kappa_m,
                        alpha=pl.alpha_m)

    # the 'overall' interior Ra covering the whole domain from bottom to surface
    pl.Ra_i = th.Ra(eta=pl.eta_m, rho=pl.rho_m, g=pl.g_sfc, kappa=pl.kappa_m, alpha=pl.alpha_m,
                    deltaT=abs(pl.dT_m) + (pl.T_l - pl.T_s),
                    l=pl.R_p - pl.R_c)

    pl.delta_rh = th.bdy_thickness(d_m=pl.d_m, Ra_crit=pl.Ra_crit_u, beta=pl.beta_u, Ra_rh=pl.Ra_i_eff, **kwargs)
    pl.q_ubl = th.q_bl(deltaT=pl.T_m - pl.T_l, k=pl.k_m, d_bl=pl.delta_rh, beta=pl.beta_u, **kwargs)
    pl.Q_ubl = th.Q_bl(q=pl.q_ubl, R=pl.R_l)
    pl.Ra_crit_c = 0.28 * pl.Ra_i_eff ** 0.21  # where is this from?

    # core

    pl.TBL_c = th.bdy_thickness(dT=abs(pl.T_c - pl.dT_m), d_m=pl.d_m, eta_m=pl.eta_cmb, g=pl.g_cmb, Ra_crit=pl.Ra_crit_c,
                                kappa_m=pl.kappa_m, alpha_m=pl.alpha_m, rho_m=pl.rho_m, **kwargs)
    pl.q_core = th.q_bl(deltaT=pl.T_c - pl.T_m, k=pl.k_lm, d_bl=pl.TBL_c, beta=pl.beta_c, **kwargs)
    pl.Q_core = th.Q_bl(q=pl.q_core, R=pl.R_c)

    pl.Ra_F = th.Ra_F(eta=pl.eta_m, kappa=pl.kappa_m, H=pl.h_rad_m, alpha=pl.alpha_m, k=pl.k_m, rho=pl.rho_m,
                      g=pl.g_sfc, l=pl.R_l - pl.R_c, F_b=pl.q_core)

    pl.q_sfc = th.sph_flux(pl.R_p, a0=pl.a0, T_l=pl.T_l, R_l=pl.R_l, k_m=pl.k_m, T_s=pl.T_s, R_p=pl.R_p,
                           **kwargs)  # W m^-2

    pl.Q_sfc = th.Q_bl(q=pl.q_sfc, R=pl.R_p)
    pl.urey = (pl.H_rad_m + pl.H_rad_lid) / pl.Q_sfc

    pl.T_avg = th.T_mean(T_m=pl.T_m, T_l=pl.T_l, R_p=pl.R_p, R_l=pl.R_l, R_c=pl.R_c, T_s=pl.T_s, k_m=pl.k_m,
                         a0=pl.a0, **kwargs)

    return pl


def steadystate(pl, age=4.5/p.sec2Gyr, T_m0=2000, tol=1e-1, verbose=False, **kwargs):
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
    h_t = th.h_rad(t=age, H_0=pl.H_0, c_n=pl.c_n, p_n=pl.p_n, lambda_n=pl.lambda_n, **kwargs)  # rad heating per unit mass
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
        Nu = (rho_m*g*alpha_m*(T_m - T_l)*d_m**3 / (kappa_m * eta_m * Ra_crit))**beta
        T_m1 = M_m/SA_p * h_t * d_m / (Nu * k_m) + T_l
        diff = abs(T_m - T_m1)
        # print('T_m1', T_m1)
    if verbose:
        print('T_m', T_m1, 'K | d_lid', d_lid*1e-3, 'km')

    #save
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
#         Q_ubl = th.Q_bl(q=q_ubl, R=R_l)
#         q_sfc = th.sph_flux(pl.R_p, a0=pl.a0, T_l=T_l, R_l=R_l, k_m=pl.k_m, T_s=pl.T_s, R_p=pl.R_p, **kwargs)
#
#         Q_sfc = th.Q_bl(q=q_sfc, R=pl.R_p)  # W
#
#     # TODO: sympy?



