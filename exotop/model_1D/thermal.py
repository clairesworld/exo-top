import numpy as np
from . import parameters as p

###### SOME TEMPERATURES AND FLUXES ######
from .parameters import R_b


def thermal_diffusivity(k, rho, C_p):
    """
    Calculate thermal diffusivity
    
    Parameters
    ----------
    k : Thermal conductivity
    C_p : Specific heat capacity in J K^-1 kg^-1
    rho : density in kg m^-3
    """
    return k / (rho * C_p)


def adiabat(T_0, R=None, g=None, R_p=None, h=None, c_v=None, alpha_m=None, f_d=0.5, **kwargs):
    """
    Calculate adiabatic temperature increase at some fractional depth f_d below top

    Parameters
    ----------
    k : Thermal conductivity
    C_p : Specific heat capacity in J K^-1 kg^-1
    rho : density in kg m^-3
    """
    R_0 = R_p - f_d * h  # depth to avg mantle temp (taken to be midpoint between surface and cmb)
    u = np.exp(-(R - R_0) * alpha_m * g / c_v)  # page 39 driscoll & bercovici 2014
    # print('adiabatic T decrease', u)
    return u * T_0


def sph_conduction(r, k_m=None, T_l=None, T_s=None, R_p=None, R_l=None,
                   a0=None, **kwargs):
    """ calculate temperature at r based on conduction in spherical geometry with total radius R_p, a0: internal heating rate """
    c1 = k_m * (T_l - T_s - a0 / (6 * k_m) * (R_p ** 2 - R_l ** 2)) / (R_l ** -1 - R_p ** -1)
    c2 = T_s + a0 / (6 * k_m) * R_p ** 2 - c1 / (k_m * R_p)
    return -a0 / (6 * k_m) * r ** 2 + c1 / (k_m * r) + c2


def sph_flux(r, a0=None, k_m=None, T_l=None, T_s=None, R_p=None, R_l=None, **kwargs):
    """ calculate heat flux at r based on conduction in spherical geometry with total radius R_p """
    c1 = k_m * (T_l - T_s - a0 / (6 * k_m) * (R_p ** 2 - R_l ** 2)) / (R_l ** -1 - R_p ** -1)
    dTdr = -a0 / (3 * k_m) * r - c1 / (k_m * r ** 2)
    return -k_m * dTdr  # for r>0 in m


def rect_flux(r, a0=None, q0=None, r0=None, **kwargs):
    """ calculate heat flux at r based on conduction in rectangular geometry """
    c0 = q0 - a0 * r0
    return a0 * r + c0


def sfc_flux_Nu(eta=None, kappa=None, alpha=None, k=None, rho=None, g_sfc=None, d=None, deltaT_visc=None, a=0.39,
                beta=0.333, ):
    # Thiriet supp S8
    gamma = 1 + beta
    return a * k * d**(3*beta - 1) * (alpha * rho * g_sfc / (kappa * eta)) ** beta * deltaT_visc



def T_mean_old(T_m=None, T_l=None, R_p=None, R_l=None, R_c=None, T_s=None, k_m=None, a0=None, **kwargs):
    """average temperature across convecting region and lid"""
    c1 = k_m * (T_l - T_s - a0 / (6 * k_m) * (R_p ** 2 - R_l ** 2)) / (R_l ** -1 - R_p ** -1)
    c2 = T_s + a0 / (6 * k_m) * R_p ** 2 - c1 / (k_m * R_p)
    return 3 / (R_p ** 3 - R_c ** 3) * (
            (T_m / 3) * (R_l ** 3 - R_c ** 3) - a0 / (30 * k_m) * (R_p ** 5 - R_l ** 5) + c1 / (2 * k_m) * (
            R_p ** 2 - R_l ** 2) + c2 / 3 * (R_p ** 3 - R_l ** 3))


def T_mean(R_c=None, R_l=None, R_p=None, k_m=None, T_m=None, T_l=None, T_s=None, a0=None, delta_rh=None, **kwargs):
    # weighted by volume of shell at r - average temperature over whole mantle
    # import matplotlib.pyplot as plt
    # convective portion

    if np.size(R_l) == 1:
        z0 = R_l - delta_rh
        T0 = T_m
        w0 = 4/3 * np.pi * (z0 - R_c)**3

        # conductive portion
        z_con = np.linspace(z0, R_p, num=1000)
        T_con = sph_conduction(z_con, k_m=k_m, T_l=T0, T_s=T_s, R_p=R_p, R_l=z0, a0=a0)

        # T_con = np.zeros_like(z_con)
        w_con = np.zeros_like(z_con)
        for ii, r in enumerate(z_con[:-1]):
            # T_con[ii] = sph_conduction(r, k_m=k_m, T_l=T0, T_s=T_s, R_p=R_l, R_l=z0, a0=a0)
            w_con[ii] = 4/3 * np.pi * (z_con[ii + 1] - r)**3

        z = np.insert(z_con, 0, z0)
        T = np.insert(T_con, 0, T0)
        w = np.insert(w_con, 0, w0)

        T_av = np.sum(T*w) / np.sum(w)

    else:
        # iterable
        T_av = []
        if np.size(a0) == 1:
            a0 = np.array([a0]*len(R_l))

        # conductive portion
        for ii in range(len(R_l)):
            z0 = R_l[ii] - delta_rh[ii]
            T0 = T_m[ii]
            w0 = 4 / 3 * np.pi * (z0 - R_c) ** 3

            z_con = np.linspace(z0, R_p, num=1000)
            T_con = sph_conduction(z_con, k_m=k_m, T_l=T0, T_s=T_s, R_p=R_p, R_l=z0, a0=a0[ii])

            # T_con = np.zeros_like(z_con)
            w_con = np.zeros_like(z_con)
            for jj, r in enumerate(z_con[:-1]):
                # T_con[jj] = sph_conduction(r, k_m=k_m, T_l=T0, T_s=T_s, R_p=R_l, R_l=z0, a0=a0)
                w_con[jj] = 4 / 3 * np.pi * (z_con[jj + 1] - r) ** 3

            z = np.insert(z_con, 0, z0)
            T = np.insert(T_con, 0, T0)
            w = np.insert(w_con, 0, w0)
            T_av.append(np.sum(T * w) / np.sum(w))

        T_av = np.array(T_av)

    # print('T av', T_av)
    # plt.plot(z[::-1], T[::-1])
    # plt.show()

    return T_av




#####################################################################
#
#   
#    
#    THERMAL MODEL for stagnant lid adapted from Thiriet+ 2019 
#
#
#
#####################################################################

def Ra(nu=None, eta=None, kappa=None, alpha=None, rho=None, g=None, deltaT=None, l=None):
    """ calculate thermal Ra"""
    if eta is None:
        eta = nu * rho
    return rho * alpha * deltaT * l ** 3 * g / (kappa * eta)


def Ra_F(nu=None, eta=None, kappa=None, H=None, alpha=None, k=None, rho=None, g=None, l=None,
         F_b=None):  # basal heating Ra
    """ calculate flux Ra, H is volumetric heating, F_b is bottom heating in W/m^2 """
    if eta is None:
        eta = nu * rho
    return rho * g * alpha * (F_b + H * l) * l ** 4 / (k * kappa * eta)


def bdy_thickness(dT=None, d_m=None, Ra_crit=None, beta=1 / 3, g=None, Ra_rh=None,
                  kappa_m=None, eta_m=None, alpha_m=None, rho_m=None, **kwargs):
    """Thickness of thermal boundary layer based on critical Ra """
    if Ra_rh is None:
        Ra_rh = Ra(eta=eta_m, kappa=kappa_m, alpha=alpha_m, rho=rho_m, g=g, deltaT=dT, l=d_m)
    return d_m * (Ra_crit / Ra_rh) ** beta


def h_rad_old_old(t, age=4.5):
    """Calculate radiogenic heating in W kg^-1"""
    # Half-lives in years from Dye (2012) in Treatise on Geophys
    t_40K_half = 1.26e9
    t_235U_half = 7.04e8
    t_238U_half = 4.46e9
    t_232Th_half = 1.4e10

    # Heating rates of radioisotopes per mass of isotope in W kg^-1 from Dye (2012) in Treatise on Geophys
    h_40K_0 = 28.47e-6
    h_235U_0 = 568.47e-6
    h_238U_0 = 95.13e-6
    h_232Th_0 = 26.3e-6

    # radioisotope abundances
    X_K = 305  #250  # initial abundance of K in wt ppm (in Treatise on Geophysics, think these are by weight but double check)
    X_U = 16e-3  #2e-2  # initial abundane of U in wt ppm ""
    X_Th = 56e-3  #7e-2  # initial abundance of Th in wt ppm ""

    K_0 = 0.0117e-2  # ratio of 40-K to total K at time 0 (in Treatise on Geophysics)
    U_0_235 = 0.0072  # ratio of 235-U to total U at time 0 (in Treatise on Geophysics)
    U_0_238 = 0.9927  # ratio of 238-U to total U at time 0 (in Treatise on Geophysics)
    Th_0 = 1  # ratio of 232-Th to total Th at time 0 (in Treatise on Geophysics)
    H_0 = np.array([h_40K_0, h_235U_0, h_238U_0, h_232Th_0])
    X_0 = np.array([X_K, X_U, X_U, X_Th])*1e-6
    el_0 = np.array([K_0, U_0_235, U_0_238, Th_0])
    t_half = np.array([t_40K_half, t_235U_half, t_238U_half, t_232Th_half])*p.years2sec

    if np.size(t) > 1:
        t_vect = True
    else:
        t_vect = False

    if not t_vect:
        return sum(H_0 * X_0 * el_0 * np.exp(-np.log(2) * (t - age/p.sec2Gyr) / t_half))
    else:
        # for a list of ages
        h = np.zeros(len(t))
        for ii, val in enumerate(t):
            h[ii] = sum(H_0 * X_0 * el_0 * np.exp(-np.log(2) * (val - age/p.sec2Gyr) / t_half))
        return h


def h_rad_old(t, tf=None, H_0=None, H_f=None, c_n=None, p_n=None, lambda_n=None, backwards_cooling=False, verbose=False,
              **kwargs):
    """Calculate radiogenic heating in W kg^-1 from Korenaga (2006)"""
    c_n = np.array(c_n)
    p_n = np.array(p_n)
    lambda_n = np.array(lambda_n)
    x_n = c_n * p_n
    x_tot = np.sum(x_n)
    h_n = x_n / x_tot  # heat produced per kg of isotope normalized to total U (hence a ratio; unitless)

    if backwards_cooling:  # H_f refers to tf (nominally 4.5 Gyr)
        # print('WARNING: possible bug in backwards cooling (confuses H_0 and H_f)')
        # if verbose and H_0 is not None:
        #     print('WARNING: radiogenic heating calculated backwards from present but H_0 is given')
        #     print('   H_0:', H_0)
        #     print('   H_f:', H_f)
        try:
            #         h = H_f*sum([x*np.exp(y*(tf-t)) for x, y in zip(h_n, lambda_n)])
            h = H_f * sum(h_n * np.exp(lambda_n * (tf - t)))
        except ValueError:
            # for a list of ages
            h = np.zeros(len(t))
            for ii, t_val in enumerate(t):
                h[ii] = H_f * sum(h_n * np.exp(lambda_n * (tf - t_val)))

    else:  # H_0 is W/kg at t0
        # if verbose and H_f is not None:
        #     print('WARNING: radiogenic heating calculated forwards from t0 but H_f is given')
        #     print('   H_0:', H_0)
        #     print('   H_f:', H_f)
        try:
            h = H_0 / sum(h_n * np.exp(lambda_n * t))
        except ValueError:
            # for a list of ages
            h = np.zeros(len(t))
            for ii, t_val in enumerate(t):
                h[ii] = H_0 / sum(h_n * np.exp(lambda_n * t_val))

    return h


def h_rad(t, c_i, h_i, tau_i, age, x_Eu=1):
    """ radiogenic heating in W/kg after Table 1 and eqn 1 in O'Neill+ 2020 (SSR) """
    # order of isotopes (IMPORTANT) is [40K, 238U, 235U, 232Th]
    # convert times to Myr to be consistent with units of tau
    t_Myr = t * p.sec2Gyr * 1e3
    h_K = np.array(c_i[0] * h_i[0] * np.exp((age * 1e3 - t_Myr) * np.log(2) / tau_i[0]))
    try:
        h_UTh = np.array(sum(c_i[1:] * h_i[1:] * np.exp((age * 1e3 - t_Myr) * np.log(2) / tau_i[1:])))
    except ValueError:
        t_Myr = np.vstack((t_Myr, t_Myr, t_Myr))
        c_i = c_i.reshape((4, 1))
        h_i = h_i.reshape((4, 1))
        tau_i = tau_i.reshape((4, 1))
        h_UTh = np.array(np.sum(c_i[1:] * h_i[1:] * np.exp((age * 1e3 - t_Myr) * np.log(2) / tau_i[1:]), axis=0))

    return h_K + x_Eu * h_UTh


def H_rad(h=None, M=None, **kwargs):
    """Calculate energy flux radioisotope decay in W"""
    return h * M


def q_bl(deltaT=None, k=None, d_bl=None, **kwargs):
    """ flux (conductive) across boundary layer """
    return k * deltaT / d_bl


def Q_bl(q=None, k=None, deltaT=None, d_bl=None, R=None, **kwargs):
    """ Calculate energy flux from conduction across thermal bdy layer in W using spherical geometry"""
    SA = 4 * np.pi * R ** 2
    if q is None:
        q = q_bl(deltaT=deltaT, k=k, d_bl=d_bl, **kwargs)
    return SA * q


def T_lid(T_m, a_rh, dT_visc, **kwargs):
    """ temperature at lid base """
    return T_m - a_rh * dT_visc  # Thiriet+ eq. 15


def lid_growth(T_m=None, q_ubl=None, h0=None, R_p=None, R_l=None, T_l=None, rho_m=None, T_s=None,
               c_m=None, k_m=None, **kwargs):
    """ change in lid thickness with spherical geometry, instantaneous adjustment,non-PDE """
    a0 = h0 * rho_m  # radiogenic heating in W/m^3
    c1 = k_m * (T_l - T_s - a0 / (6 * k_m) * (R_p ** 2 - R_l ** 2)) / (R_l ** -1 - R_p ** -1)
    return (-q_ubl + a0 / 3 * R_l + c1 / R_l ** 2) / (rho_m * c_m * (T_m - T_l))


def dTdt(Q=None, M=None, C=None, **kwargs):
    """ temperature change with heat flow
    
    Q : flux balance in W
    M : mass in kg
    C : specific heat in J kg^-1 K^-1
    """
    if M == 0:
        return 0
    else:
        return Q / (M * C)


def q_theta_scaling(theta_FK=None, Ra_i=None, T_m=None, T_s=None, d=None, k=None, c1=0.5):
    return c1*k*(T_m - T_s)/d * theta_FK**(-4/3) * Ra_i**(1/3)


def melt_fraction(cc_slope=None, p_i=None, p_f=None):
    # from Foley and Smye eq 9-11 based on Korenaga+ 2002 solidus parameterisation
    f = (p_i - p_f) / 2 * cc_slope
    return f


def convective_velocity(Ra_i=None, kappa=None, d=None, theta=None, c2=0.05, c4=2 / 3):
    # cited in Foley and Smye eq 13-14 from various scaling laws
    return c2 * kappa / d * (Ra_i / theta) ** c4


def melt_flux(phi_melt=None, v_conv=None, p_i_melt=None, rho_lith=None, g=None, d_m=None, SA_p=None, D_l=None, R_p=None):
    # volumetric melt flux in m3/s, Foley & Smye eq 12
    d_melt = p_i_melt / (rho_lith * g)  # depth at which melting begins
    D_upwell = d_m  # diameter of cylindrical upwelling region assumed to be 1:1 aspect ratio cells
    flow_cylinder = v_conv * (d_melt - D_l) * np.pi * D_upwell
    # return flow_cylinder / (np.pi * (d_m / 2) ** 2) * SA_p * phi_melt  # total melt production rate
    return 17.8 * np.pi * R_p * v_conv * (d_melt - D_l) * phi_melt  # only for Earthlike interior structure


def Q_melt(f_melt=None, rho_melt=None, cp=None, T_m=None, T_s=None, L_melt=None, p_melt_i=None, gamma_melt=None,
           gamma_m=None):
    # Foley & Smye parameterisation assumes all melt transported to surface immediately
    # to turn this off just set rho_melt to 0 I guess
    gamma = gamma_melt - gamma_m
    deltaT_melt = T_m - p_melt_i * gamma - T_s
    return f_melt * rho_melt * (cp * deltaT_melt + L_melt)


def melting(pl, **kwargs):
    # to turn this off set pl.f_melt0 = 0, to calculate melt fractions set pl.f_melt = None
    if pl.f_melt0 == 0:
        pl.Q_melt = 0
        return pl

    pl.p_i = (pl.T_m - 1423) / (120e-9 - pl.gamma_m)  # pressure where melting begins
    pl.p_f = pl.rho_lith * pl.g_sfc * pl.D_l  # pressure where melting ends (base of lid)
    try:
        if pl.p_i < pl.p_f:
            pl.p_i = pl.p_f  # don't consider melting within non-convectiong region
    except ValueError:
        for ii, p in enumerate(pl.p_i):
            if p < pl.p_f[ii]:
                pl.p_i[ii] = pl.p_f[ii]  # don't consider melting within non-convectiong region

    if pl.f_melt0 is None:
        pl.v = convective_velocity(Ra_i=pl.Ra_i, kappa=pl.kappa_m, d=pl.d_m, theta=pl.theta_FK, c2=0.05, c4=2 / 3)
        pl.phi_melt = melt_fraction(p_i=pl.p_i, p_f=pl.p_f, cc_slope=pl.cc_slope)  # melt fraction
        pl.f_melt = melt_flux(phi_melt=pl.phi_melt, v_conv=pl.v, p_i_melt=pl.p_i, rho_lith=pl.rho_lith, g=pl.g_sfc,
                              d_m=pl.d_m, SA_p=pl.SA_p, D_l=pl.D_l, R_p=pl.R_p)  # melt flux volumetric
    else:
        pl.f_melt = pl.f_melt0  # in m3/s

    pl.Q_melt = Q_melt(f_melt=pl.f_melt, rho_melt=pl.rho_melt, cp=pl.c_m, T_m=pl.T_m, T_s=pl.T_s,
                       L_melt=pl.L_melt, p_melt_i=pl.p_i,
                       gamma_melt=pl.gamma_melt,
                       gamma_m=pl.gamma_m)
    return pl


def d_lid_ss(pl=None, Tm=None, a_rh=None, k=None, Ea=None, H0=None, Ra_crit=None, eta_0=None, T_ref=None,
             kappa_m=None, alpha_m=None, g_sfc=None, rho_m=None, Ts=None, **kwargs):
    """ from sympy solution for d in steady state """
    if pl is not None:
        Tm = pl.T_m
        a_rh = pl.a_rh
        k = pl.k_m
        Ea = pl.Ea
        H0 = pl.H_0
        Ra_crit = pl.Ra_crit_u
        eta_0 = pl.eta_0  # reference visc
        T_ref = pl.T_ref
        kappa_m = pl.kappa_m
        alpha_m = pl.alpha_m
        g_sfc = pl.g_sfc
        rho_m = pl.rho_m
        Ts = pl.T_s

    return (-R_b * Tm ** 2 * a_rh * k + np.sqrt(k * (2.0 * Ea ** 2 * H0 * Tm * (
            Ea * Ra_crit * eta_0 * kappa_m * np.exp(Ea / (R_b * Tm) - Ea / (R_b * T_ref)) / (
            R_b * Tm ** 2 * a_rh * alpha_m * g_sfc * rho_m)) ** 0.666666666666667 - 2.0 * Ea ** 2 * H0 * Ts * (
                                                             Ea * Ra_crit * eta_0 * kappa_m * np.exp(
                                                         Ea / (R_b * Tm) - Ea / (R_b * T_ref)) / (
                                                                     R_b * Tm ** 2 * a_rh * alpha_m * g_sfc * rho_m)) ** 0.666666666666667 - 2.0 * Ea * H0 * R_b * Tm ** 2 * a_rh * (
                                                             Ea * Ra_crit * eta_0 * kappa_m * np.exp(
                                                         Ea / (R_b * Tm) - Ea / (R_b * T_ref)) / (
                                                                     R_b * Tm ** 2 * a_rh * alpha_m * g_sfc * rho_m)) ** 0.666666666666667 + R_b ** 2 * Tm ** 4 * a_rh ** 2 * k))) / (
                   Ea * H0 * (Ea * Ra_crit * eta_0 * kappa_m * np.exp(Ea / (R_b * Tm) - Ea / (R_b * T_ref)) / (
                   R_b * Tm ** 2 * a_rh * alpha_m * g_sfc * rho_m)) ** (1 / 3))
