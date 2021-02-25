import numpy as np
from . import parameters as p


###### SOME TEMPERATURES AND FLUXES ######


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


def T_mean(T_m=None, T_l=None, R_p=None, R_l=None, R_c=None, T_s=None, k_m=None, a0=None, **kwargs):
    '''average temperature across convecting region and lid'''
    c1 = k_m * (T_l - T_s - a0 / (6 * k_m) * (R_p ** 2 - R_l ** 2)) / (R_l ** -1 - R_p ** -1)
    c2 = T_s + a0 / (6 * k_m) * R_p ** 2 - c1 / (k_m * R_p)
    return 3 / (R_p ** 3 - R_c ** 3) * (
                (T_m / 3) * (R_l ** 3 - R_c ** 3) - a0 / (30 * k_m) * (R_p ** 5 - R_l ** 5) + c1 / (2 * k_m) * (
                    R_p ** 2 - R_l ** 2) + c2 / 3 * (R_p ** 3 - R_l ** 3))


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


def bdy_thickness(dT=None, d_m=None, Ra_crit=None, beta=1/3, g=None, Ra_rh=None,
                       kappa_m=None, eta_m=None, alpha_m=None, rho_m=None, **kwargs):
    """Thickness of thermal boundary layer based on critical Ra """
    if Ra_rh is None:
        Ra_rh = Ra(eta=eta_m, kappa=kappa_m, alpha=alpha_m, rho=rho_m, g=g, deltaT=dT, l=d_m)
    return d_m * (Ra_crit / Ra_rh) ** beta


def h_rad(t, tf=None, H_0=None, H_f=None, c_n=None, p_n=None, lambda_n=None, backwards_cooling=False, verbose=False, **kwargs):
    """Calculate radiogenic heating in W kg^-1 from Korenaga (2006)"""
    c_n = np.array(c_n)
    p_n = np.array(p_n)
    lambda_n = np.array(lambda_n)
    x_n = c_n * p_n
    x_tot = np.sum(x_n)
    h_n = x_n / x_tot  # heat produced per kg of isotope normalized to total U (hence a ratio; unitless)

    if backwards_cooling:  # H_f refers to tf (nominally 4.5 Gyr)
        if verbose and H_0 is not None:
            print('WARNING: radiogenic heating calculated backwards from present but H_0 is given')
        try:
            #         h = H_f*sum([x*np.exp(y*(tf-t)) for x, y in zip(h_n, lambda_n)])
            h = H_f * sum(h_n * np.exp(lambda_n * (tf - t)))
        except ValueError:
            # for a list of ages
            h = np.zeros(len(t))
            for ii, t_val in enumerate(t):
                h[ii] = H_f * sum(h_n * np.exp(lambda_n * (tf - t_val)))
    else:  # H_0 is W/kg at t0
        if verbose and H_f is not None:
            print('WARNING: radiogenic heating calculated forwards from t0 but H_f is given')
        try:
            h = H_0 / sum(h_n * np.exp(lambda_n * t))
        except ValueError:
            # for a list of ages
            h = np.zeros(len(t))
            for ii, t_val in enumerate(t):
                h[ii] = H_0 / sum(h_n * np.exp(lambda_n * t_val))

    return h


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


def T_lid(T_m, a_rh=None, Ea=None, **kwargs):
    """ temperature at lid base """
    return T_m - a_rh * (p.R_b * T_m ** 2 / Ea)  # Thiriet+ eq. 15


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

