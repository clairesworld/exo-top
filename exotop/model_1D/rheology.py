import numpy as np
from . import parameters as p


def dynamic_viscosity(T=None, pl=None, visc_type=None, **kwargs):  # note kwargs can contain planet object
    # if np.size(T)> 1:
    #     print('using', visc_type, 'rheology')
    if visc_type == 'constant':
        return pl.nu_0 * pl.rho_m
    elif visc_type == 'Dorn':
        return nu_Dorn(T, pl=pl, **kwargs) * pl.rho_m
    elif visc_type == 'KW':
        return eta_KW(T, pl=pl, **kwargs)
    elif visc_type == 'Driscoll':
        return nu_Driscoll(T, pl=pl, **kwargs) * pl.rho_m
    elif visc_type == 'Thi':
        return eta_Thi(T, pl=pl, **kwargs)
    elif visc_type == 'FK':
        return eta_FK(T, T_i=T, pl=pl, **kwargs)
    elif visc_type == 'Nim':
        return eta_Nim(T, eta_0=1e21, zeta=1e-2, T_ref=kwargs['T_ref'])
    else:
        print('missing viscosity method')


def nu_Driscoll(T, pl=None, nu_0=7e7, Ea=3e5, **kwargs):
    """kinematic viscosity (upper mantle) from eqn 6 in Driscoll & Bercovici"""
    if pl is not None:
        Ea = pl.Ea
        nu_0 = pl.nu_0
    return nu_0 * np.exp(Ea / (p.R_b * T))


def nu_Dorn(T, pl=None, nu_0=1.6e20, Ea=300e3, T_0=1800, **kwargs):
    # viscosity (below lithosphere) from Dorn, Noack & Rozal 2018
    if pl is not None:
        Ea = pl.Ea
        nu_0 = pl.nu_0
        T_0 = pl.T_0
    return nu_0 * np.exp(Ea / p.R_b * (T ** -1 - T_0 ** -1))


def eta_KW(T, pl=None, P=0, Ea=300e3, V_rh=6e-6, mu=80e9, A_rh=8.7e15, h_rh=2.07e-3, B_rh=0.5e-9, m_rh=2.5,
           eta_pre=None, **kwargs):
    """ Karato & Wu 1993, defaults are diffusion creep for dry olivine 
    V = activation volume
    mu = shear modulus
    A = preexponential factor
    B = Burgers vector
    h = grain size
    m = grain size exponent
    """
    if pl is not None:
        Ea = pl.Ea
        V_rh = pl.V_rh
        mu = pl.mu
        A_rh = pl.A_rh
        B_rh = pl.B_rh
        m_rh = pl.m_rh
        h_rh = pl.h_rh
        eta_pre = pl.eta_pre
    Q = Ea + P * V_rh  # activation enthalpy

    if eta_pre is None:
        b = mu / (2 * A_rh) * (h_rh / B_rh) ** m_rh
    else:
        b = eta_pre

    return b * np.exp(Q / (p.R_b * T))


def eta_Thi(T, pl=None, eta_0=1e21, T_ref=1600, Ea=300e3, **kwargs):  # diffusion creep, dry rheology (Thiriet+ 2019)
    if pl is not None:
        Ea = pl.Ea
        eta_0 = pl.eta_0
        T_ref = pl.T_ref
    return eta_0 * np.exp(Ea / p.R_b * (T ** -1 - T_ref ** -1))


def eta_Nim(T, eta_0=1e21, zeta=1e-2, T_ref=1573):  # Nimmo+ 2003 like FK
    return eta_0 * np.exp(-zeta*(T - T_ref))


def eta_FK(T, pl=None, T_s=None, eta_s=None, T_i=None, Ea=300e3, P=0, V_rh=6e-6, **kwargs):
    # Frank-Kamenetskii approximation
    """ takes surface temperature and viscosity, T_i is temperature below the stagnant lid
        Q is technically the activation enthalpy but assume it is constant and equal to activation energy (p=0)
    """
    if pl is not None:
        Ea = pl.Ea
        V_rh = pl.V_rh
        T_s = pl.T_s
        eta_s = pl.eta_s
    Q = Ea + P * V_rh
    gamma = Q / (p.R_b * T_i ** 2)
    return eta_s * np.exp(-gamma * (T - T_s))


def theta_FK(Ea, T_i, T_s, R_b=8.314):
    # Foley & Smye eqn 3
    return Ea * (T_i - T_s) / (R_b*T_i**2)


def grain_size(A_rh, eta_0, Ea, T_ref, m, B_rh, mu):
    # calculate grain size to correspond to linearization
    return B_rh * (2 * A_rh * eta_0 * np.exp(-Ea / (p.R_b * T_ref)) / mu) ** (1 / m)


def visc_factor(pl=None, T_i=None, Ea=None, dT=None, which_Ti='T_m', **kwargs):
    if pl is not None:
        if which_Ti == 'T_m':
            T_i = pl.T_m
        elif which_Ti == 'T_ubl':
            T_i = (pl.T_l + pl.T_m) / 2
        Ea = pl.Ea
        dT = pl.delta_T
    return dT / (p.R_b * T_i ** 2 / Ea)
