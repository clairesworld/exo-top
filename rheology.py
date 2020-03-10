import numpy as np
import parameters as p

def nu_Driscoll(T, nu_0=7e7, Ea=3e5, **kwargs):
    """kinematic viscosity (upper mantle) from eqn 6 in Driscoll & Bercovici"""
    return nu_0*np.exp(Ea/(p.R_b*T))/10

def nu_Dorn(T, nu_0=1.6e20, Ea=300e3, T_0=1800, **kwargs):
    # viscosity (below lithosphere) from Dorn, Noack & Rozal 2018
    return nu_0*np.exp(Ea/p.R_b*(T**-1-T_0**-1))

def eta_KW(T, pl=None, P=0, Ea=300e3, V_rh=6e-6, mu=80e9, A_rh=8.7e15, h_rh=10**-2, B_rh=0.5e-9, m_rh=2.5, **kwargs): 
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
    Q = Ea + P*V_rh # activation enthalpy
    b = mu/(2*A_rh) * (h_rh/B_rh)**m_rh
    return b*np.exp(Q/(p.R_b*T))

def eta_Thi(T, eta_0=1e21, T_ref=1600, Ea=300e3, **kwargs): # diffusion creep, dry rheology (Thiriet+ 2019)
    return eta_0*np.exp(Ea/p.R_b*(T**-1 - T_ref**-1))

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
    Q = Ea + P*V_rh
    gamma = Q/(p.R_b*T_i**2)
    return eta_s*np.exp(-gamma*(T - T_s))

def dynamic_viscosity(T=None, visc_type=None, **kwargs):
    if visc_type=='constant':
        return nu_0*rho_m
    elif visc_type=='Dorn':
        return nu_Dorn(T, **kwargs)*rho_m
    elif visc_type=='KW':
        return eta_KW(T, **kwargs)
    elif visc_type=='Driscoll':
        return nu_Driscoll(T, **kwargs)*rho_m
    elif visc_type=='Thi':
        return eta_Thi(T, **kwargs)
    elif visc_type=='FK':
        return eta_FK(T, T_i=T, **kwargs)