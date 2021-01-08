###### TOPOGRAPHY ######
from . import parameters as p
import numpy as np

def topography(pl, **kwargs):
    # get all topography parameters for a planet (given its thermal history)
    pl.heuristic_h = pl.delta_rh * pl.dT_rh * pl.alpha_m
    pl.dyn_top_rms = dyn_topo_aspect(pl, **kwargs)
    pl.dyn_top_KH = dyn_topo_KH(pl)
    pl.dyn_top_rms_isoviscous = dyn_topo_Lees(pl)
    
    return pl


def dyn_topo_aspect(pl, **kwargs):
    h_prime = 2.26 * pl.heuristic_h  # fit to chaotic regime
    return h_prime


def dyn_topo_KH(pl=None, Ra_i=None, **kwargs): # Kiefer and Hager 1992 scaling
    if pl is not None:
        Ra_i = pl.Ra_i_eff
    return 66*Ra_i**-0.121
                

def dyn_topo_Lees(pl=None, F=None, rho_m=None, rho_w=0, alpha_m=None, eta_m=None, kappa_m=None, g_sfc=None, l=None,
             k_m=None, C=5.4, deltaT_m = None, **kwargs):
    # root mean square dynamic topography from isoviscous model
    
    if pl is not None:
        rho_m = pl.rho_m
        alpha_m = pl.alpha_m
        F = pl.q_ubl
        eta_m = pl.eta_m
        kappa_m = pl.kappa_m
        g_sfc = pl.g_sfc
        k_m = pl.k_m
        l = pl.d_m
        dT_m = pl.dT_m
        T_m = pl.T_m
        Ea = pl.Ea
        Ra_crit = pl.Ra_crit_u
        dT_rh = p.R_b*T_m**2/Ea
        TBL_u = pl.delta_rh
        a_rh = pl.a_rh
        Ra = pl.Ra_i_eff
    
   # print('original value', C*rho_m/(rho_m-rho_w) * ((alpha_m*F[-1]*eta_m[-1]*kappa_m)/(rho_m*g_sfc*k_m))**(1/2))

    Ra_F = g_sfc * rho_m * alpha_m * l**4 * F / (kappa_m * k_m * eta_m)
    #Ra_F = a_rh*dT_rh/(dT_m*Ra_crit**(1/3)) * Ra**(4/3)   

    RMS = C* eta_m*kappa_m / (rho_m*g_sfc*l**2) * (Ra_F)**(1/2) # eqn 33 Parsons & Daly
    RMS = C*rho_m/(rho_m-rho_w) * ((alpha_m*F*eta_m*kappa_m)/(rho_m*g_sfc*k_m))**(1/2)
    #RMS = C * alpha_m * l * (a_rh*dT_m*dT_rh / Ra_crit**(1/3))**(1/2) * Ra**(-1/3) # Ra version
    return RMS


def convective_stress(pl=None, where='lid', rho_m=None, alpha_m=None, g_sfc=None, T_m=None, Ea=None, k_m=None, q_ubl=None,
                      d_bl=None, **kwargs): 
    # scaling law from Reese+ (2005) - warning: this is for shear stress, don't think you can apply to dyn topo which wants normal stress at the surface minus average normal stress at the surface
    if where is 'lid':
        C = 2.2
    elif (where is 'mantle') or (where is 'interior'):
        C = 0.1
    if pl is not None:
        rho_m = pl.rho_m
        alpha_m = pl.alpha_m
        q_ubl = pl.q_ubl
        g_sfc = pl.g_sfc
        k_m = pl.k_m
        T_m = pl.T_m
        Ea = pl.Ea
        d_bl = pl.delta_rh
    return C* rho_m*alpha_m*g_sfc*(p.R_b*T_m**2/Ea)**2 * k_m/q_ubl
