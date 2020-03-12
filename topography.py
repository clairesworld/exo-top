###### TOPOGRAPHY ######
import parameters as p
import numpy as np

def dyn_topo(pl=None, F=None, rho_m=None, rho_w=0, alpha_m=None, eta_m=None, kappa_m=None, g_sfc=None, 
             k_m=None, C=5.4, **kwargs):
    # root mean square dynamic topography
    
    if pl is not None:
        rho_m = pl.rho_m
        alpha_m = pl.alpha_m
        F = pl.q_sfc
        eta_m = pl.eta_m
        kappa_m = pl.kappa_m
        g_sfc = pl.g_sfc
        k_m = pl.k_m
    RMS = C*rho_m/(rho_m-rho_w) * ((alpha_m*F*eta_m*kappa_m)/(rho_m*g_sfc*k_m))**(1/2) # eqn 33 Parsons & Daly
#     print('F', np.array([F[0], F[-1]])*1e3, 'mW/m2......', 'RMS', np.array([RMS[0], RMS[-1]]), 'm')
    return RMS

def convective_stress(pl=None, where='lid', rho_m=None, alpha_m=None, g_sfc=None, T_m=None, Ea=None, k_m=None, q_ubl=None,
                      **kwargs):
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
    return C* rho_m*alpha_m*g_sfc*(p.R_b*T_m**2/Ea)**2 * k_m/q_ubl

def max_topo(pl=None, s=None, rho_m=3500, rho_w=0, g_sfc=9.8, **kwargs):
    if pl is not None:
        try:
            s = pl.sigma_interior
        except:
            s = convective_stress(pl=pl, **kwargs)
        rho_m = pl.rho_m
        g_sfc = pl.g_sfc
    if s is None:
        kwargs.update({'where':'lid'}) #  lid stress is what matters for topography
        s = convective_stress(rho_m=rho_m, g_sfc=g_sfc, **kwargs)
    delta_rho = rho_m - rho_w
    return s/(delta_rho*g_sfc)