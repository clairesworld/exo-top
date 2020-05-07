###### TOPOGRAPHY ######
import parameters as p
import numpy as np

def topography(pl, C=None, **kwargs):
    # get all topography parameters for a planet (given its thermal history)
    if C=='a_rh':
        C = pl.a_rh
    
    pl.dyn_top_KH = dyn_topo_KH(pl)
    pl.dyn_top_rms = dyn_topo(pl)
    pl.dyn_top_stress = dyn_topo_stress(pl, C=C)
    
    return pl

def dyn_topo_KH(pl=None, Ra_i=None, **kwargs): # Kiefer and Hager 1992 scaling
    if pl is not None:
        Ra_i = pl.Ra_i
    return 66*Ra_i**-0.121
                

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

def dyn_topo_stress(pl=None, C=1, T_m=None, alpha_m=None, Ea=None, l=None, Ra_crit_u=None, Ra_i=None, R_l=None, R_c=None,
                    k=None, H_rad=None, q_ubl=None, a_rh=None, **kwargs): 
    # equivalent to stress scaling in Reese 04 or equatn 34 in PD83
    # scaling for Arrhenius or FK rheology
    if pl is not None:
        alpha_m = pl.alpha_m
        T_m = pl.T_m
        Ea = pl.Ea
        l = pl.d_m
        k = pl.k_m
        Ra_crit_u = pl.Ra_crit_u
        Ra_i = pl.Ra_i
        R_l = pl.R_l
        R_c = pl.R_c
        H_rad = pl.H_rad_m
        q_ubl = pl.q_ubl
        a_rh = pl.a_rh
        
    #Nu = 0.67 * theta**(4/3) * dT*Ea/(p.R_b*T_i**2)
    #F = k*dT/l*Nu 
    
#     r_t = R_l
#     r_b = R_c
    #F = 1/3*H_rad*r_t*(1 - (r_b**3/r_t**3)) #equation 10, H in W m^-3
    F = q_ubl
    dT_rh = a_rh*p.R_b*T_m**2/Ea
    #return C*alpha_m*dT_rh**2*k/F # equation 32
    return C*alpha_m*dT_rh*l*(Ra_crit_u/Ra_i)**(1/3)

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
        d_bl = pl.TBL_u
    return C* rho_m*alpha_m*g_sfc*(p.R_b*T_m**2/Ea)**2 * k_m/q_ubl

def max_topo(pl=None, s=None, rho_m=3500, rho_w=0, g_sfc=9.8, **kwargs): # dynamic topography based on stress scalings from Reese
    if pl is not None:
        try:
            s = pl.sigma_interior
        except:
            s = convective_stress(pl=pl, **kwargs)
        rho_m = pl.rho_m
        g_sfc = pl.g_sfc
    elif s is None:
        kwargs.update({'where':'lid'}) #  lid stress is what matters for topography
        s = convective_stress(rho_m=rho_m, g_sfc=g_sfc, **kwargs)
    delta_rho = rho_m - rho_w
    return s/(delta_rho*g_sfc)