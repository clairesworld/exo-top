###### TOPOGRAPHY ######

def dyn_topo(pl, F=None, rho_m=None, rho_w=0, alpha_m=None, eta_m=None, kappa_m=None, g_sfc=None, 
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
    return C*rho_m/(rho_m-rho_w) * ((alpha_m*F*eta_m*kappa_m)/(rho_m*g_sfc*k_m))**(1/2) # eqn 33 Parsons & Daly