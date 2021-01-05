import numpy as np
from scipy import integrate
import six
import math
from . import parameters as p
# import astroenvironment as ast
# import geometry as geom
from . import rheology as rh
# import terrestrialplanet as tp
from collections.abc import Iterable



###### SOME THERMODYNAMICS ######
def thermal_diffusivity(k, rho, C_p):
    """
    Calculate thermal diffusivity
    
    Parameters
    ----------
    k : Thermal conductivity
    C_p : Specific heat capacity in J K^-1 kg^-1
    rho : density in kg m^-3
    """
    return k/(rho*C_p)

def adiabat(T_0, R=None, g=None, R_p=None, h=None, c_v=None, alpha_m=None, adiabatic=True, **kwargs):
    if adiabatic:
        R_0 = R_p - 0.5*h # depth to avg mantle temp (taken to be midpoint between surface and cmb)
        u = np.exp(-(R - R_0)*alpha_m*g/c_v) # page 39 driscoll & bercovici 2014
        #print('adiabatic T decrease', u)
        return u*T_0
    else:
        return T_0

def Ra(nu=None, eta=None, kappa=None, alpha=None, rho=None, g=None, deltaT=None, l=None):
    if (nu is None) and (eta is not None):
        return rho*alpha*deltaT*l**3*g/(kappa*eta)
    elif (nu is not None) and (eta is None):
        return alpha*deltaT*l**3*g/(kappa*nu)
    
def Ra_F(pl=None, nu=None, eta=None, kappa=None, H=None, alpha=None, k=None, rho=None, g=None, l=None, F_b=None): # basal heating Ra
    # H is volumetric heating, F_b is bottom heating in W/m^2
    if pl is not None:
        eta = pl.eta_m
        kappa = pl.kappa_m
        H = 0#pl.h_rad*pl.rho_m
        alpha = pl.alpha_m
        k = pl.k_m
        rho = pl.rho_m
        g = pl.g_sfc
        l = pl.d_m
        F_b = pl.q_ubl
    elif eta is None:
        eta = nu*rho
    return rho*g*alpha*(F_b + H*l)*l**4 / (k*kappa*eta)

def sph_conduction(r, k_m=None, T_l=None, T_s=None, R_p=None, R_l=None,
                   a0=None, **kwargs):
    c1 = k_m*(T_l - T_s - a0/(6*k_m)*(R_p**2 - R_l**2))/(R_l**-1 - R_p**-1)
    c2 = T_s + a0/(6*k_m)*R_p**2 - c1/(k_m*R_p)
    return -a0/(6*k_m)*r**2 + c1/(k_m*r) + c2

def sph_flux(r, a0=None, k_m=None, T_l=None, T_s=None, R_p=None, R_l=None, **kwargs):
    c1 = k_m*(T_l - T_s - a0/(6*k_m)*(R_p**2 - R_l**2))/(R_l**-1 - R_p**-1)
    dTdr = -a0/(3*k_m)*r - c1/(k_m*r**2)
    return -k_m*dTdr # for r>0 in m

def rect_flux(r, a0=None, q0=None, r0=None, **kwargs):
    c0 = q0 - a0*r0
    return a0*r + c0

def T_mean(T_m=None, T_l=None, R_p=None, R_l=None, R_c=None, T_s=None, k_m=None, a0=None, **kwargs):
    '''average temperature across convecting region and lid'''
    c1 = k_m*(T_l - T_s - a0/(6*k_m)*(R_p**2 - R_l**2))/(R_l**-1 - R_p**-1)
    c2 = T_s + a0/(6*k_m)*R_p**2 - c1/(k_m*R_p)
    return 3/(R_p**3 - R_c**3)*((T_m/3)*(R_l**3 - R_c**3) - a0/(30*k_m)*(R_p**5 - R_l**5) + c1/(2*k_m)*(R_p**2 - R_l**2)  + c2/3*(R_p**3 - R_l**3))


 #####################################################################
#
#   
#    
#    THERMAL MODEL for stagnant lid adapted from Thiriet+ 2019 
#
#
#
 #####################################################################
    
    
def bdy_thickness_beta(dT=None, d_m=None, Ra_crit=None, beta=None, g=None, Ra_rh=None,
                       kappa_m=None, eta_m=None, alpha_m=None, rho_m=None, **kwargs):
    """Thickness of thermal boundary layer """
    if beta is None:
        beta = 1/3

    if Ra_rh is None:
        Ra_rh = alpha_m*rho_m*g*dT*(d_m)**3 / (kappa_m*eta_m)
    return (d_m) * (Ra_crit/Ra_rh)**beta

def inv_bdy_thickness(dT=None, Ra_crit=None, g=None, kappa_m=None, eta_m=None, alpha_m=None, rho_m=None, 
                  **kwargs):
    """Thickness of thermal boundary layer """
    return ((alpha_m*rho_m*g*np.absolute(dT))/(Ra_crit*kappa_m*eta_m))**(1/3)
    
def h_rad(t, tf=None, H_0=None, c_n=None, p_n=None, lambda_n=None, **kwargs):
    """Calculate radiogenic heating in W kg^-1 from Korenaga (2006)"""
    c_n = np.array(c_n)
    p_n = np.array(p_n)
    lambda_n = np.array(lambda_n)
    x_n = c_n*p_n
#     x_n = [x*y for x, y in zip(c_n, p_n)]
    x_tot = np.sum(x_n)
    h_n = x_n/x_tot
    
    try:
#         h = H_0*sum([x*np.exp(y*(tf-t)) for x, y in zip(h_n, lambda_n)])
        h = H_0*sum(h_n*np.exp(lambda_n*(tf-t)))
    except ValueError:
        # for a list of ages
        h = np.zeros(len(t))
        for ii, t_val in enumerate(t):
            h[ii] = H_0*sum(h_n*np.exp(lambda_n*(tf-t_val)))
    return h
        
def H_rad(h=None, M=None, **kwargs):
    """Calculate energy flux radioisotope decay in W"""
    return h*M 

def q_bl(deltaT, k=None, d_bl=None, beta=None, d_bl_inv=None, **kwargs):
#     print('k_m', k)
    if d_bl_inv is None:
        return k*deltaT/d_bl #a_BL*Ra_rh**beta_BL * k*deltaT/h
    else:
        return k*deltaT*d_bl_inv

def Q_bl(q=None, k=None, deltaT=None, d_bl=None, beta=None, R=None, **kwargs):
    """Calculate energy flux from conduction across thermal bdy layer in W""" 
    SA = 4*np.pi*R**2
    if q is None:
        return SA*q_bl(deltaT, k=k, d_bl=d_bl, beta=beta, **kwargs)
    else:
        return SA*q

def T_lid(T_m, a_rh=None, Ea=None, **kwargs):
#     print('a_rh', a_rh)
    return T_m - a_rh*(p.R_b*T_m**2/Ea) # temperature at base of stagnant lid, Thiriet+ eq. 15

def lid_growth(T_m=None, q_ubl=None, h0=None, R_p=None, R_l=None, T_l=None, rho_m=None, T_s=None,
               c_m=None, k_m=None, **kwargs):    
    a0 = h0*rho_m # radiogenic heating in W/m^3
    c1 = k_m*(T_l - T_s - a0/(6*k_m)*(R_p**2 - R_l**2))/(R_l**-1 - R_p**-1)
    return (-q_ubl + a0/3*R_l + c1/R_l**2)/(rho_m*c_m*(T_m - T_l)) # spherical

def dTdt(Q, M, C, **kwargs):
    """ temperature change 
    
    Q : flux balance in W
    M : mass in kg
    C : specific heat in J kg^-1 K^-1
    """
    return Q/(M*C)


def LHS(t, y, pl=None, adiabats=0, complexity=3, Tlid_ini=None, **kwargs):
#     print('kwargs in LHS', kwargs)
    pl.T_m = y[0]
    pl.T_c = y[1]
    pl.D_l = y[2]
    pl = recalculate(t, pl, adiabats=adiabats, complexity=complexity, Tlid_ini=Tlid_ini, **kwargs)
    if pl.SA_c>0:
        dTdt_c = dTdt(-pl.Q_core, pl.M_c, pl.c_c)
    else:
        dTdt_c = 0
    dTdt_m =  dTdt(-pl.Q_ubl + pl.H_rad_m + pl.Q_core, pl.M_conv, pl.c_m)
    dDdt = lid_growth(T_m=pl.T_m, q_ubl=pl.q_ubl, h0=pl.h_rad_m, R_p=pl.R_p, R_l=pl.R_l, T_l=pl.T_l, rho_m=pl.rho_m,T_s=pl.T_s,
                      c_m=pl.c_m, k_m=pl.k_m, **kwargs)

    if complexity==3:
        return [dTdt_m, dTdt_c, dDdt]
    elif (complexity==2) or hasattr(pl, 'D_l_const'):
        return [dTdt_m, dTdt_c, 0]
    elif complexity==1:
        return [dTdt_m, dTdt_c, 0]


def solve(pl, t0=0, tf=4.5, T_m0=1750, T_c0=2250, D_l0=100e3, complexity=3, t_eval=None, verbose=False, **kwargs):
    # scale any model input values as necessary
    if verbose:
        print('Solving 1D thermal history with T_m0 =', T_m0, 'K, T_c0 =', T_c0, 'K, D_l0 =', D_l0, 'm, tf =', tf, 'Gyr')
    t0 = t0*1e9*p.years2sec  # input times in Gyr
    tf = tf*1e9*p.years2sec
    if complexity==1:
        T_c0 = T_m0
    f = integrate.solve_ivp(fun=lambda t, y: LHS(t, y, **dict(pl=pl, tf=tf, complexity=complexity, **kwargs)), 
                            t_span=(t0,tf), y0=[T_m0, T_c0, D_l0], t_eval=t_eval, max_step=100e6*p.years2sec,
                            method='RK45', dense_output=False)
    
    # return planet object with iteratives for evolving variables
    pl.T_m = f.y[0] 
    pl.T_c = f.y[1]
    pl.D_l = f.y[2]
    pl.t = f.t
    pl = recalculate(f.t, pl, tf=tf, **kwargs)
    return pl


def recalculate(t, pl, adiabats=0, complexity=3, Tlid_ini=None, **kwargs):
    if complexity == 1:
        pl.T_c = pl.T_m # no distinction between core and mantle
    pl.h_rad_m = h_rad(t, H_0=pl.H_0, c_n=pl.c_n, p_n=pl.p_n, lambda_n=pl.lambda_n, **kwargs) # W kg^-1
    pl.a0 = pl.h_rad_m*pl.rho_m # radiogenic heating in W m^-3
    # if complexity<3: # lid adjusts instantaneously
    #     pl.D_l = d_lid_ss(Tm=pl.T_m, a_rh=pl.a_rh, k=pl.k_m, Ea=pl.Ea, H0=pl.H_0, Ra_crit=pl.Ra_crit_u, eta_0=pl.eta_0,
    #                       T_ref=pl.T_ref, kappa_m=pl.kappa_m, alpha_m=pl.alpha_m, g_sfc=pl.g_sfc, rho_m=pl.rho_m,
    #                       Ts=pl.T_s)
    if hasattr(pl, 'D_l_const'):
        pl.D_l = [pl.D_l_const]*np.ones_like(t) # keep lid constant for testing purposes
    
    pl.R_l = pl.R_p - pl.D_l
    pl.T_l = T_lid(T_m=pl.T_m, a_rh=pl.a_rh, Ea=pl.Ea)
    pl.T_rh = pl.a_rh*(p.R_b*pl.T_m**2/pl.Ea)
    V_lid = 4/3*np.pi*(pl.R_p**3 - pl.R_l**3)
    pl.M_lid = V_lid*pl.rho_m # should use another density?
    pl.M_conv = pl.M_m - pl.M_lid
    if not hasattr(pl, 'd_m_const'):  # if exploring a constant convecting region height
        pl.d_m = pl.R_l - pl.R_c
    else:
        pl.d_m = [pl.d_m_const]*np.ones_like(t)
    
#     print('mantle eta')
    pl.eta_l = rh.dynamic_viscosity(T=pl.T_l, pl=pl, **kwargs)
    pl.eta_m = rh.dynamic_viscosity(T=pl.T_m, pl=pl, **kwargs)
#     print('cmb eta')
    pl.eta_cmb = rh.dynamic_viscosity(T=(pl.T_c+pl.T_m)/2, pl=pl, **kwargs)
    pl.nu_m = pl.eta_m/pl.rho_m
    pl.nu_cmb = pl.eta_cmb/pl.rho_m
    pl.delta_eta = pl.eta_l - pl.eta_cmb
#     ('upper bl')

    if hasattr(pl, 'dT_m_const'):
        pl.deltaT_m = pl.dT_m_const
    else:
        pl.deltaT_m = pl.T_c - pl.T_l
    if (not isinstance(pl.deltaT_m, Iterable)) and (abs(pl.deltaT_m) < 1e-9):
#         print('catching small dT', pl.T_c, '-', pl.T_l)
        pl.deltaT_m = pl.T_m - pl.T_l
    if (not isinstance(pl.deltaT_m, Iterable)) and (pl.deltaT_m < 0):
        pl.deltaT_m = 1e-5
        
    pl.Ra_i = Ra(eta=pl.eta_m, rho=pl.rho_m, g=pl.g_sfc, deltaT=abs(pl.deltaT_m), l=pl.d_m, kappa=pl.kappa_m, 
                 alpha=pl.alpha_m)
    
    if (not isinstance(pl.Ra_i, Iterable)) and (abs(pl.Ra_i) < 1e-9):
        print('Ra', pl.Ra_i, 'eta', pl.eta_m, 'deltaT', pl.deltaT_m, 'l', pl.d_m, 'T_c - T-l', pl.T_c, '-', pl.T_l)
    elif (not isinstance(pl.Ra_i, Iterable)) and (pl.Ra_i < 0):
        print('NEGATIVE: Ra', pl.Ra_i, 'eta', pl.eta_m, 'deltaT', pl.deltaT_m, 'l', pl.d_m, 'T_c - T-l', pl.T_c, '-', pl.T_l)
            
    pl.TBL_u = bdy_thickness_beta(d_m=pl.d_m, Ra_crit=pl.Ra_crit_u, Ra_rh = pl.Ra_i, beta=pl.beta_u, **kwargs)
    
    if (not isinstance(pl.TBL_u, Iterable)) and (abs(pl.TBL_u) > 1e9):
        print('Ra', pl.Ra_i, 'eta', pl.eta_m, 'deltaT', pl.deltaT_m, 'l', pl.d_m, 'T_c - T-l', pl.T_c, '-', pl.T_l)
    
    pl.q_ubl = q_bl(deltaT=pl.T_m-pl.T_l, k=pl.k_m, d_bl=pl.TBL_u, beta=pl.beta_u, **kwargs)
    pl.Q_ubl = Q_bl(q=pl.q_ubl, R=pl.R_l)
    pl.Ra_crit_c = 0.28*pl.Ra_i**0.21  
    
    if pl.SA_c>0: 
        TBL_c_inv = inv_bdy_thickness(dT=pl.deltaT_m, eta_m=pl.eta_cmb, g=pl.g_cmb, Ra_crit=pl.Ra_crit_c, 
                                      kappa_m=pl.kappa_m, alpha_m=pl.alpha_m, rho_m=pl.rho_m, **kwargs)  
        if hasattr(pl, 'q_core_const'):
            pl.q_core = pl.q_core_const
        else:    
            pl.q_core = q_bl(deltaT=pl.T_c-pl.T_m, k=pl.k_lm, d_bl_inv=TBL_c_inv, beta=pl.beta_c, **kwargs)
        pl.Q_core = Q_bl(q=pl.q_core, R=pl.R_c) 
        pl.TBL_c = TBL_c_inv**-1
    else:
        pl.TBL_c = None
        pl.Q_core = 0
    
    pl.Ra_F = Ra_F(eta=pl.eta_m, kappa=pl.kappa_m, H=pl.h_rad_m, alpha=pl.alpha_m, k=pl.k_m, rho=pl.rho_m, 
                   g=pl.g_sfc, l=pl.R_l - pl.R_c, F_b=pl.q_core)
    
    pl.H_rad_m = H_rad(h=pl.h_rad_m, M=pl.M_conv, **kwargs) # mantle radiogenic heating in W
    pl.H_rad_lid = H_rad(h=pl.h_rad_m, M=pl.M_lid, **kwargs) # lid radiogenic heating in W
    
    pl.q_sfc = sph_flux(pl.R_p, a0=pl.a0, T_l=pl.T_l, R_l=pl.R_l, k_m=pl.k_m, T_s=pl.T_s, R_p=pl.R_p, **kwargs) # W m^-2
    if Tlid_ini=='linear':
        pl.q_sfc = sph_flux(pl.R_p, a0=0, T_l=pl.T_l, R_l=pl.R_l, k_m=pl.k_m, T_s=pl.T_s, R_p=pl.R_p, **kwargs)
    
    pl.Q_sfc = Q_bl(q=pl.q_sfc, R=pl.R_p)
    pl.urey = (pl.H_rad_m + pl.H_rad_lid)/pl.Q_sfc
    
    pl.T_avg = T_mean(T_m=pl.T_m, T_l=pl.T_l, R_p=pl.R_p, R_l=pl.R_l, R_c=pl.R_c, T_s=pl.T_s, k_m=pl.k_m, a0=pl.a0, 
                      **kwargs)
    if Tlid_ini=='linear':
        pl.T_avg = T_mean(T_m=pl.T_m, T_l=pl.T_l, R_p=pl.R_p, R_l=pl.R_l, R_c=pl.R_c, T_s=pl.T_s, k_m=pl.k_m, a0=0, 
                          **kwargs)

    return pl





# def d_lid_ss(Tm, a_rh=None, k=None, Ea=None, H0=None, Ra_crit=None, eta_0=None, T_ref=None,
#           kappa_m=None, alpha_m=None, g_sfc=None, rho_m=None, Ts=None, visc_type=None, **kwargs):
#     """ from sympy solution for d in steady state """
#     R_b = p.R_b
#     if visc_type is 'Thi':
#         return ((-R_b*Tm**2*a_rh*k + np.sqrt(k*(2.0*Ea**2*H0*Tm*(Ea*Ra_crit*eta_0*kappa_m*np.exp(Ea/(R_b*Tm) - Ea/(R_b*T_ref))/(R_b*Tm**2*a_rh*alpha_m*g_sfc*rho_m))**0.666666666666667 - 2.0*Ea**2*H0*Ts*(Ea*Ra_crit*eta_0*kappa_m*np.exp(Ea/(R_b*Tm) - Ea/(R_b*T_ref))/(R_b*Tm**2*a_rh*alpha_m*g_sfc*rho_m))**0.666666666666667 - 2.0*Ea*H0*R_b*Tm**2*a_rh*(Ea*Ra_crit*eta_0*kappa_m*np.exp(Ea/(R_b*Tm) - Ea/(R_b*T_ref))/(R_b*Tm**2*a_rh*alpha_m*g_sfc*rho_m))**0.666666666666667 + R_b**2*Tm**4*a_rh**2*k)))/(Ea*H0*(Ea*Ra_crit*eta_0*kappa_m*np.exp(Ea/(R_b*Tm) - Ea/(R_b*T_ref))/(R_b*Tm**2*a_rh*alpha_m*g_sfc*rho_m))**(1/3)))
#
#     elif visc_type is 'KW':
#         return R_p - 0.111111111111111*(9.0*R_p**2 - 6.0e-14*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(1/3) - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(1/3))/(H_0*rho_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(1/3)))/(-R_p**3 + (-(R_p**2 - 6.66666666666667e-15*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(-0.333333333333333)*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333 - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333)/(H_0*rho_m))**3 + (-R_p**3 + 1.0e-14*R_p*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(-0.333333333333333)*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333 - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333)/(H_0*rho_m) + 0.5*(-6.0*R_p*T_l*k_m + 6.0*R_p*T_s*k_m)/(H_0*rho_m))**2)**0.5 + 1.0e-14*R_p*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(-0.333333333333333)*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333 - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333)/(H_0*rho_m) + 0.5*(-6.0*R_p*T_l*k_m + 6.0*R_p*T_s*k_m)/(H_0*rho_m))**(1/3) - 1.0*(-R_p**3 + (-(R_p**2 - 6.66666666666667e-15*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(-0.333333333333333)*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333 - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333)/(H_0*rho_m))**3 + (-R_p**3 + 1.0e-14*R_p*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(-0.333333333333333)*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333 - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333)/(H_0*rho_m) + 0.5*(-6.0*R_p*T_l*k_m + 6.0*R_p*T_s*k_m)/(H_0*rho_m))**2)**0.5 + 1.0e-14*R_p*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(-0.333333333333333)*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333 - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333)/(H_0*rho_m) + 0.5*(-6.0*R_p*T_l*k_m + 6.0*R_p*T_s*k_m)/(H_0*rho_m))**(1/3), R_p - 0.111111111111111*(-0.5 - 0.866025403784439*I)*(9.0*R_p**2 - 6.0e-14*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(1/3) - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(1/3))/(H_0*rho_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(1/3)))/(-R_p**3 + (-(R_p**2 - 6.66666666666667e-15*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(-0.333333333333333)*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333 - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333)/(H_0*rho_m))**3 + (-R_p**3 + 1.0e-14*R_p*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(-0.333333333333333)*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333 - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333)/(H_0*rho_m) + 0.5*(-6.0*R_p*T_l*k_m + 6.0*R_p*T_s*k_m)/(H_0*rho_m))**2)**0.5 + 1.0e-14*R_p*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(-0.333333333333333)*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333 - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333)/(H_0*rho_m) + 0.5*(-6.0*R_p*T_l*k_m + 6.0*R_p*T_s*k_m)/(H_0*rho_m))**(1/3) - 1.0*(-0.5 + 0.866025403784439*I)*(-R_p**3 + (-(R_p**2 - 6.66666666666667e-15*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(-0.333333333333333)*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333 - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333)/(H_0*rho_m))**3 + (-R_p**3 + 1.0e-14*R_p*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(-0.333333333333333)*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333 - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333)/(H_0*rho_m) + 0.5*(-6.0*R_p*T_l*k_m + 6.0*R_p*T_s*k_m)/(H_0*rho_m))**2)**0.5 + 1.0e-14*R_p*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(-0.333333333333333)*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333 - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333)/(H_0*rho_m) + 0.5*(-6.0*R_p*T_l*k_m + 6.0*R_p*T_s*k_m)/(H_0*rho_m))**(1/3), R_p - 0.111111111111111*(-0.5 + 0.866025403784439*I)*(9.0*R_p**2 - 6.0e-14*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(1/3) - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(1/3))/(H_0*rho_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(1/3)))/(-R_p**3 + (-(R_p**2 - 6.66666666666667e-15*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(-0.333333333333333)*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333 - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333)/(H_0*rho_m))**3 + (-R_p**3 + 1.0e-14*R_p*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(-0.333333333333333)*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333 - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333)/(H_0*rho_m) + 0.5*(-6.0*R_p*T_l*k_m + 6.0*R_p*T_s*k_m)/(H_0*rho_m))**2)**0.5 + 1.0e-14*R_p*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(-0.333333333333333)*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333 - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333)/(H_0*rho_m) + 0.5*(-6.0*R_p*T_l*k_m + 6.0*R_p*T_s*k_m)/(H_0*rho_m))**(1/3) - 1.0*(-0.5 - 0.866025403784439*I)*(-R_p**3 + (-(R_p**2 - 6.66666666666667e-15*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(-0.333333333333333)*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333 - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333)/(H_0*rho_m))**3 + (-R_p**3 + 1.0e-14*R_p*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(-0.333333333333333)*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333 - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333)/(H_0*rho_m) + 0.5*(-6.0*R_p*T_l*k_m + 6.0*R_p*T_s*k_m)/(H_0*rho_m))**2)**0.5 + 1.0e-14*R_p*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**(-0.333333333333333)*(-377976314968461.0*R_p*T_l*k_m + 377976314968461.0*R_p*T_m*k_m + 300000000000000.0*T_l*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333 - 300000000000000.0*T_s*k_m*(Ra_crit*kappa_m*mu*(h_rh/B_rh)**m_rh*exp(Ea/(R_b*T_m))/(A_rh*T_c*alpha_m*g_sfc*rho_m - A_rh*T_m*alpha_m*g_sfc*rho_m))**0.333333333333333)/(H_0*rho_m) + 0.5*(-6.0*R_p*T_l*k_m + 6.0*R_p*T_s*k_m)/(H_0*rho_m))**(1/3)
#
