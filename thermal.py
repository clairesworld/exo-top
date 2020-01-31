import numpy as np
import six
import math
import astroenv as astro
import structure 

###### PHYSICAL CONSTANTS #####~
M_E = 5.972e24 # earth mass in kg
R_E = 6371e3 # earth radius in m
years2sec = 31557600
sec2Gyr = 1e-9/years2sec
AU2m = 1.5e11
R_b = 8.3144598 # universal gas constant in J mol −1 K −1

# Decay constant in 1/sec from Dye (2012) in Treatise on Geophys
lambda_n = np.array([0.15541417, 0.98458406, 0.04951051, 0.55011681])*1e-9/years2sec  # [238U, 235U, 232Th, 40K]

# Heating rates of radioisotopes per mass of isotope in W kg^-1 from Dye (2012) in Treatise on Geophys
p_n = [95.13e-6, 568.47e-6, 26.3e-6, 28.47e-6] # [238U, 235U, 232Th, 40K]

# radioisotope abundances
# TODO: are the below values in moles or mass? 
K_0 = 0.0117e-2 # ratio of 40-K to total K at time 0 (in Treatise on Geophysics)
U_0_235 = 0.0072 # ratio of 235-U to total U at time 0 (in Treatise on Geophysics)
U_0_238 = 0.9927 # ratio of 238-U to total U at time 0 (in Treatise on Geophysics)
Th_0 = 1 # ratio of 232-Th to total Th at time 0 (in Treatise on Geophysics)




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

def nu_Driscoll(T, nu_0=7e7, Ea=3e5, **kwargs):
    """kinematic viscosity (upper mantle) from eqn 6 in Driscoll & Bercovici"""
    return nu_0*np.exp(Ea/(R_b*T))/10

def nu_Dorn(T, nu_0=1.6e20, Ea=300e3, T_0=1800, **kwargs):
    # viscosity (below lithosphere) from Dorn, Noack & Rozal 2018
    return nu_0*np.exp(Ea/R_b*(T**-1-T_0**-1))

def nu_KW(T, p=0, **kwargs): # Karato & Wu 1993, diffusion law for dry olivine
    return 2.6e10*np.exp((3e5 + (6e3*p))/(R_b*T))   

def eta_Thi(T, eta_0=1e21, T_ref=1600, Ea=300e3, **kwargs): # diffusion creep, dry rheology (Thiriet+ 2019)
    return eta_0*np.exp(Ea/R_b*(T**-1 - T_ref**-1))

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
    
def d_lid_ss(Tm, a_rh=None, k=None, Ea=None, H0=None, Ra_crit=None, eta_0=None, T_ref=None, 
          kappa_m=None, alpha_m=None, g_sfc=None, rho_m=None, Ts=None, **kwargs):
    """ from sympy solution for d in steady state """
    return ((-R_b*Tm**2*a_rh*k + np.sqrt(k*(2.0*Ea**2*H0*Tm*(Ea*Ra_crit*eta_0*kappa_m*np.exp(Ea/(R_b*Tm) - Ea/(R_b*T_ref))/(R_b*Tm**2*a_rh*alpha_m*g_sfc*rho_m))**0.666666666666667 - 2.0*Ea**2*H0*Ts*(Ea*Ra_crit*eta_0*kappa_m*np.exp(Ea/(R_b*Tm) - Ea/(R_b*T_ref))/(R_b*Tm**2*a_rh*alpha_m*g_sfc*rho_m))**0.666666666666667 - 2.0*Ea*H0*R_b*Tm**2*a_rh*(Ea*Ra_crit*eta_0*kappa_m*np.exp(Ea/(R_b*Tm) - Ea/(R_b*T_ref))/(R_b*Tm**2*a_rh*alpha_m*g_sfc*rho_m))**0.666666666666667 + R_b**2*Tm**4*a_rh**2*k)))/(Ea*H0*(Ea*Ra_crit*eta_0*kappa_m*np.exp(Ea/(R_b*Tm) - Ea/(R_b*T_ref))/(R_b*Tm**2*a_rh*alpha_m*g_sfc*rho_m))**(1/3)))
    
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









###### THERMAL MODEL for stagnant lid adapted from Thiriet+ 2019 ######

def init(k_m=None, rho_m=None, rho_c=None, c_m=None, CMF=None, M_p=None, R_p0=None, R_c0=None, 
         T_s=None, **kwargs):
    if R_p0 is None:
        R_p = radius_zeng(M_p, CMF)*R_E # in m
    else:
        R_p = R_p0
    if R_c0 is None:
        M_m = M_p*(1 - CMF) # mass of mantle
        CRF = CMF**0.5 # Zeng & Jacobsen 2017
        M_c = M_p*CMF # M_p - M_m
        R_c = R_p*CRF
        #R_c = radius_seager(M_p*CMF, CMF=0, m1=4.34*M_E, r1=2.23*R_E) # EoS for iron... is this consistent?
    else:
        R_c = R_c0
        CRF = R_c0/R_p0
        M_c = 4/3*np.pi*R_c**3 * rho_c
        CMF = M_c/M_p
        M_m = 4/3*np.pi*(R_p**3 - R_c**3)*rho_m  #M_p - M_c
    SA_p = SA(R=R_p)
    SA_c = SA(R=R_c) # core surface area 
    h = R_p - R_c # mantle thickness
    g_sfc = grav(M_p, R_p)
    if CMF>0:
        g_cmb = grav(M_c, R_c)
    else:
        g_cmb = 0
    kappa_m = thermal_diffusivity(k_m, rho_m, c_m)
    
    if T_s is None:
        q_out = q_sfc_outgoing(R_p=R_p, SA_p=SA_p, **kwargs)
        T_s = T_sfc(q_out)
    return dict(kappa_m=kappa_m, SA_p=SA_p, SA_c=SA_c, M_m=M_m, g_sfc=g_sfc, R_p=R_p, R_c=R_c,
                CRF=CRF, g_cmb=g_cmb, h=h, M_c=M_c, T_s=T_s)

def dyn_visc(T=None, nu_0=None, visc_type=None, rho_m=None, **kwargs):
    if visc_type=='constant':
        return nu_0*rho_m
    elif visc_type=='Dorn':
        return nu_Dorn(T, **kwargs)*rho_m
    elif visc_type=='KW':
        return nu_KW(T, **kwargs)*rho_m
    elif visc_type=='Driscoll':
        return nu_Driscoll(T, **kwargs)*rho_m
    elif visc_type=='Thi':
        return eta_Thi(T, **kwargs)

def bdy_thickness_beta(dT=None, R_l=None, R_c=None, Ra_crit=None, beta=None, g=None,
                       kappa_m=None, eta_m=None, Ea=None, alpha_m=None, rho_m=None, a_rh=None,**kwargs):
    """Thickness of thermal boundary layer """
    if beta is None:
        beta = 1/3
    Ra_rh = alpha_m*rho_m*g*dT*(R_l - R_c)**3 / (kappa_m*eta_m)
    return (R_l - R_c) * (Ra_crit/Ra_rh)**beta
    
def bdy_thickness(dT=None, Ra_crit=None, g=None, kappa_m=None, eta_m=None, alpha_m=None, rho_m=None, 
                  **kwargs):
    """Thickness of thermal boundary layer """
    return (Ra_crit*kappa_m*eta_m/(alpha_m*rho_m*g*dT))**(1/3)

def inv_bdy_thickness(dT=None, Ra_crit=None, g=None, kappa_m=None, eta_m=None, alpha_m=None, rho_m=None, 
                  **kwargs):
    """Thickness of thermal boundary layer """
    return ((alpha_m*rho_m*g*np.absolute(dT))/(Ra_crit*kappa_m*eta_m))**(1/3)
    
def h_rad(t, tf=None, H_0=None, c_n=None, p_n=None, lambda_n=None, t_vect=False, **kwargs):
    """Calculate radiogenic heating in W kg^-1 from Korenaga (2006)"""
    x_n = np.array(c_n)*np.array(p_n)
    x_tot = sum(x_n)
    h_n = x_n/x_tot
    
    try:
        h = H_0*sum(h_n*np.exp(lambda_n*(tf-t)))
    except ValueError:
        # for a list of ages
        h = np.zeros(len(t))
        for ii, t_val in enumerate(t):
            h[ii] = H_0*sum(h_n*np.exp(lambda_n*(tf-t_val)))
    return h
        
def H_rad(t=None, M=None, **kwargs):
    """Calculate energy flux radioisotope decay in W"""
    h = h_rad(t, **kwargs)
    return h*M 

def q_bl(deltaT, k=None, d_bl=None, beta=None, d_bl_inv=None, **kwargs):
    if d_bl_inv is None:
        return k*deltaT/d_bl #a_BL*Ra_rh**beta_BL * k*deltaT/h
    else:
        return k*deltaT*d_bl_inv

def Q_bl(q=None, deltaT=None, SA=None, h=None, d_bl=None,
         adiabatic=False, beta=None, Ra_rh=None, **kwargs):
    """Calculate energy flux from conduction across thermal bdy layer in W""" 
    if q is None:
        return SA*q_bl(deltaT, k=k, d_bl=d_bl, beta=beta, **kwargs)
    else:
        return SA*q

def T_lid(T_m, a_rh=None, Ea=None, **kwargs):
    return T_m - a_rh*(R_b*T_m**2/Ea) # temperature at base of stagnant lid, Thiriet+ eq. 15

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

def LHS(t, y, M_c=None, c_m=None, c_c=None, SA_c=None, adiabats=0, complexity=3, **kwargs):
    T_m = y[0]
    T_c = y[1]
    D_l = y[2]
    [R_l, T_l, eta_m, eta_cmb, nu_m, nu_cmb, TBL_u, h_rad_m, q_ubl, Q_ubl, Ra_i, Ra_crit_c, TBL_c, 
     q_core, Q_core, M_mi, M_lid, H_rad_m, a0, 
     H_rad_lid, q_sfc, Q_sfc, D_l, T_avg, T_c] = outputs(t, T_m=T_m, T_c=T_c, D_l=D_l, M_c=M_c, c_m=c_m, 
                                                         SA_c=SA_c, complexity=complexity, **kwargs)
    if SA_c>0:
        dTdt_c = dTdt(-Q_core, M_c, c_c)
    else:
        dTdt_c = 0
    dTdt_m =  dTdt(-Q_ubl + H_rad_m + Q_core, M_mi, c_m)
    dDdt = lid_growth(T_m=T_m, q_ubl=q_ubl, h0=h_rad_m, R_l=R_l, T_l=T_l, c_m=c_m, **kwargs)
    if complexity==3:
        return [dTdt_m, dTdt_c, dDdt]
    elif complexity==2:
        return [dTdt_m, dTdt_c, 0]
    elif complexity==1:
        return [dTdt_m, dTdt_c, 0]

def solve(t0=0, tf=None, T_m0=None, T_c0=None, D_l0=None, complexity=3, **kwargs):
    kwargs = {**kwargs, **init(complexity=complexity, **kwargs)} # get derived parameters
    if complexity==1:
        T_c0 = T_m0
    f = integrate.solve_ivp(fun=lambda t, y: LHS(t, y, **dict(tf=tf*1e9*years2sec, complexity=complexity, **kwargs)), 
                            t_span=(t0*1e9*years2sec,tf*1e9*years2sec), y0=[T_m0, T_c0, D_l0], max_step=100e6*years2sec,
                            method='RK45', dense_output=False)
    return f

def outputs(t=None, T_m=None, T_c=None, D_l=None, T_s=None, M_m=None, M_c=None, c_m=None, c_c=None, 
            SA_p=None, SA_c=None, g_sfc=None, g_cmb=None, R_p=None, R_c=None, k_m=None, k_lm=None, 
            Ra_crit_u=None, beta_u=None, beta_c=None, rho_lith=None, adiabats=0, complexity=3, Tlid_ini=None,
            **kwargs):
    if complexity == 1:
        T_c = T_m
    h_rad_m = h_rad(t, **kwargs) # W kg^-1
    a0 = h_rad_m*rho_m # radiogenic heating in W m^-3
    if complexity<3: # lid adjusts instantaneously
        D_l = d_lid_ss(T_m, k=k_m, H0=a0, Ra_crit=Ra_crit_u, g_sfc=g_sfc, Ts=T_s, **kwargs)
    
    R_l = R_p - D_l
    T_l = T_lid(T_m, **kwargs)
    V_lid = 4/3*np.pi*(R_p**3 - R_l**3)
    M_lid = V_lid*rho_m # should use another density?
    M_mi = M_m - M_lid
    if (beta_u is None) or (beta_c is None):
        h = None
    else:
        pass
    eta_m = dyn_visc(T=T_m, **kwargs)
    eta_cmb = dyn_visc(T=(T_c+T_m)/2, **kwargs)
    nu_m = eta_m/rho_m
    nu_cmb = eta_cmb/rho_m
    TBL_u = bdy_thickness_beta(dT=T_c-T_l, eta_m=eta_m, g=g_sfc, Ra_crit=Ra_crit_u, 
                               R_l=R_l, R_c=R_c, beta=beta_u, **kwargs)
    #TBL_u = bdy_thickness(dT=T_c-T_l, eta_m=eta_m, g=g_sfc, Ra_crit=Ra_crit_u, **kwargs)
    q_ubl = q_bl(deltaT=T_m-T_l, k=k_m, d_bl=TBL_u, beta=beta_u, **kwargs)
    Q_ubl = Q_bl(q_ubl, SA=4*np.pi*R_l**2) 
    Ra_i = Ra(eta=eta_m, kappa=kwargs['kappa_m'], alpha=kwargs['alpha_m'], rho=rho_m, 
              g=g_sfc, deltaT=T_m-T_l, l=R_l-R_c)
    Ra_crit_c = 0.28*Ra_i**0.21  
    
    if SA_c>0:
        TBL_c_inv = inv_bdy_thickness(dT=T_c-T_m, eta_m=eta_cmb, g=g_cmb, Ra_crit=Ra_crit_c, **kwargs)  
        q_core = q_bl(deltaT=T_c-T_m, k=k_lm, d_bl_inv=TBL_c_inv, beta=beta_c, **kwargs)
        Q_core = Q_bl(q_core, SA=SA_c) 
        TBL_c = TBL_c_inv**-1
    else:
        TBL_c=None
        Q_core = 0
        
    H_rad_m = H_rad(t, M=M_mi, **kwargs) # mantle radiogenic heating in W
    H_rad_lid = H_rad(t, M=M_lid, **kwargs) # lid radiogenic heating in W
    
    #q_sfc = sfc_flux(q_bl=q_ubl, R_p=R_p, R_l=R_l, m=2, **kwargs)
    q_sfc = sph_flux(R_p, a0=a0, T_l=T_l, R_l=R_l, k_m=k_m, T_s=T_s, R_p=R_p, **kwargs) # surface heat flux in W m^-2
    if Tlid_ini=='linear':
        q_sfc = sph_flux(R_p, a0=0, T_l=T_l, R_l=R_l, k_m=k_m, T_s=T_s, R_p=R_p, **kwargs)
    
    Q_sfc = q_sfc*4*np.pi*R_p**2
    
    T_avg = T_mean(T_m=T_m, T_l=T_l, R_p=R_p, R_l=R_l, R_c=R_c, T_s=T_s, k_m=k_m, a0=a0, **kwargs)
    if Tlid_ini=='linear':
        T_avg = T_mean(T_m=T_m, T_l=T_l, R_p=R_p, R_l=R_l, R_c=R_c, T_s=T_s, k_m=k_m, a0=0, **kwargs)
    
    return R_l, T_l, eta_m, eta_cmb, nu_m, nu_cmb, TBL_u, h_rad_m, q_ubl, Q_ubl, Ra_i, Ra_crit_c, TBL_c, q_core, Q_core, M_mi, M_lid, H_rad_m, a0, H_rad_lid, q_sfc, Q_sfc, D_l, T_avg, T_c 
