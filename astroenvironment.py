import numpy as np
import parameters as p


###### ASTRONOMY ######

def luminosity(tau, L=None, **kwargs):
    return L*p.L_sun
    
def q_star(t=None, Alb=None, sma=None, **kwargs):
    """Calculate incident stellar flux density (over entire heliocentric sphere with radius a) in W m^-2"""
    return luminosity(t, **kwargs)*(1-Alb)/(4*np.pi*(sma*p.AU2m)**2) # 4pi is the total solid angle subtended by a sphere

def q_sfc_outgoing(R_p=None, SA_p=None, **kwargs):
    return q_star(**kwargs)*(np.pi*R_p**2)/SA_p # assume no geothermal contribution, pi*R^2 cancels out

def T_sfc(q_out=None, **kwargs):
    """Calculate equilibrium surface temperature given outgoing radiation in W m^-2"""
    return (q_out/p.sb)**(1/4)