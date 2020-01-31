import numpy as np
import six
import math

###### PHYSICAL CONSTANTS #####~
M_E = 5.972e24 # earth mass in kg
R_E = 6371e3 # earth radius in m
L_sun =  3.9e26 # solar luminosity in W
years2sec = 31557600
sec2Gyr = 1e-9/years2sec
AU2m = 1.5e11
sb = 5.67e-8 # Stefan Boltzmann constant in W m^-2 K^-4

###### ASTRONOMY ######

def luminosity(tau, L=None, **kwargs):
    return L*L_sun
    
def q_star(t=None, Alb=None, sma=None, **kwargs):
    """Calculate incident stellar flux density (over entire heliocentric sphere with radius a) in W m^-2"""
    return luminosity(t, **kwargs)*(1-Alb)/(4*np.pi*(sma*AU2m)**2) # 4pi is the total solid angle subtended by a sphere

def q_sfc_outgoing(R_p=None, SA_p=None, **kwargs):
    return q_star(**kwargs)*(np.pi*R_p**2)/SA_p # assume no geothermal contribution, pi*R^2 cancels out

def T_sfc(q_out=None, **kwargs):
    """Calculate equilibrium surface temperature given outgoing radiation in W m^-2"""
    return (q_out/sb)**(1/4)