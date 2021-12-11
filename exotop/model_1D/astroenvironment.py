import numpy as np
from . import parameters as p


###### ASTRONOMY ######

def luminosity(tau, L=None, **kwargs):
    return L*p.L_sun
    
def q_star(t=None, Alb=None, sma=None, **kwargs):
    """Calculate incident stellar flux density (over entire heliocentric sphere with radius a) in W m^-2"""
    return luminosity(t, **kwargs)*(1-Alb)/(4*np.pi*(sma*p.AU2m)**2) # 4pi is the total solid angle subtended by a sphere

def q_sfc_outgoing(R_p=None, SA_p=None, **kwargs):
    return q_star(**kwargs)*(np.pi*R_p**2)/SA_p # assume no geothermal contribution, pi*R_b^2 cancels out

def T_sfc(q_out=None, **kwargs):
    """Calculate equilibrium surface temperature given outgoing radiation in W m^-2"""
    return (q_out/p.sb)**(1/4)

def radius_seager(M_p, CMF=0.3, k1=-0.20945, k2=0.0804, k3=0.394, m1=None, r1=None):
    if (m1 is None) and (r1 is None):
        if CMF==0.3:
            m1 = 6.41*p.M_E
            r1 = 2.84*p.R_E
        elif CMF==0.675:
            m1 = 6.41*p.M_E
            r1 = 3.19*p.R_E
        elif CMF==0: # all perovskite
            m1 = 7.38*p.M_E
            r1 = 3.58*p.R_E
    M_s = M_p/m1
    R_s = 10**(k1 + 1/3*np.log10(M_s) - k2*M_s**k3)
    return R_s*r1


def radius_otegi(M_p):
    return np.exp(1/3.45*np.log(M_p/0.9))


def radius_zeng(M_p, CMF=None):
    # applicable to M_E <= 8 and CMF <= 0.4
#     print('using Zeng radius model')
    return (1.07 - 0.21*CMF)*(M_p/p.M_E)**(1/3.7)


def grav(M, R):
    """Calculate acceleration due to gravity on a point mass in m s^-2"""
    return 6.674e-11*M/R**2