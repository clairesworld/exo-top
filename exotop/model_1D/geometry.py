""" Set up planet physical structure """

import numpy as np
from . import parameters as p
from .astroenvironment import radius_zeng


###### GEOMETRY ######

def SA(M=None, rho=None, R=None, V=None):
    if R is None:
        if (V is None) and (M is not None) and (rho is not None):
            V = vol(M, rho)
        else:
            print('Missing mass or density to calculate volume')
        R = radius_by_vol(V)
    return 4*np.pi*R**2


def radius_by_vol(V):
    return (3*V/(4*np.pi))**(1/3)


def vol(M, rho):
    try:
        return M/rho
    except TypeError:
        return np.array(M)/np.array(rho)
    

def core_density(M, CMF):  # Zeng & Jacobsen eq 10
    G = 6.674e-11
    R = radius_zeng(M, CMF)*p.R_E
    R_c = R*CMF**0.5
    g_s = G*M/R**2
    return 3*g_s / (4*np.pi*G*R_c)


def mantle_density(r, M, CMF): # Zeng & Jacobsen eq 10
    G = 6.674e-11
    R = radius_zeng(M, CMF)*p.R_E
    g_s = G*M/R**2
    return g_s/(2*np.pi*G*r)


def mantle_mass(r, M, CMF):
    G = 6.674e-11
    R = radius_zeng(M, CMF)*p.R_E
    return (r/R)**2 * M



