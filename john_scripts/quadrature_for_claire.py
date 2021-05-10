import pyshtools as pysh
import numpy as np

# Generate a random spherical harmonic model for testing
l = 100 # max spherical harmonic degree
degrees = np.arange(l+1, dtype=float)
degrees[0] = np.inf
power = degrees**(-2)
clm = pysh.SHCoeffs.from_random(power, seed=12345)

# Expand function on to grid points. Grid has l+1 points in latitude, 2l + 1 points in longitude
f = clm.expand(grid='GLQ').to_array()

# calculate Gauss-Legendre weights for quadrature rule
nodes, lat_weights = pysh.expand.SHGLQ(l)  # gives weights for latitude
lon_weights = 2.0*np.pi/(2*l +1) # weights for longitude are uniform
w = lon_weights * np.tile(lat_weights, (2*l +1,1)).T  # 2d grid of lat-lon weights for quadrature

# Some example integrals
integral1 = np.sum(w)   # Integral of 1 -- should give 4 pi, surface area of sphere
integral2 = np.sum(w*f) # Integral of function -- for above example, mean is zero 
integral3 = (1.0/(4.0*np.pi)) * np.sum(w*f*f)  # Mean square of function, check of Parseval's theorem

print("integral1", integral1, "expected", 4.0*np.pi)
print("integral2", integral2, "expected", 0.0)
print("integral3", integral3, "expected", np.sum(clm.spectrum()))
