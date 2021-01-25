import numpy as np
import matplotlib.pyplot as plt
import pyshtools

def RMS_from_l_powerspectrum(clm_topo, lmax=None):
    # Calculate RMS topography from power spectrum
    spec = clm_topo.spectrum(unit='per_l', lmax=lmax)
    RMS_topo = np.sqrt(np.sum(spec))
    print("RMS_topo", RMS_topo, "m")
    return RMS_topo


def read_coeffs(filename='VenusTopo719.shape', lmax=719, path='benchmarks/'):
    # Read in spherical harmonic coefficients of topography
    lmax_read = lmax
    cilm_topo, lmax_out = pyshtools.shio.shread(path+filename, lmax=lmax_read)
    clm_topo = pyshtools.SHCoeffs.from_array(cilm_topo)

    # Radius of Venus in m
    a0 = pyshtools.constants.Venus.r.value

    # Get rid of degree 0 corresponding to mean radius
    clm_topo.coeffs[:,0] = 0.0
    return clm_topo, cilm_topo


def plot_Venus(filename='VenusTopo719.shape', lmax=719, path='benchmarks/'):
    clm_topo, cilm_topo = read_coeffs(filename=filename, lmax=lmax, path=path)

    # Some different ways of plotting power spectrum
    fig, ax = clm_topo.plot_spectrum(unit='per_l') # per degree
    plt.ylabel("Power (m$^2$)")
    fig, ax = clm_topo.plot_spectrum(unit='per_lm') # per coefficient
    plt.ylabel("Power (m$^2$)")

    # Calculate RMS topography from power spectrum
    spec = clm_topo.spectrum(unit='per_l')
    RMS_topo = np.sqrt(np.sum(spec))
    print("RMS_topo", RMS_topo, "m")

    # Make a plot like Rappaport et al 1999
    rms_l_spec = (1.0/r0) * np.sqrt(clm_topo.spectrum(unit='per_lm'))
    plt.figure()
    plt.loglog(np.arange(0,lmax_read+1), rms_l_spec)
    plt.xlim([1.0, 1e3])
    plt.ylim([1e-8, 1e-4])
    plt.xlabel("Degree")
    plt.ylabel("RMS (dimensionless)")
    plt.title("Like Figure 7 of Rappaport et al 1999")
    plt.show()
   
