import pyshtools
# import cartopy.crs as ccrs
import numpy as np
import pandas as pd


def hpeak_from_spectrum(power, norm='4pi', lmax=40, n=10, **kwargs):
    ii = 0
    h_peak = []
    while ii < n:
        coeffs_global = pyshtools.SHCoeffs.from_random(power, normalization=norm)
        topo = coeffs_global.expand(lmax=lmax)
        h_peak.append(np.max(topo.data))  # in m
        ii += 1
    return np.mean(h_peak)  # average of a bunch of random spectra consistent with given


def load_spectrum(fpath='', fname='', **kwargs):
    df = pd.read_csv(fpath + fname, header=None, names=['l', 'S_l'], index_col=False)
    degrees = df['l']
    power = df['S_l']
    # todo : find nearest l for imperfect digitization
    return power


def scale_spectrum(h_rms, h_rms0=None, phi0=None, pl=None, pl0=None, h_func=None, age=4.5, **kwargs):
    # scale power spectrum phi0 such that it has new rms equal to h_rms

    # if h_func is None:
        # here you are trying to scale dynamic topography spectrum by ratio of h_rms derived from known 1D params
        # get time index nearest to desired snap given in Gyr
        # it = min(enumerate(pl_baseline.t), key=lambda x: abs(age - x[1] * parameters.sec2Gyr))[0]
        # ratio = (pl.d_m * pl.T_m ** 2 * pl.Ra_i ** (-1 / 3)) / (
        #         pl0.d_m * pl0.T_m ** 2 * pl0.Ra_i ** (-1 / 3))  # = phi / phi0
        # ratio = ratio[it]
    if h_rms0 is None:
        h_rms0 = powerspectrum_RMS(power_lm=phi0)  # I hope this is a power spectrum, but this probably needs to be multipled by R or something
    ratio = h_rms/h_rms0
    sqrtphi = np.sqrt(phi0) * ratio
    return sqrtphi ** 2


def vol_from_peak(r0, h_peak, vol_ref=1, **kwargs):
    # getting volume in between trough and peak maxima - as troughs fill in peaks, net volume is 0 and total volume is
    # the spherical annulus between r0 and h_peak
    # R_p = 6051.88e3
    EO = 1.35e9 * 1000**3  # e.g. reference volume based on R_E
    vol = 4 / 3 * np.pi * ((r0 + h_peak) ** 3 - r0 ** 3)
    # print('R_p:', r0*1e-3, 'km, h_peak', h_peak*1e-3, 'km, vol:', vol/EO, 'TO')
    return vol / vol_ref


def vol_from_spectrum(phi0=None, fname='', fpath='', r0=1, n_stats=10, **kwargs):
    # volume between min trough and max peak
    d_ocean_Earth = 3.7e3
    if phi0 is None:
        phi0 = load_spectrum(fpath=fpath, fname=fname)
    h_peak0 = hpeak_from_spectrum(phi0, n=n_stats)
    V0 = vol_from_peak(r0, h_peak0)
    return V0


def powerspectrum_RMS(path=None, power_lm=None, amplitude=False, n=None): # try to calcuate RMS from digitized power spectrum
    if path is not None:
        df = pd.read_csv(path, header=None, names=['degree', 'value'], index_col=False)
        ls = np.array(df['degree'])
        S = np.array(df['value'])
    elif power_lm is not None:
        ls = np.arange(1, len(power_lm)-1)
        S = power_lm

    if n is None:
        n = len(ls)

    RMS_l = []
    for ii, l in enumerate(ls[0:n]):
        Slm = S[ii]
        if amplitude:
            Slm = Slm**2
        RMS_l.append(np.sqrt(Slm/(2*l + 1)))
    return sum(RMS_l)
