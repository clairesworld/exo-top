from . import parameters
import numpy as np
from useful_and_bespoke import age_index
import sh_things as sh
from model_1D.topography import dimensionalise


def simple_vol_scaling(pl, verbose=False, **kwargs):
    # use ASPECT fit to h_peak and scale volume to entire surface area
    try:
        L = pl.L*8



def max_ocean(pl, n_stats=10, at_age=None, name_rms='dyn_top_aspect_prime', phi0=None, plot=False, verbose=False, **kwargs):

    h_rms1 = eval('pl.' + name_rms)
    if verbose:
        h_rms1_dim = eval('pl.' + 'dyn_top_rms')
        print('nondimensional h rms',  h_rms1[-1])
        print('dimensional h rms', h_rms1_dim[-1])

    l = np.arange(len(phi0))
    k = sh.l_to_k(l, R=2)  # original model spectrum uses R = 2d = 2
    # phi0_dim = phi0 * pl.d_m[-1] ** 3 * pl.dT_m[-1] ** 2 * pl.alpha_m ** 2  # dimensionalise power of model spec
    # k_dim = k * pl.d_m[-1] ** -1  # dimensionalise wavenumber of model spec
    h_rms0 = sh.parseval_rms(phi0, k)
    # h_rms0_dim = sh.parseval_rms(phi0_dim, k_dim)
    # print('mass', pl.M_p/parameters.M_E, 'M_E, dimensional rms of phi0', h_rms0_dim)

    # get time index nearest to desired snap given in Gyr
    if at_age is not None:
        it = age_index(pl.t, at_age, parameters.sec2Gyr)
        it = range(it, it + 1)
    else:
        it = range(len(pl.t))
        pl.max_ocean = np.zeros_like(it, dtype=np.float64)
    for ii in it:
        h_rms = h_rms1[ii]
        h_ratio = h_rms / h_rms0
        nn = 0
        vols = []
        # print('forcing h_ratio = 1')
        while nn < n_stats:
            clm = sh.random_harms_from_psd(phi0, l, R=2, h_ratio=h_ratio, plot=plot, verbose=verbose)
            grid = sh.coeffs_to_grid(clm, R=2, plot_grid=False, plot_spectrum=False, verbose=verbose,
                                     d=1, alpha_m=1, dT=1, scale_to_1D=False)
            grid_dim = dimensionalise(grid, pl, i=ii)
            if verbose:
                print('RMS of map dimensionalised', np.sqrt(np.mean(grid_dim**2)), 'm')

            vol = sh.integrate_to_peak(grid_dim, R=pl.R_p, fudge_to_rms=pl.dyn_top_rms[ii], verbose=verbose)
            vols.append(vol)
            nn = nn + 1
        basin_capacity = np.mean(vols)

        # phi = sh.scale_spectrum(h_rms=h_rms1[ii], **kwargs)  # old verj
        # h_peak_grid = sh.hpeak_from_spectrum(phi, n=n_stats, **kwargs)
        # basin_capacity = sh.vol_from_peak(r0=pl.R_p, h_peak=h_peak_grid, **kwargs)

        if at_age is None:
            pl.max_ocean[ii] = basin_capacity
        else:
            pl.max_ocean = basin_capacity

        # print('t', ii, 'R_p:', pl.R_p * 1e-3, 'km, h_rms', h_rms1[ii]*1e-3, 'km, h_peak', h_spectrum_max * 1e-3, 'km, ocn vol:', basin_capacity / parameters.TO, 'TO')
    # print('final h_ratio', h_ratio)

    if plot:
        from matplotlib.pyplot import show as pltshow
        pltshow()
    return pl
