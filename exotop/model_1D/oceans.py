from . import parameters
import numpy as np
from useful_and_bespoke import age_index
import sh_things as sh


def max_ocean(pl, n_stats=10, age=None, name_rms='dyn_top_rms', phi0=None, **kwargs):
    # get time index nearest to desired snap given in Gyr
    h_rms1 = eval('pl.' + name_rms)
    l = np.arange(len(phi0))
    k = sh.l_to_k(l, R=2)  # original model spectrum uses R = 2d = 2
    phi0_dim = phi0 * pl.d_m[-1]**3 * pl.dT_m[-1]**2 * pl.alpha_m**2  # dimensionalise power of model spec
    k_dim = k * pl.d_m[-1]**-1  # dimensionalise wavenumber of model spec
    h_rms0 = sh.parseval_rms(phi0_dim, k_dim)

    if age is not None:
        it = age_index(pl.t, age, parameters.sec2Gyr)
        it = range(it, it + 1)
    else:
        it = range(len(pl.t))
    pl.max_ocean = np.zeros_like(it, dtype=np.float64)
    for ii in it:
        h_rms = h_rms1[ii]
        h_ratio = h_rms/h_rms0
        nn = 0
        peaks = []
        while nn < n_stats:
            clm = sh.random_harms_from_psd(phi0, l, R=pl.R_p, h_ratio=h_ratio, plot=False, verbose=False)
            h_rms_grid, h_peak_grid = sh.coeffs_to_grid(clm, R=pl.R_p, plot_grid=False, plot_spectrum=False, verbose=False,
                                                        d=pl.d_m[ii], alpha_m=pl.alpha_m, dT=pl.dT_m[ii])
            peaks.append(h_peak_grid)
        h_peak_grid = np.mean(peaks)

        # phi = sh.scale_spectrum(h_rms=h_rms1[ii], **kwargs)
        # h_peak_grid = sh.hpeak_from_spectrum(phi, n=n_stats, **kwargs)

        basin_capacity = sh.vol_from_peak(r0=pl.R_p, h_peak=h_peak_grid, **kwargs)
        pl.max_ocean[ii] = basin_capacity
        # print('t', ii, 'R_p:', pl.R_p * 1e-3, 'km, h_rms', h_rms1[ii]*1e-3, 'km, h_peak', h_spectrum_max * 1e-3, 'km, ocn vol:', basin_capacity / parameters.TO, 'TO')
    return pl
