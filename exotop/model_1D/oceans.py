from . import parameters
import numpy as np
import sys

sys.path.append("..")
import exotop.asharms as harm


def max_ocean(pl, n_stats=10, age=None, name_rms='dyn_top_rms', **kwargs):
    # get time index nearest to desired snap given in Gyr
    h_rms1 = eval('pl.' + name_rms)
    if age is not None:
        it = min(enumerate(pl.t), key=lambda x: abs(age - x[1] * parameters.sec2Gyr))[0]
        it = range(it, it + 1)
    else:
        it = range(len(pl.t))
    pl.max_ocean = np.zeros_like(it, dtype=np.float64)
    for ii in it:
        phi = harm.scale_spectrum(h_rms=h_rms1[ii], **kwargs)
        h_spectrum_max = harm.hpeak_from_spectrum(phi, n=n_stats, **kwargs)
        basin_capacity = harm.vol_from_peak(r0=pl.R_p, h_peak=h_spectrum_max, **kwargs)
        pl.max_ocean[ii] = basin_capacity
        # print('t', ii, 'R_p:', pl.R_p * 1e-3, 'km, h_rms', h_rms1[ii]*1e-3, 'km, h_peak', h_spectrum_max * 1e-3, 'km, ocn vol:', basin_capacity / parameters.TO, 'TO')
    return pl
