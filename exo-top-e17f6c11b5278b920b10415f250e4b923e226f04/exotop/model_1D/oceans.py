from . import parameters
import numpy as np
from exotop.useful_and_bespoke import age_index
import exotop.sh_things as sh
from exotop.model_1D.topography import dimensionalise, dyn_topo_peak_prime_aspect


def simple_vol_scaling(pl, verbose=False, at_age=None, it=None, **kwargs):
    # use ASPECT fit to h_peak and scale volume to entire surface area
    h_peak_prime = dyn_topo_peak_prime_aspect(pl)
    R = pl.R_p

    # get time index nearest to desired snap given in Gyr
    if at_age is not None:
        i = age_index(pl.t, at_age, parameters.sec2Gyr)
        it = range(i, i + 1)
    else:
        it = range(len(pl.t))
        pl.simple_ocean = np.zeros_like(it, dtype=np.float64)

    for ii in it:
        h_peak = dimensionalise(h_peak_prime[ii], pl=pl, i=ii)
        if at_age is None:
            pl.simple_ocean[ii] = h_peak * (4 * np.pi * R**2)
        else:
            if verbose:
                print('h peak from scaling', h_peak, 'm')
                print('h rms', pl.dyn_top_rms[ii], 'm')
            pl.simple_ocean = h_peak * (4 * np.pi * R**2)

    return pl


def max_ocean(pl, n_stats=10, at_age=None, name_rms='dyn_top_aspect_prime', phi0=None, plot=False, verbose=False,
              spectrum_fname='base_spectrum_l1.pkl', spectrum_fpath='/home/claire/Works/exo-top/exotop/figs_scratch/',
              **kwargs):
    if phi0 is None:
        degree, phi0 = sh.load_model_spectrum_pkl(fname=spectrum_fname, path=spectrum_fpath)

    h_rms1 = eval('pl.' + name_rms)
    if verbose:
        h_rms1_dim = eval('pl.' + 'dyn_top_rms')
        print('nondimensional h rms',  h_rms1[-1])
        print('dimensional h rms', h_rms1_dim[-1])

    l = np.arange(len(phi0))
    k = sh.l_to_k(l, R=2)  # original model spectrum uses R = 2d = 2
    h_rms0 = sh.parseval_rms(phi0, k)
    if np.isnan(h_rms0):
        print('phi0', phi0)
        print('k', k)
        raise Exception('rms of model spectrum is nan!')

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
        peaks = []
        rms_dims = []
        rms_nondims = []
        while nn < n_stats:
            clm = sh.random_harms_from_psd(phi0, l, R=2, h_ratio=h_ratio, plot=plot, verbose=verbose)
            # shgrid = sh.coeffs_to_grid(clm, R=2, plot_grid=False, plot_spectrum=False, verbose=verbose,)
            # lmax = shgrid.lmax
            # data = shgrid.data
            # lats = shgrid.lats()
            # lons = shgrid.lons()
            # vol = sh.integrate_to_peak(grid_dim, lats, lons, R=pl.R_p, lmax=shgrid.lmax, verbose=verbose)
            data = clm.expand(grid='GLQ', extend=False).to_array()
            lmax = clm.lmax
            rms_nondim = np.sqrt(np.mean(data**2))
            # if verbose:
            #     print('RMS of map nondimensional', rms_nondim)
            grid_dim = dimensionalise(data, pl, i=ii)
            vol = sh.integrate_to_peak_GLQ(grid_dim, R=pl.R_p, lmax=lmax, verbose=verbose)
            if np.isnan(vol):
                print('phi0', phi0, 'h_ratio', h_ratio)
                raise Exception('ocean vol nan!')
            rms_dim = np.sqrt(np.mean(grid_dim ** 2))
            vols.append(vol)
            peaks.append(np.max(grid_dim))
            rms_nondims.append(rms_nondim)
            rms_dims.append(rms_dim)
            nn = nn + 1
            if verbose:
                # print('RMS of map dimensionalised', rms_dim, 'm')
                if nn > 0:
                    verbose = False  # don't repeat output
        basin_capacity = np.mean(vols)

        # phi = sh.scale_spectrum(h_rms=h_rms1[ii], **kwargs)  # old verj
        # h_peak_grid = sh.hpeak_from_spectrum(phi, n=n_stats, **kwargs)
        # basin_capacity = sh.vol_from_peak(r0=pl.R_p, h_peak=h_peak_grid, **kwargs)

        if at_age is None:
            pl.max_ocean[ii] = basin_capacity
        else:
            pl.max_ocean = basin_capacity
            print(pl.M_p/parameters.M_E, 'M_E |', basin_capacity/1.4e18, 'EO | V_shell =', 4/3*np.pi*((pl.R_p + np.mean(peaks))**3 - pl.R_p**3)/1.4e18, 'EO | h_peak =', np.mean(peaks)*1e-3, 'km | h_rms =', np.mean(rms_dims)*1e-3, 'km =', np.mean(rms_nondims))

        # print('t', ii, 'R_p:', pl.R_p * 1e-3, 'km, h_rms', h_rms1[ii]*1e-3, 'km, h_peak', h_spectrum_max * 1e-3, 'km, ocn vol:', basin_capacity / parameters.TO, 'TO')
    # print('final h_ratio', h_ratio)

    if plot:
        from matplotlib.pyplot import show as pltshow
        pltshow()
    return pl

def min_topo(x_h2o, R_p, M_p, n_stats=50, rms_1=1000, tol=0.5, phi0=None, rho_w=1000, spectrum_fname='base_spectrum_l1.pkl',
             spectrum_fpath='/home/claire/Works/exo-top/exotop/figs_scratch/', verbose=False, **kwargs):
    """ predict min rms topography for land given sfc water mass frac and radius"""
    vol_EO = 1.4e21/rho_w

    if phi0 is None:
        degree, phi0 = sh.load_model_spectrum_pkl(fname=spectrum_fname, path=spectrum_fpath)

    l = np.arange(len(phi0))
    k = sh.l_to_k(l, R=2)  # original model spectrum uses R = 2d = 2
    h_rms0 = sh.parseval_rms(phi0, k)
    M_w = M_p * parameters.M_E * x_h2o  # mass of sfc water in kg
    vol_w = M_w / rho_w  # corresponding volume

    print('vol_w', vol_w/vol_EO, 'EO')

    h_rms = rms_1
    flag = True
    while flag:
        # get vol from rms
        vols = []
        n = 0
        while n < n_stats:
            h_ratio = h_rms / h_rms0
            clm = sh.random_harms_from_psd(phi0, l, R=2, h_ratio=h_ratio, plot=False, verbose=verbose)
            grid = clm.expand(grid='GLQ', extend=False).to_array()
            lmax = clm.lmax
            vol = sh.integrate_to_peak_GLQ(grid, R=R_p, lmax=lmax, verbose=verbose)
            vols.append(vol)
            n = n + 1
        vol = np.mean(vols)
        diff = vol_w - vol
        if abs(diff) < tol*vol_EO:
            flag = False
            print('h_rms', h_rms, 'm | vol', vol / vol_EO, 'EO')
        else:
            h_rms = h_rms * (vol_w/vol)
    print('ANS:', 'h_rms', h_rms, 'm | vol', vol / vol_EO, 'EO')
    return h_rms