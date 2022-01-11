from . import parameters
import numpy as np
from useful_and_bespoke import age_index, dark_background
import sh_things as sh
from model_1D.topography import dimensionalise, dyn_topo_peak_prime_aspect


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
            pl.simple_ocean[ii] = h_peak * (4 * np.pi * R ** 2)
        else:
            if verbose:
                print('h peak from scaling', h_peak, 'm')
                print('h rms', pl.dyn_top_rms[ii], 'm')
            pl.simple_ocean = h_peak * (4 * np.pi * R ** 2)

    return pl


def max_ocean_fast(pl, peak_ratio=3.5, at_age=None, name_rms='dyn_top_aspect_prime', rho_w=1000, verbose=False, **kwargs):
    water_load_ratio = pl.rho_m / (pl.rho_m - rho_w)
    if verbose:
        print('including water loading |', 'peak ratio', peak_ratio)

    h_rms_prime = eval('pl.' + name_rms)

    if at_age is not None:
        it = age_index(pl.t, at_age, parameters.sec2Gyr)  # get time index nearest to desired snap given in Gyr
        it = range(it, it + 1)
    else:  # all times
        it = range(len(pl.t))
        pl.max_ocean = np.zeros_like(it, dtype=np.float64)
        pl.h_peak_spectral = np.zeros_like(it, dtype=np.float64)
    for ii in it:
        h_rms_ii = dimensionalise(h_rms_prime, pl, i=ii) * water_load_ratio
        h_peak_ii = h_rms_ii * peak_ratio
        V_cap_ii = 4 * np.pi * pl.R_p**2 * h_peak_ii
        # print('V_cap_ii', V_cap_ii, 'R_p', pl.R_p, 'h_peak_ii', h_peak_ii)

        if at_age is None:
            pl.max_ocean[ii] = V_cap_ii
            pl.h_peak_spectral[ii] = h_peak_ii
        else:
            pl.max_ocean = V_cap_ii
            pl.h_peak_spectral = h_peak_ii
    return pl


def max_ocean(pl, n_stats=10, at_age=None, name_rms='dyn_top_aspect_prime', phi0=None, plot=False, verbose=False,
              spectrum_fname='base_spectrum_l1.pkl', spectrum_fpath='/home/claire/Works/exo-top/exotop/top_spectra/',
              rho_w=1000, **kwargs):
    water_load_ratio = pl.rho_m / (pl.rho_m - rho_w)
    if verbose:
        print('including water loading')

    if phi0 is None:
        degree, phi0 = sh.load_model_spectrum_pkl(fname=spectrum_fname, path=spectrum_fpath)

    h_rms1 = eval('pl.' + name_rms)
    if verbose:
        h_rms1_dim = eval('pl.' + 'dyn_top_rms')
        print('nondimensional h rms', h_rms1[-1])
        print('dimensional h rms', h_rms1_dim[-1])

    l = np.arange(len(phi0))
    k = sh.l_to_k(l, R=2)  # original model spectrum uses R_b = 2d = 2
    h_rms0 = sh.parseval_rms(phi0, k)
    if np.isnan(h_rms0):
        print('phi0', phi0)
        print('k', k)
        raise Exception('rms of model spectrum is nan!')

    if at_age is not None:
        it = age_index(pl.t, at_age, parameters.sec2Gyr)  # get time index nearest to desired snap given in Gyr
        it = range(it, it + 1)
    else:  # all times
        it = range(len(pl.t))
        pl.max_ocean = np.zeros_like(it, dtype=np.float64)
        pl.h_peak_spectral = np.zeros_like(it, dtype=np.float64)
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
            # shgrid = sh.coeffs_to_grid(clm, R_b=2, plot_grid=False, plot_spectrum=False, verbose=verbose,)
            # lmax = shgrid.lmax
            # data = shgrid.data
            # lats = shgrid.lats()
            # lons = shgrid.lons()
            # vol = sh.integrate_to_peak(grid_dim, lats, lons, R_b=pl.R_p, lmax=shgrid.lmax, verbose=verbose)
            data = clm.expand(grid='GLQ', extend=False).to_array()
            lmax = clm.lmax
            rms_nondim = np.sqrt(np.mean(data ** 2))
            # if verbose:
            #     print('RMS of map nondimensional', rms_nondim)
            grid_dim = dimensionalise(data, pl, i=ii) * water_load_ratio
            vol = sh.integrate_to_peak_GLQ(grid_dim, R=pl.R_p, lmax=lmax, verbose=verbose)
            if np.isnan(vol):
                print('phi0', phi0, 'h_ratio', h_ratio)
                raise Exception('ocean vol nan!')
            vols.append(vol)
            rms_dim = np.sqrt(np.mean(grid_dim ** 2))
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
            pl.h_peak_spectral[ii] = np.mean(peaks)
        else:
            pl.max_ocean = basin_capacity
            pl.h_peak_spectral = np.mean(peaks)
            if verbose:
                print('\n')
                print('      scaled h_rms', name_rms, '=', dimensionalise(h_rms1[-1], pl=pl, i=-1), 'm without water loading')
                print('      water load ratio', water_load_ratio)
                print("%.2f" % (pl.M_p / parameters.M_E), 'M_E |', "%.2f" % (pl.R_p * 1e-3), 'km radius | ocean vol:',
                      "%.2f" % (basin_capacity / 1.4e18), 'EO | V_shell =',
                      "%.2f" % (4 / 3 * np.pi * ((pl.R_p + np.mean(peaks)) ** 3 - pl.R_p ** 3) / 1.4e18), 'EO | h_peak =',
                      "%.2f" % (np.mean(peaks) * 1e-3), 'km | h_rms =', "%.2f" % (np.mean(rms_dims) * 1e-3), 'km')
                print('h_rms_prime =', "%.5f" % np.mean(rms_nondims), '| b =', "%.2f" % pl.b[-1], '| log Ra_i =',
                      "%.2f" % np.log10(pl.Ra_i[-1]), '| d =', "%.2f" % (pl.d * 1e-3), '| dT =', "%.2f" % pl.delta_T[-1],
                      '| alpha_m =', "%.2f" % pl.alpha_m)

        # print('t', ii, 'R_p:', pl.R_p * 1e-3, 'km, h_rms', h_rms1[ii]*1e-3, 'km, h_peak', h_spectrum_max * 1e-3, 'km, ocn vol:', basin_capacity / parameters.TO, 'TO')
    # print('final h_ratio', h_ratio)

    if plot:
        from matplotlib.pyplot import show as pltshow
        pltshow()
    return pl


def min_topo(x_h2o, R_p, M_p, n_stats=50, rms_1=1000, tol=0.5, phi0=None, rho_m=3500, rho_w=1000,
             spectrum_fname='base_spectrum_l1.pkl',
             spectrum_fpath='/home/claire/Works/exo-top/exotop/figs_scratch/', verbose=False, **kwargs):
    """ predict min rms topography for land given sfc water mass frac and radius"""
    vol_EO = 1.4e21 / rho_w
    water_load_ratio = rho_m / (rho_m - rho_w)

    if phi0 is None:
        degree, phi0 = sh.load_model_spectrum_pkl(fname=spectrum_fname, path=spectrum_fpath)

    l = np.arange(len(phi0))
    k = sh.l_to_k(l, R=2)  # original model spectrum uses R_b = 2d = 2
    h_rms0 = sh.parseval_rms(phi0, k)
    M_w = M_p * parameters.M_E * x_h2o  # mass of sfc water in kg
    vol_w = M_w / rho_w  # corresponding volume

    print('vol_w', vol_w / vol_EO, 'EO')

    h_rms = rms_1
    flag = True
    while flag:
        # get vol from rms
        vols = []
        n = 0
        while n < n_stats:
            h_ratio = h_rms / h_rms0
            clm = sh.random_harms_from_psd(phi0, l, R=2, h_ratio=h_ratio, plot=False, verbose=verbose)
            grid = clm.expand(grid='GLQ', extend=False).to_array() * water_load_ratio
            lmax = clm.lmax
            vol = sh.integrate_to_peak_GLQ(grid, R=R_p, lmax=lmax, verbose=verbose)
            vols.append(vol)
            n = n + 1
        vol = np.mean(vols)
        diff = vol_w - vol
        if abs(diff) < tol * vol_EO:
            flag = False
            print('h_rms', h_rms, 'm | vol', vol / vol_EO, 'EO')
        else:
            h_rms = h_rms * (vol_w / vol)
    print('ANS:', 'h_rms', h_rms, 'm | vol', vol / vol_EO, 'EO')
    return h_rms

# def no_h_scaling(pl):


def plot_map(pl, at_age=4.5, phi0=None, name_rms='dyn_top_aspect_prime', spectrum_fname='base_spectrum_l1.pkl',
             spectrum_fpath='/home/claire/Works/exo-top/exotop/top_spectra/', rho_w=1000, verbose=False,
             fig_path='/home/claire/Works/exo-top/exotop/figs_scratch/', cbar=True, clabel='Dynamic topography (m)',
             cmap='gist_earth', labelsize=16, ticksize=14, clabelpad=20, fig_fmt='.pdf', dark=False, save=True,
             axlabels=True, bg_colour='k', fg_colour='xkcd:off white', fname='topo_grid', texfont=True):
    import matplotlib.pyplot as plt
    from matplotlib import rc
    from matplotlib.pyplot import rcParams
    rc('text', usetex=True)  # turn off for running over ssh
    if texfont:
        rcParams['font.family'] = 'serif'
        rcParams['font.serif'] = 'CMU Serif'

    water_load_ratio = pl.rho_m / (pl.rho_m - rho_w)
    print('including water loading')

    if phi0 is None:
        degree, phi0 = sh.load_model_spectrum_pkl(fname=spectrum_fname, path=spectrum_fpath)

    h_rms1 = eval('pl.' + name_rms)
    if verbose:
        h_rms1_dim = eval('pl.' + 'dyn_top_rms')
        print('nondimensional h rms', h_rms1[-1])
        print('dimensional h rms', h_rms1_dim[-1])

    l = np.arange(len(phi0))
    k = sh.l_to_k(l, R=2)  # original model spectrum uses R_b = 2d = 2
    h_rms0 = sh.parseval_rms(phi0, k)
    if np.isnan(h_rms0):
        print('phi0', phi0)
        print('k', k)
        raise Exception('rms of model spectrum is nan!')

    ii = age_index(pl.t, at_age, parameters.sec2Gyr)  # get time index nearest to desired snap given in Gyr

    h_rms = h_rms1[ii]
    h_ratio = h_rms / h_rms0
    clm = sh.random_harms_from_psd(phi0, l, R=2, h_ratio=h_ratio, plot=False, verbose=verbose)
    data = sh.coeffs_to_grid(clm, R=2, plot_grid=False, plot_spectrum=False, verbose=verbose)
    lmax = clm.lmax
    print('lmax', lmax)
    lats = data.lats()
    lons = data.lons()
    grid_dim = dimensionalise(data.data, pl, i=ii) * water_load_ratio
    fig, ax = plt.subplots(1, 1)
    mappable = ax.imshow(grid_dim, extent=(0, 360, -90, 90), cmap=cmap, interpolation='gaussian')
    ax.set(yticks=np.arange(-90, 120, 30), xticks=np.arange(0, 390, 30))
    ax.tick_params(axis='both', labelsize=0)
    if axlabels:
        ax.set_xlabel('Latitude', fontsize=labelsize)
        ax.set_ylabel('Longitude', fontsize=labelsize)
    if cbar:
        cb = plt.colorbar(mappable, orientation='horizontal', location='top', fraction=0.07)
        cb.ax.tick_params(axis='x', labelsize=ticksize)
        if dark:
            textc = fg_colour
            # set colorbar tick color
            cb.ax.xaxis.set_tick_params(color=textc)

            # set colorbar edgecolor
            cb.outline.set_edgecolor(textc)

            # set colorbar ticklabels
            plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color=textc)
        else:
            textc = 'k'
        cb.set_label(label=clabel, fontsize=labelsize, labelpad=clabelpad, c=textc)
    if dark:
        fig, ax = dark_background(fig, ax, bgc=bg_colour,  fgc=fg_colour)
    if save:
        plt.savefig(fig_path + fname + fig_fmt, bbox_inches='tight',
                    facecolor=fig.get_facecolor())
    return fig, ax




