import sh_things as sh
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from postaspect.plt_aspect import plot_save
from postaspect.setup_postprocessing import data_path_bullard, fig_path_bullard
from useful_and_bespoke import dark_background

data_path = '/home/claire/Works/aspect/runs/model-output/'
fig_path = '/home/claire/Works/exo-top/exotop/figs_scratch/'
benchmark_path = '/home/cmg76/Works/exo-top/benchmarks/lees_topo_grids/'


def cmp_dimensions(benchmark_path='/home/cmg76/Works/exo-top/benchmarks/lees_topo_grids/',
                   data_path='/raid1/cmg76/aspect/model-output/', save=True, fig_path=fig_path_bullard, **kwargs):
    """ load 3D spectrum from digitised plot """
    df3 = pd.read_csv(benchmark_path + 'psd_freefree.csv', index_col=False, comment='#', header=0)
    psd3 = df3.y.to_numpy()  # km^2 km^2
    wl3 = df3.x.to_numpy()  # km
    k3 = 2 * np.pi / wl3

    d = 600  # km
    alpha = 4e-5
    dT = 442  # potential temperature difference K

    """ load 2D spectrum from aspect data """
    psd2, k2 = sh.load_model_spectrum_pkl(path='', data_path=data_path_bullard, case='Lees-Ra1e6-2D')
    psd2 = psd2 * d ** 3 * alpha ** 2 * dT ** 2
    k2 = k2 * d ** -1

    # convert from km^3 to km^4
    psd2_iso = 1 / k2 * psd2  # Jacobs eqn 5 but pi changed to 1 in numerator says JFR

    """ plot """
    fig, ax = plt.subplots()
    ax.loglog(k2, psd2_iso, 'o-', lw=1, alpha=0.5, c='xkcd:slate', label='2D')
    ax.loglog(k3, psd3, '^-', lw=1, alpha=0.5, c='xkcd:forest green', label='3D')
    ax.set_xlabel(r'Wavenumber (km$^{-1}$)')
    ax.set_ylabel(r'PSD (km$^{2}$ km$^{2}$)')
    ax.legend()

    if save:
        plot_save(fig, fname='2D_benchmark', fig_path=fig_path)


def Venus_correction(baseline_fname='base_spectrum.pkl', fig_path=fig_path_bullard, R_base=2, lmin=1,
                     save=True, plot=True, units='km4', scale_to='Venus', labelsize=16, alpha=0.5, **kwargs):
    R_Venus = sh.get_pysh_constants('Venus', 'r') * 1e-3  # in km

    # get 2D PSDs
    l, phi = sh.load_model_spectrum_pkl(fname=baseline_fname, path=fig_path, **kwargs)
    k = sh.l_to_k(l, R_Venus)
    phi_iso = 1 / k * phi
    lmax = np.max(l)

    lV, phiV = sh.get_psd_Venus(unit='per_lm', to_1D=False, lmax=lmax, to_km=True,
                                verbose=False)  # power per lm in km^2 km^2
    kV = sh.l_to_k(lV, R_Venus)

    # remove 0 degree
    l, phi_iso, phi, k = l[lmin:], phi_iso[lmin:], phi[lmin:], k[lmin:]
    lV, phiV, kV = lV[lmin:], phiV[lmin:], kV[lmin:]

    # scale to Venus RMS at power
    rms_V = sh.parseval_rms_2D(phiV, kV)
    rms_base = sh.parseval_rms_2D(phi_iso, k)
    print('rms Venus', rms_V, 'km')
    print('rms base', rms_base, 'nondim')

    S = np.array(phi_iso) * (2 * l + 1) / (4 * np.pi * R_Venus ** 2)  # factor of 4piR^2 from Lees eq A7
    SV = np.array(phiV) * (2 * lV + 1) / (4 * np.pi * R_Venus ** 2)

    if scale_to == 'base':
        rms0 = rms_base
    elif scale_to == 'Venus':
        rms0 = rms_V

    S_sc = S * (rms0 / rms_base) ** 2
    phi_sc = 4 * np.pi * R_Venus ** 2 * S_sc / (2 * l + 1)
    SV_sc = SV * (rms0 / rms_V) ** 2
    phiV_sc = 4 * np.pi * R_Venus ** 2 * SV_sc / (2 * lV + 1)
    print('rms base new', sh.parseval_rms_2D(phi_sc, k))
    print('rms Venus new', sh.parseval_rms_2D(phiV_sc, kV))

    if units == 'km3':
        # convert back to 1D
        phiV_sc = phiV_sc * kV
        phi_sc = phi_sc * k

    if plot:
        fig, ax = plt.subplots()
        ax.loglog(lV, phiV_sc, 'o-', lw=1, alpha=alpha, c='xkcd:slate', label='Venus (Wieczorek 2015)')
        ax.loglog(l, phi_sc, '^-', lw=1, alpha=alpha, c='xkcd:lime green', label='2D ASPECT @ Venus rms')
        ax.set_xlabel(r'Degree', fontsize=labelsize)
        if units == 'km4':
            ylabel = r'2D PSD (km$^{2}$ km$^{2}$)'
        elif units == 'km3':
            ylabel = r'1D PSD (km$^{2}$ km)'
        ax.set_ylabel(ylabel, fontsize=labelsize)
        ax.legend(frameon=False, fontsize=labelsize-2)
        if save:
            plot_save(fig, fname='Venus_correction', fig_path=fig_path)

    if lmin > 0:
        lV = np.arange(0, lmax + 1)
        phiV_sc = np.insert(np.array(phiV_sc), 0, [0.0] * lmin)  # no power below lmin

    if plot:
        return lV, phiV_sc, fig, ax
    return lV, phiV_sc


def make_Venus_reference(newfname='Venus_spectrum.pkl', baseline_fname='base_spectrum.pkl', fig_path='',
                         lmin=1, plot=True):
    import pickle as pkl

    # want to 'nondimensionalise' observed Venus spectrum by scaling it such that it has same rms as baseline and
    # corresponds to R=2
    l, Sl = Venus_correction(baseline_fname=baseline_fname, fig_path=fig_path, R_base=2,
                             save=False, plot=plot, lmin=lmin, units='km3', scale_to='base')

    pkl.dump((l, Sl), open(fig_path + newfname, "wb"))
    return l, Sl


lV, phiV_sc, fig, ax = Venus_correction(baseline_fname='base_spectrum_l1.pkl', fig_path=fig_path, R_base=2, lmin=1,
                     save=True, plot=True, units='km4', scale_to='Venus', alpha=0.9, labelsize=16)

fig, ax = dark_background(fig, ax)
# make_Venus_reference(newfname='Venus_spectrum_l1.pkl', baseline_fname='base_spectrum_l1.pkl',
#                      fig_path=fig_path)
plt.show()
