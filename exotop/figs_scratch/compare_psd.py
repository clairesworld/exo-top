import sh_things as sh
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from postaspect.plt_aspect import plot_save
from postaspect.setup_postprocessing import data_path_bullard, fig_path_bullard

def cmp_dimensions(benchmark_path='/home/cmg76/Works/exo-top/benchmarks/lees_topo_grids/',
                   data_path='/raid1/cmg76/aspect/model-output/', save=True, fig_path=fig_path_bullard, **kwargs):
    """ load 3D spectrum from digitised plot """
    df3 = pd.read_csv(benchmark_path + 'psd_freefree.csv', index_col=False, comment='#', header=0)
    psd3 = df3.y.to_numpy()  # km^2 km^2
    wl3 = df3.x.to_numpy()  # km
    k3 = 2*np.pi/wl3

    d = 600  # km
    alpha = 4e-5
    dT = 442  # potential temperature difference K

    """ load 2D spectrum from aspect data """
    psd2, k2 = sh.load_model_spectrum_pkl(path='', data_path=data_path_bullard, case='Lees-Ra1e6-2D')
    psd2 = psd2 * d ** 3 * alpha ** 2 * dT ** 2
    k2 = k2 * d**-1

    # convert from km^3 to km^4
    psd2_iso = 1 / k2 * psd2  # Jacobs eqn 5 but pi changed to 1 in numerator says JFR

    """ plot """
    fig, ax = plt.subplots()
    ax.loglog(k2, psd2_iso,  '.-', lw=0.5, alpha=0.5, c='xkcd:slate', label='2D')
    ax.loglog(k3, psd3,  '^-', lw=0.5, alpha=0.5, c='xkcd:forest green', label='3D')
    ax.set_xlabel(r'Wavenumber (km$^{-1}$)')
    ax.set_ylabel(r'PSD (km$^{2}$ km$^{2}$)')
    ax.legend()

    if save:
        plot_save(fig, fname='2D_benchmark', fig_path=fig_path)


def Venus_correction(baseline_fname='base_spectrum.pkl', fig_path=fig_path_bullard, save=True, **kwargs):
    R_Venus = sh.get_pysh_constants('Venus', 'r')*1e-3  # in km
    lV, SV = sh.get_psd_Venus()  # power per lm in km^2 km^2
    l, S =  sh.load_model_spectrum_pkl(fname=baseline_fname, path=fig_path, **kwargs)

    # scale to Venus RMS
    rms_V = sh.parseval_rms_2D(SV, sh.l_to_k(lV, R_Venus))
    rms_base = sh.parseval_rms(S, sh.l_to_k(l, R=2))
    S = S * (rms_V/rms_base)**2
    print('rms Venus', rms_V)

    # plot
    fig, ax = plt.subplots()
    ax.loglog(lV, SV,  '.-', lw=0.5, alpha=0.5, c='xkcd:slate', label='Venus (Wieczorek 2015)')
    ax.loglog(l, S,  '^-', lw=0.5, alpha=0.5, c='xkcd:forest green', label='2D ASPECT')
    ax.set_xlabel(r'Degree')
    ax.set_ylabel(r'PSD (km$^{2}$ km$^{2}$)')
    ax.legend()

    if save:
        plot_save(fig, fname='Venus_correction', fig_path=fig_path)

Venus_correction()
