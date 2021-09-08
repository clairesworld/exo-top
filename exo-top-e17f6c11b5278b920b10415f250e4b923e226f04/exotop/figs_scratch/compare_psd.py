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



lV, phiV_sc, fig, ax = sh.Venus_correction(baseline_fname='base_spectrum_l1.pkl', fig_path=fig_path, R_base=2, lmin=1,
                     save=True, plot=True, units='km4', scale_to='Venus', alpha=0.9, labelsize=16)

fig, ax = dark_background(fig, ax)
# sh.make_Venus_reference(newfname='Venus_spectrum_l1.pkl', baseline_fname='base_spectrum_l1.pkl',
#                      fig_path=fig_path)
plt.show()
