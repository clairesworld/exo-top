import sh_things as sh
import numpy as np
import postaspect.plt_aspect as plat
from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path_home, fig_path_home, \
    fig_fmt, regime_grid_td, load_grid, p_Earth, postprocess_kwargs, benchmark_path, data_path_bullard, fig_path_bullard
from useful_and_bespoke import dark_background, colourised_legend, colorize
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.pyplot import rcParams
from useful_and_bespoke import dark_background, cmap_from_ascii

rc('text', usetex=True)  # turn off for running over ssh
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'CMU Serif'

data_path = data_path_home
fig_path = fig_path_home
psd_path = '/home/claire/Works/exo-top/exotop/top_spectra/'
labelsize = 18
ticksize = 14
legsize = 14

xlim = (1e-7, 1e-5)  #None # (6e-1, 6e1)
ylim0 = (1e-6, 1e1)
ylim1 = (3e-4, 1e1)
model_c = 'xkcd:bordeaux'  # (0.02156863, 0.68274886, 0.93022931)
regimes_use = ['chaotic']
norm = 'rms'  # 'rms'

cmap_f = 'mars_1'  #'os250k-metres'
cmap = cmap_from_ascii(name=cmap_f, path='/home/claire/Works/exo-top/exotop/plot/cmaps/')

" plot spectrum "
def spec_plot(ax=None, fig=None, lmax=None, c_spec='k'):
    l, S = sh.load_model_spectrum_pkl(fname='base_spectrum.pkl', path='/home/claire/Works/exo-top/exotop/top_spectra/',
                                      data_path=data_path, case=None)

    fig, ax = sh.plot_norm_psd(baseline_fname='base_spectrum_l1.pkl', fig_path=fig_path, data_path=data_path,
                               psd_path=psd_path, lmax=lmax,
                               R=6000e3, lmin=1, c=c_spec,  # 'xkcd:bordeaux',
                               x_name='wavenumber', ticksize=ticksize, xlim=xlim, ylim=None,
                               save=False, labelsize=labelsize, legend=True, label='Model dynamic topography',
                               show_degrees=True, x2label='Spherical harmonic degree, $l$', marker='o',
                               ylabel='Power spectral density (m$^3$)',
                               norm=norm, rms1=1000, legsize=legsize, fig=fig, ax=ax)
    ax.set_xlabel('Wavenumber (m$^{-1}$)', fontsize=labelsize)
    return fig, ax


def map_plot(ax=None, fig=None, lmax=None):
    degree, phi0 = sh.load_model_spectrum_pkl(fname='base_spectrum.pkl', path='/home/claire/Works/exo-top/exotop/top_spectra/',)
    if lmax is None:
        lmax = len(phi0)
    l = np.arange(lmax)
    phi = phi0[:lmax]
    k = sh.l_to_k(l, R=2)  # original model spectrum uses R_b = 2d = 2
    h_rms0 = sh.parseval_rms(phi, k)
    h_ratio = 1000 / h_rms0
    clm = sh.random_harms_from_psd(phi, l, R=2, h_ratio=h_ratio, plot=False, verbose=False)
    dat, fig, ax = sh.coeffs_to_grid(clm, R=2, scale_to_1D=False, zscale=1e-3, plot_grid=True, plot_spectrum=False,
                                     cmap=cmap, cline='k', lw=3, figsize=(5, 3), labelsize=labelsize, ticksize=16,
                                     vmin=-1, vmax=1, fig=fig, ax=ax,
                                     # flood=fl,
                                     verbose=False, save=False, cbar=False, clabel='Dynamic topography (km)')
    return fig, ax


def animate_expansion(l_max=40, l_min=2):
    for ii in range(2, l_max):
        fig, axes = plt.subplots(1, 2, figsize=(10, 3))
        fig, axes[0] = spec_plot(ax=axes[0], fig=fig, c_spec='k', lmax=ii)
        fig, axes[1] = map_plot(ax=axes[1], fig=fig, lmax=ii)

        plt.subplots_adjust(wspace=0.4)
        fig.savefig(fig_path + 'clever_map' + str(ii) + '.png', bbox_inches='tight')
# plt.show()

animate_expansion(l_min=10, l_max=40)