import numpy as np
import sh_things as sh
import model_1D.evolve as evol
from datetime import date
from matplotlib import rc
import matplotlib.pyplot as plt
from useful_and_bespoke import dark_background, cmap_from_ascii
import pyshtools


rc('text', usetex=True)  # turn off for running over ssh
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'CMU Serif'
today = date.today().strftime("%b-%d-%Y")
spectrum_fname='base_spectrum_l1.pkl'
spectrum_fpath='/home/claire/Works/exo-top/exotop/top_spectra/'
rho_w=1000
cmap_f = 'mars_1'  #'os250k-metres'
cmap = cmap_from_ascii(name=cmap_f, path='/home/claire/Works/exo-top/exotop/plot/cmaps/')  #.reversed()
cmap.set_bad('xkcd:light teal')  # flooded val

pl = evol.build_planet_from_id(ident='baseline', run_kwargs=None, update_kwargs={'x_Eu': 0.35})
labelsize = 16 #20
water_load_ratio = pl.rho_m / (pl.rho_m - rho_w)
degree, phi0 = sh.load_model_spectrum_pkl(fname=spectrum_fname, path=spectrum_fpath)
l = np.arange(len(phi0))
k = sh.l_to_k(l, R=2)  # original model spectrum uses R_b = 2d = 2
h_rms0 = sh.parseval_rms(phi0, k)
h_ratio = pl.dyn_top_rms[-1] / h_rms0
clm = sh.random_harms_from_psd(phi0, l, R=2, h_ratio=h_ratio, plot=False, verbose=False)

clm = pyshtools.datasets.Venus.VenusTopo719(
    lmax=719)  # 719 degree and order spherical harmonic model of the shape of the planet Venus (Wieczorek 2015).



dat, fig, ax = sh.coeffs_to_grid(clm, R=2, scale_to_1D=False, zscale=1e-3, plot_grid=True, plot_spectrum=False,
                                     cmap=cmap, cline='k', lw=3, figsize=(5, 3), labelsize=labelsize, ticksize=16, vmin=-1, vmax=1,
                                     # flood=fl,
                                     verbose=False, save=True, cbar=True, clabel='Dynamic topography (km)')

print('max', np.max(dat.data)*1e-3, 'km')

# n = 1 #20
# for ii, fl in enumerate(np.linspace(-1, 1, num=n)):  #, -0.5, 0, 0.5, 0.8]):
#     dat, fig, ax = sh.coeffs_to_grid(clm, R=2, scale_to_1D=False, zscale=1e-3, plot_grid=True, plot_spectrum=False,
#                                      cmap=cmap, cline='k', lw=3, figsize=(5, 3), labelsize=labelsize, ticksize=16, vmin=-1, vmax=1,
#                                      # flood=fl,
#                                      verbose=False, save=False, cbar=True, clabel='Dynamic topography (km)')
#     # Get the images on an axis
#     im = ax.images
#
#     # Assume colorbar was plotted last one plotted last
#     cb = im[-1].colorbar
#     # cb.outline.set_edgecolor('xkcd:off white')
#
#     ax.contour(dat.lons(), dat.lats(), dat.data, colors='0.1', linewidths=1)
#
#     plt.subplots_adjust(wspace=0.3)
#     # fig, *axes = dark_background(fig, [ax, cb.ax])
#     fname = 'topo_gridw_'
#     fig.savefig('/home/claire/Works/exo-top/exotop/figs_scratch/' + fname + str(ii) + '.png', dpi=400, bbox_inches='tight',
#                 facecolor=fig.get_facecolor())
plt.show()