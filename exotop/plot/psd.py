import sh_things as sh
import numpy as np
import postaspect.plt_aspect as plat
from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path_home, fig_path_home, \
    fig_fmt, regime_grid_td, load_grid, p_Earth, postprocess_kwargs, benchmark_path, data_path_bullard, fig_path_bullard
from useful_and_bespoke import dark_background, colourised_legend, colorize
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.pyplot import rcParams

rc('text', usetex=True)  # turn off for running over ssh
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'CMU Serif'

data_path = data_path_home
fig_path = fig_path_home
psd_path = '/home/claire/Works/exo-top/exotop/top_spectra/'
labelsize = 18
ticksize = 14
legsize = 14
cmap = 'gray'  #'rainbow'
norm = 'rms'  # 'rms'
ylabel = 'Normalised power spectral density (m$^3$)'  # \n' + '(% relative to total)',

""" set dimensionalisation factors """
R_p = 6371  # Earth
# d, dT, alpha = 600, 442, 4e-5 # Lees table 1-2: Ra=1e6
# d, dT, alpha = 2700, 3000, 2e-5  # Venus
# d, dT, alpha = 2890, 3000, 3e-5  # Hoggard AGU Monograph
# d, dT, alpha = 2700, 3000, 3e-5  # test
d, dT, alpha = 1, 1, 1

fig, axes = plt.subplots(2, 1, figsize=(7, 10))
xlim = (6e-1, 6e1)
ylim0 = (1e-6, 1e1)
ylim1 = (3e-4, 1e1)
model_c = 'xkcd:bordeaux'  # (0.02156863, 0.68274886, 0.93022931)
regimes_use = ['chaotic']

""" manu - all norm spectra with fit and Venus (and Hoggard?) """
fig, *axs = plat.plot_norm_spectra(Ra_ls, eta_ls, end_grid=end_grid, regime_grid=regime_grid_td,
                                   include_regimes=regimes_use, save=False, show_natscales=False, cmap=cmap,
                                   data_path=data_path, pend='_sph', fend='.pkl', test=False,
                                   figsize=(8, 5), z_name='case', cbar=False,
                                   show_beta_guide=False, clabelpad=30, fig=fig, ax=axes[0],
                                   labelsize=labelsize, ticksize=ticksize, marker=None, lw=1, alpha=0.4, labelpad=16,
                                   xlim=xlim, ylim=ylim0,
                                   max_dscale=2, bl_fudge=5, legsize=legsize,  # c_guide='xkcd:off white',
                                   xlabel='', ylabel='',
                                   x2label='Spherical harmonic degree', clabel='Case',
                                   # clabel=r'log(Ra$_{i, {\rm eff}})$',
                                   norm=norm, whole=False, dim=False, d=d, dT=dT, alpha_m=alpha, R_p=2 * d,
                                   # xlim_l=(0.3, 130),
                                   x1_name='wavenumber', show_degrees=False, legend=True,
                                   # vmin=1, vmax=11,
                                   )
# cases = np.arange(1, 11 + 1)
# ax = colourised_legend(axes[0], clist=colorize(cases, cmap=cmap)[0],
#                        cleglabels=cases, lw=1, ls='-', marker=None,
#                        markersize=0,
#                        legsize=legsize, ncol=1, title='Case')


fig, ax = sh.plot_norm_psd(baseline_fname='base_spectrum_l1.pkl', fig_path=fig_path, data_path=data_path, psd_path=psd_path,
                           R=2, lmin=1, c=model_c ,  # 'xkcd:bordeaux',
                           x_name='wavenumber', ticksize=ticksize, xlim=xlim, ylim=ylim0,
                           save=False, labelsize=labelsize, legend=True, label='Model dynamic topography',
                           show_degrees=True, x2label='Spherical harmonic degree, $l$', marker='^',
                           ylabel='',
                           norm=norm, legsize=legsize, fig=fig, ax=axes[0])

fig, ax = sh.plot_norm_psd(baseline_fname='base_spectrum_l1.pkl', fig_path=fig_path, data_path=data_path, psd_path=psd_path,
                           R=2, lmin=1, c=model_c ,  # 'xkcd:bordeaux',
                           x_name='wavenumber', ticksize=ticksize, xlim=(6e-1, 6e1), ylim=ylim1,
                           show_degrees=False, save=False, labelsize=labelsize, legend=False,
                           label='Model dynamic topography', marker='^', lmax=53,
                           norm=norm, legsize=legsize, fig=fig, ax=axes[1])

fig, ax = sh.plot_norm_psd(baseline_fname='Venus', fig_path=fig_path, data_path=data_path, psd_path=psd_path,
                           R=2, lmin=1, c='xkcd:squash', lmax=53,
                           x_name='wavenumber', ticksize=ticksize, xlim=xlim, ylim=ylim1,
                           show_degrees=True, save=False, labelsize=labelsize, legend=False,
                           label='Venus (Wieczorek 2015)', marker='v',
                           norm=norm, legsize=legsize, fig=fig, ax=axes[1])

fig, ax = sh.plot_norm_psd(baseline_fname='spectrum_-2.pkl', fig_path=fig_path, data_path=data_path, psd_path=psd_path,
                           R=2, lmin=1, c='xkcd:reddish orange', lmax=53,
                           x_name='wavenumber', ticksize=ticksize, xlim=xlim, ylim=ylim1,
                           show_degrees=True, save=False, labelsize=labelsize, legend=True,
                           label=r'Red noise, $k^{-2}$', marker='o', ylabel='',
                           norm=norm, legsize=legsize, fig=fig, ax=axes[1])

# _, _, fig, _ = sh.Venus_correction(baseline_fname='base_spectrum_l1.pkl', fig_path=fig_path, data_path=data_path,
#                                     R_base=2, lmin=1, set_axlabels=False, c_fit='xkcd:bordeaux', c_Ve='xkcd:squash',#'xkcd:dark',
#                                     x_name='wavenumber', ticksize=ticksize, xlim=(6e-1, 3e1),
#                                     save=False, plot=True, units='m3', scale_to=1.0, alpha=0.9, labelsize=labelsize,
#                                     legsize=legsize, fig=fig, ax=axes[1])  # axs[0] if no secondary ax; this plots degrees
#
# _, _, fig, _ = sh.Venus_correction(baseline_fname='base_spectrum_l1.pkl', fig_path=fig_path, data_path=data_path,
#                                     load_fname='spectrum_-2.pkl', is_1D=True, show_orig=False, V_label=r'$k^{-2}$',
#                                     R_base=2, lmin=1, set_axlabels=False, c_fit='xkcd:dark', c_Ve='xkcd:reddish orange', #'xkcd:bubblegum pink',
#                                     x_name='wavenumber', marker_Ve='v', legsize=legsize, ticksize=ticksize, xlim=(6e-1, 3e1),
#                                     save=False, plot=True, units='m3', scale_to=1.0, alpha=0.9, labelsize=labelsize,
#                                     fig=fig, ax=axes[1])  # axs[0] if no secondary ax; this plots degrees

axes[1].set_xlabel('Nondimensional wavenumber, $kd$', fontsize=labelsize)
fig.supylabel(ylabel, fontsize=labelsize)


# plt.show()
fig.savefig(fig_path + 'psd_stacked_k.pdf', bbox_inches='tight')
#

















""" plot normalised spectra relative power on single axis - full spectrum norm rms"""

# fig, *ax = plat.plot_norm_spectra(Ra_ls, eta_ls, cmap='rainbow', end_grid=end_grid, regime_grid=regime_grid_td,
#                                   include_regimes=regimes_use, save=False, show_natscales=True,
#                                   data_path=data_path, pend='_sph', fend='.pkl', figname='h_spectra_stacked',
#                                   fig=None, ax=None, figsize=(8, 5), z_name='Ra_i_eff', cbar=True,
#                                   show_beta_guide=False,
#                                   labelsize=16, ticksize=12, marker='.', lw=0.5, alpha=0.8, labelpad=16,
#                                   # xlim=(1e-3, 3e-2),
#                                   max_dscale=2, bl_fudge=5,  # c_guide='xkcd:off white',
#                                   xlabel=None, ylabel='Normalised power spectral density',
#                                   x2label='Spherical harmonic degree', clabel=r'log(Ra$_{i, eff}$)',
#                                   norm=None, whole=True, dim=False, d=d, dT=dT, alpha_m=alpha, R_p=2 * d,
#                                   relative_power=True,
#                                   )
# fig.savefig(fig_path + 'psd_relative_stacked.png', bbox_inches='tight')


""" plot normalised spectra on single axis - full spectrum norm rms"""

# fig, *ax = plat.plot_norm_spectra(Ra_ls, eta_ls, cmap='rainbow', end_grid=end_grid, regime_grid=regime_grid_td,
#                                   include_regimes=regimes_use, save=False, show_natscales=True,
#                                   data_path=data_path, pend='_sph', fend='.pkl', figname='h_spectra_stacked',
#                                   fig=None, ax=None, figsize=(8, 5), z_name='Ra_i_eff', cbar=True,
#                                   show_beta_guide=False,
#                                   labelsize=16, ticksize=12, marker='.', lw=0.5, alpha=0.8, labelpad=16,
#                                   # xlim=(1e-3, 3e-2),
#                                   max_dscale=2, bl_fudge=5,  # c_guide='xkcd:off white',
#                                   xlabel=None, ylabel='Normalised power spectral density',
#                                   x2label='Spherical harmonic degree', clabel=r'log(Ra$_{i, eff}$)',
#                                   norm='rms', whole=True, dim=False, d=d, dT=dT, alpha_m=alpha, R_p=2 * d,
#                                   )
# fig.savefig(fig_path + 'psd_stacked.png', bbox_inches='tight')

""" plot normalised spectra on single axis"""

# fig, *ax = plat.plot_norm_spectra(Ra_ls, eta_ls, cmap='rainbow', end_grid=end_grid, regime_grid=regime_grid_td,
#                                   include_regimes=regimes_use, save=True, show_natscales=True,
#                                   data_path=data_path, pend='_sph', fend='.pkl', figname='psd_stacked_sub',
#                                   fig=None, ax=None, figsize=(8, 5), z_name='Ra', cbar=True,
#                                   show_beta_guide=False,
#                                   labelsize=16, ticksize=12, marker='.', lw=0.5, alpha=0.4, labelpad=16,
#                                   # xlim=(1e-3, 3e-2),
#                                   max_dscale=1, bl_fudge=5,  # c_guide='xkcd:off white',
#                                   xlabel=None, ylabel='Normalised power spectral density',
#                                   x2label='Spherical harmonic degree', clabel=r'log(Ra$_{i, eff}$)',
#                                   norm='intercept', whole=False, dim=False, d=d, dT=dT, alpha_m=alpha, R_p=2 * d,
#                                   )


""" plot normalised spectra on single axis - slides, norm intercept """

# fig, *ax = plat.plot_norm_spectra(Ra_ls, eta_ls, cmap='rainbow', end_grid=end_grid, regime_grid=regime_grid_td,
#                                  include_regimes=regimes_use, save=False,
#                                  data_path=data_path, pend='_sph', fend='.pkl', figname='h_spectra_stacked_slides',
#                                  fig=None, ax=None, figsize=(8, 5), z_name='Ra_i_eff', cbar=True, show_beta_guide=True,
#                                  labelsize=16, ticksize=12, marker='.', lw=0.5, alpha=0.8, labelpad=16,
#                                  xlim=(1e-3, 3e-2), max_dscale=1, bl_fudge=5, # c_guide='xkcd:off white',
#                                  xlabel=None, ylabel='Normalised power spectral density',
#                                  x2label='Spherical harmonic degree', clabel=r'log(Ra$_{i, eff}$)',
#                                  norm='intercept', dim=True, d=d, dT=dT, alpha_m=alpha, R_p=d,
#                                  # add_files=[benchmark_path + 'lees_topo_grids/psd_hoggard.csv'], add_label=['Hoggard+ (2016)']
#                                  )
#
# fig, *axes = dark_background(fig, ax)
# fig.savefig(fig_path + 'psd_stacked_slides.png', bbox_inches='tight', transparent=True)

""" what does a cartesian projection look like? """

# base_case = 'Ra1e8-eta1e7-wide'
# l, S = sh.make_baseline_spectrum(base_case, R=1, data_path=data_path, fig_path=fig_path, newfname='base_spectrum',
#                                  max_dscale=1, bl_fudge=5,)


