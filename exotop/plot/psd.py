import sh_things as sh
import numpy as np
import postaspect.plt_aspect as plat
from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path_home, fig_path_home, \
    fig_fmt, regime_grid_td, load_grid, p_Earth, postprocess_kwargs, benchmark_path, data_path_bullard, fig_path_bullard
from useful_and_bespoke import dark_background
import matplotlib.pyplot as plt

data_path = data_path_home  #bullard
fig_path = fig_path_home #bullard
labelsize = 18
ticksize = 14
legsize = 14

""" set dimensionalisation factors """
R_p = 6371  # Earth
# d, dT, alpha = 600, 442, 4e-5 # Lees table 1-2: Ra=1e6
# d, dT, alpha = 2700, 3000, 2e-5  # Venus
# d, dT, alpha = 2890, 3000, 3e-5  # Hoggard AGU Monograph
# d, dT, alpha = 2700, 3000, 3e-5  # test
d, dT, alpha = 1, 1, 1
regimes_use = ['chaotic']

fig, axes = plt.subplots(2, 1, figsize=(7, 10))

""" manu - all norm spectra with fit and Venus (and Hoggard?) """
fig, *axs = plat.plot_norm_spectra(Ra_ls, eta_ls, cmap='rainbow', end_grid=end_grid, regime_grid=regime_grid_td,
                                  include_regimes=regimes_use, save=False, show_natscales=False,
                                  data_path=data_path, pend='_sph', fend='.pkl', test=False,
                                  figsize=(8, 5), z_name='case', cbar=True,
                                  show_beta_guide=False, clabelpad=30, fig=fig, ax=axes[0],
                                  labelsize=labelsize, ticksize=ticksize, marker=None, lw=1, alpha=0.4, labelpad=16,
                                  # xlim=(1e-3, 3e-2),
                                  max_dscale=2, bl_fudge=5, legsize=legsize, # c_guide='xkcd:off white',
                                  xlabel='', ylabel='Power spectral density\n'+'(\% relative to total)',
                                  x2label='', clabel='Case', #clabel=r'log(Ra$_{i, {\rm eff}})$',
                                  norm='rel_power', whole=False, dim=False, d=d, dT=dT, alpha_m=alpha, R_p=2 * d,
                                 xlim_l=(0.3, 130), x1_name='wavenumber', show_degrees=True,
                                   #vmin=6, vmax=7.2
                                )

_, _, fig, _ = sh.Venus_correction(baseline_fname='base_spectrum_l1.pkl', fig_path=fig_path, data_path=data_path,
                                    R_base=2, lmin=1, set_axlabels=False, c_fit='xkcd:bordeaux', c_Ve='xkcd:squash',#'xkcd:dark',
                                    x_name='wavenumber', ticksize=ticksize,
                                    save=False, plot=True, units='m3', scale_to=1.0, alpha=0.9, labelsize=labelsize,
                                    legsize=legsize, fig=fig, ax=axes[1])  # axs[0] if no secondary ax; this plots degrees

_, _, fig, _ = sh.Venus_correction(baseline_fname='base_spectrum_l1.pkl', fig_path=fig_path, data_path=data_path,
                                    load_fname='spectrum_-2.pkl', is_1D=True, show_orig=False, V_label=r'$k^{-2}$',
                                    R_base=2, lmin=1, set_axlabels=False, c_fit='xkcd:dark', c_Ve='xkcd:reddish orange', #'xkcd:bubblegum pink',
                                    x_name='wavenumber', marker_Ve='v', legsize=legsize, ticksize=ticksize,
                                    save=False, plot=True, units='m3', scale_to=1.0, alpha=0.9, labelsize=labelsize,
                                    fig=fig, ax=axes[1])  # axs[0] if no secondary ax; this plots degrees

axes[1].set_xlabel('Nondimensional wavenumber', fontsize=labelsize)
fig.savefig(fig_path + 'psd_stacked_k.png', bbox_inches='tight')



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


""" get all time-averaged spectra and store """

# for ii, eta in enumerate(eta_ls):  # across eta_ls
#     cases_ii = ['Ra' + Ra + '-eta' + eta + e for Ra, e in zip(Ra_ls, end_grid[ii])]
#     labels_ii = ['Ra=' + Ra for Ra in Ra_ls]
#     for jj, Ra in enumerate(Ra_ls):
#         if regime_grid_td[ii, jj] in regimes_use:
#             case = cases_ii[jj]
#             t1 = t1_grid[ii, jj]
#             print('Calculating spectrum for', case)
#             fig, ax = sh.dct_spectrum_avg(case, L_x=8,
#                                           dim=False, R_p=d, d=d, dT=dT, alpha=alpha,
#                                           t0=t1, x_res=1, t_res=1,
#                                           test=False, data_path=data_path, fig_path=fig_path,
#                                           check_norm=False,
#                                           plot=True, load=False, dump=True, save=True, y0_guide=1e0,
#                                           )
#             print('    ...finished!')
