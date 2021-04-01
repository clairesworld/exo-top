import sh_things as sh
import numpy as np
import postaspect.plt_aspect as plat
from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, \
    fig_fmt, regime_grid_td, load_grid, p_Earth, postprocess_kwargs, benchmark_path

""" set dimensionalisation factors """
R_p = 6371
# d, dT, alpha = 600, 442, 4e-5 # Lees table 1-2: Ra=1e6
# d, dT, alpha = 2700, 3000, 2e-5  # Venus
d, dT, alpha = 2890, 3000, 3e-5  # Hoggard AGU Monograph
# d, dT, alpha = 2700, 3000, 3e-5  # test

regimes_use = ['chaotic']

""" plot normalised spectra on single axis """

fig, ax = plat.plot_norm_spectra(Ra_ls, eta_ls, cmap='rainbow', end_grid=end_grid, regime_grid=regime_grid_td,
                                 include_regimes=regimes_use, save=True,
                                 data_path=data_path, pend='_sph', fend='.pkl', figname='h_spectra_stacked_dim',
                                 fig=None, ax=None, figsize=(8, 5), z_name='Ra_i_eff', cbar=True,
                                 labelsize=16, ticksize=12, marker='.', lw=0.5,
                                 xlabel=r'Wavenumber (km$^{-1}$)', ylabel='Normalised power spectral density',
                                 x2label='Spherical harmonic degree', clabel=r'log(Ra$_{i, {\rm eff}}$)',
                                 norm='intercept', dim=True, d=d, dT=dT, alpha_m=alpha, R_p=d,
                                 # add_files=[benchmark_path + 'lees_topo_grids/psd_hoggard.csv'], add_label=['Hoggard+ (2016)']
                                 )

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
