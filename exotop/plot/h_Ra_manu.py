""" ASPECT runs: plot scalings of h vs. Ra or T-heuristic, with subplots by eta """

import matplotlib.pyplot as plt
from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    c_regimes_td, fig_fmt, regime_grid_td, regime_names_td, load_grid, p_Earth, postprocess_kwargs
from postaspect import plt_aspect as plat

"""setup"""

load = load_grid
labelsize = 16
ms = 10
legsize = 12
vmin = None  # 1
vmax = None  # 3
cmap = 'rainbow'  # 'art-nouveau-03'
c_rms = ['k', 'xkcd:lime green', 'xkcd:lilac', 'xkcd:orange', 'xkcd:yellow']  # update: must correspond to entire z range i.e. eta_ls
cleglabels = [r'$\Delta \eta = 10^{5}$', r'$\Delta \eta = 10^{6}$', r'$\Delta \eta = 10^{7}$',
              r'$\Delta \eta = 10^{8}$',  r'$\Delta \eta = 10^{9}$']  # these correspond to entire z range as well

include_regimes = ['chaotic']  # , 'not ready'
averagescheme = 'timefirst'
which_x = 'Ra_i_eff'

fig, axes = plt.subplots(1, 2, figsize=(10, 5))

"""all eta on one axis"""

fig, ax0 = plat.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=load, data_path=data_path,
                   fig_path=fig_path, averagescheme=averagescheme, p_dimensionals=None, which_x=which_x,
                   beta0=[0.1, -0.1], sigma=1, fiterror=False, legend=True, fig=fig, ax=axes[0],
                   include_regimes=include_regimes, save=True, fname='h_Raieff_chaotic_timeavg', labelsize=labelsize,
                   xlabel=r'Ra$_{i,eff}$', ylabel='rms dynamic topography', legsize=legsize, cleglabels=cleglabels,
                   title='', showpeak=False, vmin=vmin, vmax=vmax,
                   cmap=None, c_rms=c_rms,
                   fit=True, logx=True, logy=True, ms=ms,
                   show_isoviscous=False, ylim=None, xlim=None, postprocess_kwargs=postprocess_kwargs,
                   regime_grid=regime_grid_td, errortype='standard', cbar=False)

"""model versus data"""

fig, ax1 = plat.plot_model_data_errorbars(Ra_ls, eta_ls, regime_grid=regime_grid_td, t1_grid=t1_grid, load_grid=load,
                                         end_grid=end_grid, legend=False, ms=ms, fig=fig, ax=axes[1],
                                         postprocess_kwargs=postprocess_kwargs, averagescheme=averagescheme,
                                         ylim=[6e-3, 2e-2], which_x=which_x, which_h='rms', data_path=data_path,
                                         clist=c_rms,
                                         cmap=None, z_name='eta', fname='model-data-chaotic_timeavg',
                                         save=True, include_regimes=include_regimes, errs=[0.5, 0.2, 0.1, 0.05],
                                         fig_fmt=fig_fmt, vmin=vmin, vmax=vmax,
                                         show_cbar=False, errortype='standard',
                                         title='',
                                         ylabel=r'Model $\Delta h_{rms}^\prime$',
                                         xlabel=r'Data $\Delta h_{rms}^\prime$')

fig.suptitle(r'Fit to $C$ Ra$_{i,eff}^p$', fongsize=labelsize)
fig.savefig(fig_path + 'h_Ra_scalings.png', bbox_inches='tight')

