""" ASPECT runs: plot scalings of h vs. Ra or T-heuristic, with subplots by eta """

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
figsize = (7, 5)  # (5, 5)
include_regimes = ['chaotic']  # , 'not ready'
averagescheme = 'timefirst'
which_x = 'Ra_i_eff'


"""model versus data"""

fig, ax = plat.plot_model_data_errorbars(Ra_ls, eta_ls, regime_grid=regime_grid_td, t1_grid=t1_grid, load_grid=load,
                                         end_grid=end_grid, literature_file=None, legend=False, ms=ms,
                                         postprocess_kwargs=postprocess_kwargs, averagescheme=averagescheme,
                                         ylim=[6e-3, 2e-2], which_x=which_x, which_h='rms', data_path=data_path,
                                         # clist=['xkcd:lime green', 'xkcd:lilac', 'xkcd:orange', 'xkcd:yellow'],
                                         cmap=cmap, z_name='eta', fname='model-data-chaotic_timeavg',
                                         save=True, include_regimes=include_regimes, errs=[0.5, 0.2, 0.1, 0.05],
                                         fig_fmt=fig_fmt, vmin=vmin, vmax=vmax,
                                         show_cbar=False, figsize=figsize, errortype='standard',
                                         title=r'Fit to $C$ Ra$_{i,eff}^p$',
                                         ylabel=r'Model $\Delta h_{rms}^\prime$',
                                         xlabel=r'Data $\Delta h_{rms}^\prime$')

"""all eta on one axis"""

_ = plat.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=load, data_path=data_path,
                   fig_path=fig_path, averagescheme=averagescheme, p_dimensionals=None, which_x=which_x,
                   beta0=[0.1, -0.15], sigma=1, fiterror=False, legend=True,
                   include_regimes=include_regimes, save=True, fname='h_Raieff_chaotic_timeavg', labelsize=labelsize,
                   xlabel=r'Ra$_{i,eff}$', ylabel='rms dynamic topography', legsize=legsize, cleglabels=cleglabels,
                   title=r'Fit to $C$ Ra$_{i,eff}^p$', showpeak=False, vmin=vmin, vmax=vmax,
                   cmap=None, c_rms=c_rms,
                   fit=True, logx=True, logy=True, ms=ms,
                   show_isoviscous=False, ylim=None, xlim=None, postprocess_kwargs=postprocess_kwargs,
                   regime_grid=regime_grid_td, figsize=figsize, errortype='standard', cbar=False)

_ = plat.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=load, data_path=data_path,
                   fig_path=fig_path, averagescheme=averagescheme, p_dimensionals=None, which_x=which_x,
                   beta0=[0.1, -0.15], sigma=1, fiterror=False, legend=True,
                   include_regimes=include_regimes, save=True, fname='h_Raieff_chaotic_timeavg', labelsize=labelsize,
                   xlabel=r'Ra$_{i,eff}$', ylabel='peak dynamic topography', legsize=legsize, cleglabels=cleglabels,
                   title=r'Fit to $C$ Ra$_{i,eff}^p$', showpeak=False, vmin=vmin, vmax=vmax,
                   cmap=None, c_rms=c_rms,
                   fit=True, logx=True, logy=True, ms=ms,
                   show_isoviscous=False, ylim=None, xlim=None, postprocess_kwargs=postprocess_kwargs,
                   regime_grid=regime_grid_td, figsize=figsize, errortype='standard', cbar=False,
                   which_h='peak')

# _ = plat.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=load, data_path=data_path,
#                    fig_path=fig_path, averagescheme='timefirst', p_dimensionals=None, which_x='Ra',
#                    beta0=[0.1, -0.15], sigma=1, fiterror=False,
#                    include_regimes=['chaotic'], save=True, fname='h_Ra_chaotic_timeavg', labelsize=labelsize, legend=True,
#                    xlabel=r'Ra$', ylabel='dynamic topography',
#                    title=r'Fit to $C$ Ra$^p$', showpeak=False, vmin=vmin, vmax=vmax,
#                    cmap=cmap, c_rms=None, fit=True, logx=True, logy=True, hscale=1, ms=ms,
#                    show_isoviscous=False, ylim=None, xlim=None, postprocess_kwargs=postprocess_kwargs,
#                    regime_grid=regime_grid_td, figsize=(5, 5), errortype='standard')


# Ra_F
# _ = plat.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=load, data_path=data_path,
#                    fig_path=fig_path, averagescheme='timefirst', p_dimensionals=None, which_x='Ra_F_eff',
#                    beta0=[0.1, -0.15], sigma=1, fiterror=False,
#                    include_regimes=['chaotic'], save=True, fname='h_RaF_chaotic_timeavg', labelsize=16, legend=True,
#                    xlabel=r'Ra$_{F,eff}$', ylabel='dynamic topography',
#                    title=r'fit to CRa$_{F,eff}^p$, averaging time first', showpeak=False,
#                    c_peak='xkcd:forest green', c_rms='xkcd:periwinkle', fit=True, logx=True, logy=True, hscale=1,
#                    show_isoviscous=False, ylim=None, xlim=None, postprocess_kwargs=postprocess_kwargs,
#                    regime_grid=regime_grid_td)

# try h scalings to two-component power law for chaotic
#
# _ = plat.plot_h_vs_2component(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=load, data_path=data_path,
#                  fig_path=fig_path, averagescheme='timefirst', p_dimensionals=None, which_xs=('Ra_i_eff', 'eta'),
#                  include_regimes=['chaotic'], save=True, fname='h_Ra-eta_chaotic_timeavg', labelsize=16, legend=True,
#                  title=r'fit to C Ra$^p \Delta\eta^q$', xlabel=r'Ra$_i$', ylabel='dynamic topography',
#                  fit=True, logx=True, logy=True, hscale=1, clabel=r'$\Delta\eta$', showpeak=False,
#                  show_isoviscous=False, ylim=(5e-3, 3e-2), xlim=(1e7, 3e8),  postprocess_kwargs=postprocess_kwargs,
#                  regime_grid=regime_grid_td)


# scalings with various Ra, average time first

# fig, ax = plat.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=load, data_path=data_path,
#                  fig_path=fig_path, averagescheme='timefirst', p_dimensionals=None, which_x='Ra_i', ms=7,
#                  beta0=[0.1, -0.15],  sigma=2, showpeak=False,
#                  include_regimes=['steady', 'trans', 'chaotic'], save=False, fname='h_Rai_chaotic_timeavg', labelsize=16,
#                  xlabel=r'Ra$_i$', ylabel='dynamic topography',
#                  title=r'fit to CRa$_i^n$, averaging time first', legend=False,
#                  c_rms='xkcd:periwinkle', fit=True, logx=True, logy=True, hscale=1,
#                  show_isoviscous=False, ylim=None, xlim=None, postprocess_kwargs=postprocess_kwargs,
#                  regime_grid=regime_grid_td)
#
# # add literature comparisons
# alpha = 3e-5
# dT = 2390 - 255
# d = 2890e3
# h = [x / (alpha*dT*d) for x in [7.28e3, 6.53e3, 4.41e3, 3.47e3]]
# Ra_i = [1e5, 1e5, 1e6, 1e7]
# ax.plot(Ra_i, h, 'o', c='g', label='Arnould+ 2018')
# ax.legend()
# fig.savefig(fig_path+'h_Rai_arnould'+fig_fmt, bbox_inches='tight')


#
## h scalings with heuristic all chaotic cases
# #
# _ = plat.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=load, data_path=data_path,
#                  fig_path=fig_path, averagescheme='timefirst', p_dimensionals=None, which_x='h_components',
#                  include_regimes=['chaotic'], save=True, fname='h_T_chaotic_timeavg', labelsize=16, legend=True,
#                  xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography',
#                  title=r'fit to $C\alpha \delta_{rh} \Delta T_{rh}$, averaging time first', showpeak=False,
#                  c_peak='xkcd:forest green', c_rms='xkcd:periwinkle', fit=True, logx=True, logy=True, hscale=1,
#                  show_isoviscous=False, ylim=None, xlim=None, postprocess_kwargs=postprocess_kwargs,
#                  regime_grid=regime_grid_td, fiterror=False,)
# # #
# _ = plat.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=load, data_path=data_path,
#                  fig_path=fig_path, averagescheme='timelast', p_dimensionals=None, which_x='h_components',
#                  include_regimes=['chaotic'], save=True, fname='h_T_chaotic', labelsize=16, legend=True,
#                  xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography',
#                  title=r'fit to $C\alpha \delta_{rh} \Delta T_{rh}$, averaging time last',
#                  c_peak='xkcd:forest green', c_rms='xkcd:periwinkle', fit=True, logx=True, logy=True, hscale=1,
#                  show_isoviscous=False, ylim=None, xlim=None, postprocess_kwargs=postprocess_kwargs,
#                  regime_grid=regime_grid_td)

# _ = plat.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=load, data_path=data_path,
#                  fig_path=fig_path, averagescheme=None, p_dimensionals=None, which_x='h_components',
#                  include_regimes=['chaotic'], save=True, fname='h_T_chaotic_all', labelsize=16, legend=True,
#                  xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography',
#                     title=r'fit to $C\alpha \delta_{rh} \Delta T_{rh}$, no time-averaging',
#                  c_peak='xkcd:forest green', c_rms='xkcd:periwinkle', fit=True, logx=True, logy=True, hscale=1,
#                  show_isoviscous=False, ylim=None, xlim=None, postprocess_kwargs=postprocess_kwargs,
#                  regime_grid=regime_grid_td)

# subplots across delta eta

# plot h scalings - with dT_m*delta*alpha

# _ = plat.subplots_topo_regimes(Ra_ls, eta_ls, regime_grid_td, regime_names_td, c_regimes=c_regimes_td, save=True, t1_grid=t1_grid,
#                          averagescheme='timefirst', legloc='upper right', which_x='h_components',
#                          load_grid=load, fig_path=fig_path, fname='h_T_timeavg', fig_fmt=fig_fmt, end_grid=end_grid,
#                          labelsize=14, xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography, $h^\prime$',
#                          xlabelpad=8, ylabelpad=20, fit=True, showallscatter=False,
#                          #xlim=(1e-8, 0.9e-6), ylim=(6e-3, 10e-2),
#                          logx=True, logy=True,
#                          regimes_title='Stationarity', leftleg_bbox=(-0.1, 0.95), data_path=data_path,
#                          postprocess_kwargs=postprocess_kwargs, #include_regimes=['chaotic']
#                               )

# heuristic with no averaging

# _ = plat.subplots_topo_regimes(Ra_ls, eta_ls, regime_grid_td, regime_names_td, c_regimes=c_regimes_td, save=True, t1_grid=t1_grid,
#                          which_x='h_components', averagescheme=None, legloc='upper right',
#                          load_grid=load,
#                          fig_path=fig_path, fname='h_T_all_scatter', fig_fmt=fig_fmt, end_grid=end_grid,
#                          labelsize=14, xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography $h^\prime$',
#                          xlabelpad=8, ylabelpad=20, fit=True, showallscatter=True,
#                          #xlim=(1e-8, 0.9e-6), ylim=(6e-3, 10e-2),
#                          logx=True, logy=True,
#                          regimes_title='Stationarity', leftleg_bbox=(-0.1, 0.95), data_path=data_path,
#                          postprocess_kwargs=postprocess_kwargs,)

# # plot h scalings with Ra
#
# _ = plat.subplots_topo_regimes(Ra_ls, eta_ls, regime_grid_td, regime_names_td, c_regimes=c_regimes_td, save=True, t1_grid=t1_grid,
#                          load_grid=load, show_isoviscous=True, averagescheme='timefirst',
#                          fig_path=fig_path, fname='h_Ra_timeavg', fig_fmt=fig_fmt, end_grid=end_grid, labelsize=14,
#                          xlabel='Ra_i', which_x='Ra_i',
#                          ylabel='dynamic topography $h^\prime$', y2label='dynamic topography $h$ (km)',
#                          xlabelpad=8, ylabelpad=5, fit=True, showallscatter=False,
#                          xlim=(0.6e6, 5e8), #ylim=(1, 12),
#                          logx=True, logy=True,
#                          p_dimensionals=p_Earth,
#                          regimes_title='Stationarity', leftleg_bbox=(-0.01, 0.95), data_path=data_path,
#                          postprocess_kwargs=postprocess_kwargs,)
