""" ASPECT runs: plot scalings of h vs. Ra or T-heuristic, with subplots by eta """

import sys

sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    c_regimes_td, fig_fmt, \
    regime_grid_td, regime_names_td, load_grid, p_Earth, postprocess_kwargs  # noqa: E402
from exotop.postaspect import aspect_scalings as sc  # noqa: E402

load = True #load_grid

# try h scalings to two-component power law for chaotic
#
# _ = sc.plot_h_vs_2component(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=load, data_path=data_path,
#                  fig_path=fig_path, averagescheme='timefirst', p_dimensionals=None, which_xs=('Ra_i', 'eta'),
#                  include_regimes=['chaotic'], save=True, fname='h_Ra-eta_chaotic_timeavg', labelsize=16, legend=True,
#                  title=r'fit to C Ra$^m \Delta\eta^n$', xlabel=r'Ra$_i$', ylabel='dynamic topography',
#                  fit=True, logx=True, logy=True, hscale=1, clabel=r'$\Delta\eta$',
#                  show_isoviscous=False, ylim=(5e-3, 3e-2), xlim=(1e7, 3e8),  postprocess_kwargs=postprocess_kwargs,
#                  regime_grid=regime_grid_td)
# #
# #
_ = sc.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=load, data_path=data_path,
                 fig_path=fig_path, averagescheme='timefirst', p_dimensionals=None, which_x='Ra_i',
                 beta0=[0.1, -0.15],  sigma=2,
                 include_regimes=['chaotic'], save=True, fname='h_Rai_chaotic_timeavg', labelsize=16, legend=True,
                 xlabel=r'Ra$_i$', ylabel='dynamic topography',
                 title=r'fit to CRa$_i^n$, averaging time first',
                 c_peak='xkcd:forest green', c_rms='xkcd:periwinkle', fit=True, logx=True, logy=True, hscale=1,
                 show_isoviscous=False, ylim=None, xlim=None, postprocess_kwargs=postprocess_kwargs,
                 regime_grid=regime_grid_td)
#
_ = sc.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=load, data_path=data_path,
                 fig_path=fig_path, averagescheme='timefirst', p_dimensionals=None, which_x='Ra_i_eff',
                 beta0=[0.1, -0.15], sigma=2,
                 include_regimes=['chaotic'], save=True, fname='h_Raieff_chaotic_timeavg', labelsize=16, legend=True,
                 xlabel=r'Ra$_{i,eff}$', ylabel='dynamic topography',
                 title=r'fit to CRa$_{i,eff}^n$, averaging time first',
                 c_peak='xkcd:forest green', c_rms='xkcd:periwinkle', fit=True, logx=True, logy=True, hscale=1,
                 show_isoviscous=False, ylim=None, xlim=None, postprocess_kwargs=postprocess_kwargs,
                 regime_grid=regime_grid_td)

#
# # try h scalings with heuristic all chaotic cases
#
# _ = sc.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=load, data_path=data_path,
#                  fig_path=fig_path, averagescheme='timefirst', p_dimensionals=None, which_x='h_components',
#                  include_regimes=['chaotic'], save=True, fname='h_T_chaotic_timeavg', labelsize=16, legend=True,
#                  xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography',
#                  title=r'fit to $C\alpha \delta_{rh} \Delta T_{rh}$, averaging time first',
#                  c_peak='xkcd:forest green', c_rms='xkcd:periwinkle', fit=True, logx=True, logy=True, hscale=1,
#                  show_isoviscous=False, ylim=None, xlim=None, postprocess_kwargs=postprocess_kwargs,
#                  regime_grid=regime_grid_td)
#
# _ = sc.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=load, data_path=data_path,
#                  fig_path=fig_path, averagescheme='timelast', p_dimensionals=None, which_x='h_components',
#                  include_regimes=['chaotic'], save=True, fname='h_T_chaotic', labelsize=16, legend=True,
#                  xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography',
#                  title=r'fit to $C\alpha \delta_{rh} \Delta T_{rh}$, averaging time last',
#                  c_peak='xkcd:forest green', c_rms='xkcd:periwinkle', fit=True, logx=True, logy=True, hscale=1,
#                  show_isoviscous=False, ylim=None, xlim=None, postprocess_kwargs=postprocess_kwargs,
#                  regime_grid=regime_grid_td)

# _ = sc.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=load, data_path=data_path,
#                  fig_path=fig_path, averagescheme=None, p_dimensionals=None, which_x='h_components',
#                  include_regimes=['chaotic'], save=True, fname='h_T_chaotic_all', labelsize=16, legend=True,
#                  xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography',
#                     title=r'fit to $C\alpha \delta_{rh} \Delta T_{rh}$, no time-averaging',
#                  c_peak='xkcd:forest green', c_rms='xkcd:periwinkle', fit=True, logx=True, logy=True, hscale=1,
#                  show_isoviscous=False, ylim=None, xlim=None, postprocess_kwargs=postprocess_kwargs,
#                  regime_grid=regime_grid_td)

# plot h scalings - with dT_m*delta*alpha

# _ = sc.subplots_topo_regimes(Ra_ls, eta_ls, regime_grid_td, regime_names_td, c_regimes=c_regimes_td, save=True, t1_grid=t1_grid,
#                          averagescheme='timefirst', legloc='upper right', which_x='h_components',
#                          load_grid=load, fig_path=fig_path, fname='h_T_timeavg', fig_fmt=fig_fmt, end_grid=end_grid,
#                          labelsize=14, xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography, $h^\prime$',
#                          xlabelpad=8, ylabelpad=20, fit=True, showallscatter=False,
#                          #xlim=(1e-8, 0.9e-6), ylim=(6e-3, 10e-2),
#                          logx=True, logy=True,
#                          regimes_title='Stationarity', leftleg_bbox=(-0.1, 0.95), data_path=data_path,
#                          postprocess_kwargs=postprocess_kwargs, #include_regimes=['chaotic']
#                               )
load = True

# heuristic with no averaging

# _ = sc.subplots_topo_regimes(Ra_ls, eta_ls, regime_grid_td, regime_names_td, c_regimes=c_regimes_td, save=True, t1_grid=t1_grid,
#                          which_x='h_components', averagescheme=None, legloc='upper right',
#                          load_grid=load,
#                          fig_path=fig_path, fname='h_T_all_scatter', fig_fmt=fig_fmt, end_grid=end_grid,
#                          labelsize=14, xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography $h^\prime$',
#                          xlabelpad=8, ylabelpad=20, fit=True, showallscatter=True,
#                          #xlim=(1e-8, 0.9e-6), ylim=(6e-3, 10e-2),
#                          logx=True, logy=True,
#                          regimes_title='Stationarity', leftleg_bbox=(-0.1, 0.95), data_path=data_path,
#                          postprocess_kwargs=postprocess_kwargs,)
load = True

# # plot h scalings with Ra
#
# _ = sc.subplots_topo_regimes(Ra_ls, eta_ls, regime_grid_td, regime_names_td, c_regimes=c_regimes_td, save=True, t1_grid=t1_grid,
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
