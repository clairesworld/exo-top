""" ASPECT runs: plot scalings of h vs. Ra or T-heuristic, with subplots by eta """

import sys

sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    c_regimes_td, fig_fmt, \
    regime_grid_td, regime_names_td, load_grid, p_Earth, postprocess_kwargs  # noqa: E402
from exotop.postaspect import aspect_scalings as sc  # noqa: E402

load = True #load_grid

# plot time-evolution of list of keys for all cases in given regime

# plot h scalings - with dT_m*delta*alpha

_ = sc.subplots_topo_regimes(Ra_ls, eta_ls, regime_grid_td, regime_names_td, c_regimes=c_regimes_td, save=True, t1_grid=t1_grid,
                         T_components=True, averagescheme='timefirst', legloc='upper right',
                         load_grid=load, fig_path=fig_path, fname='h_T_timeavg', fig_fmt=fig_fmt, end_grid=end_grid,
                         labelsize=14, xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography, $h^\prime$',
                         xlabelpad=8, ylabelpad=20, fit=True, showallscatter=False,
                         #xlim=(1e-8, 0.9e-6), ylim=(6e-3, 10e-2),
                         logx=True, logy=True,
                         regimes_title='Stationarity', leftleg_bbox=(-0.1, 0.95), data_path=data_path,
                         postprocess_kwargs=postprocess_kwargs, include_regimes=['chaotic'])
load = True

# same but just looking at each time point

# _ = sc.subplots_topo_regimes(Ra_ls, eta_ls, regime_grid_td, regime_names_td, c_regimes=c_regimes_td, save=True, t1_grid=t1_grid,
#                          T_components=True, averagescheme=None, legloc='upper right',
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
#                          xlabel='Ra_i', Ra_i=True,
#                          ylabel='dynamic topography $h^\prime$', y2label='dynamic topography $h$ (km)',
#                          xlabelpad=8, ylabelpad=5, fit=True, showallscatter=False,
#                          xlim=(0.6e6, 5e8), #ylim=(1, 12),
#                          logx=True, logy=True,
#                          p_dimensionals=p_Earth,
#                          regimes_title='Stationarity', leftleg_bbox=(-0.01, 0.95), data_path=data_path,
#                          postprocess_kwargs=postprocess_kwargs,)
