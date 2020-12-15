import sys

sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, c_regimes_td, fig_fmt, \
    regime_grid_td, regime_names_td, load_grid, alpha_m, p_Earth
from exotop import aspect_scalings as sc


load = load_grid

# heuristic scalings master

sc.plot_heuristic_scalings(Ra_ls, eta_ls, regime_grid=regime_grid_td, t1=t1_grid, load=True, end=end_grid, literature_file=None, legend=True,
                            c='k', ylim=None, xlim=None, which_h='rms', data_path=data_path,
                            save=True, fname='model-data', ylim=[6e-3, 4e-2])


#
# # plot h scalings - with dT_m*delta*alpha
# #
# sc.subplots_topo_regimes(Ra_ls, eta_ls, regime_grid_td, regime_names_td, c_regimes=c_regimes_td, save=True, t1=t1_grid,
#                          T_components=True, averagefirst=True, legloc='upper right',
#                          load=load, fig_path=fig_path, fname='h_T_all', fig_fmt=fig_fmt, end=end_grid,
#                          labelsize=14, xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography, $h^\prime$',
#                          xlabelpad=8, ylabelpad=-2, fit=True, alpha_m=alpha_m, showallscatter=False,
#                          xlim=(1e-8, 0.9e-6), ylim=(6e-3, 10e-2), logx=True, logy=True,
#                          regimes_title='Stationarity', leftleg_bbox=(-0.01, 0.95), data_path=data_path)
# load = True
#
# # # same but just looking at each time point
# sc.subplots_topo_regimes(Ra_ls, eta_ls, regime_grid_td, regime_names_td, c_regimes=c_regimes_td, save=True, t1=t1_grid,
#                          T_components=True, averagefirst=False, legloc='upper right',
#                          load=load,
#                          fig_path=fig_path, fname='h_T_all_scatter', fig_fmt=fig_fmt, end=end_grid,
#                          labelsize=14, xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography $h^\prime$',
#                          xlabelpad=8, ylabelpad=-2, fit=True, alpha_m=alpha_m, showallscatter=True,
#                          xlim=(1e-8, 0.9e-6), ylim=(6e-3, 10e-2), logx=True, logy=True,
#                          regimes_title='Stationarity', leftleg_bbox=(-0.01, 0.95), data_path=data_path)
# load = True
# #
# # # plot h scalings with Ra
# sc.subplots_topo_regimes(Ra_ls, eta_ls, regime_grid_td, regime_names_td, c_regimes=c_regimes_td, save=True, t1=t1_grid,
#                          load=load, show_isoviscous=True, averagefirst=True, alpha_m=alpha_m,
#                          fig_path=fig_path, fname='h_Ra_all', fig_fmt=fig_fmt, end=end_grid, labelsize=14,
#                          xlabel='Ra_i', Ra_i=True,
#                          ylabel='dynamic topography $h^\prime$', y2label='dynamic topography $h$ (km)',
#                          xlabelpad=8, ylabelpad=5, fit=True, showallscatter=False,
#                          xlim=(0.6e6, 5e8), #ylim=(1, 12),
#                          logx=True, logy=True,
#                          p_dimensionals=p_Earth,
#                          regimes_title='Stationarity', leftleg_bbox=(-0.01, 0.95), data_path=data_path)

#
# load = 'auto'
#
# # compare scalings of other output parameters with Ra
#
# sc.subplots_Ra_scaling(Ra_ls, eta_ls, t1=t1_grid, end=end_grid, keys=['Nu', 'delta_0', 'T_i'], data_path=data_path,
#                        fig_path=fig_path, load=load, save=True, fname='delta-Nu-Ti', xlim=(1e5, 2e8),
#                        ylim=[(None),(None),(0.8,1)], labelsize=14,
#                        ylabels=['Nu', r'$\delta_0$', r'$T_i$'], psuffixes=['_T', '_Nu'],
#                        postprocess_functions=[sc.T_parameters_at_sol, sc.Nu_at_ts], Ra_i=True,
#                        compare_label='Moresi & Solomatov 1995', compare_pub=sc.moresi95,
#                        fig_fmt=fig_fmt, cmap='winter', fit=True)
#
# # plot scalings of chaotic time-dependence T parameters - effective Ra_i
#
# sc.subplots_Ra_scaling(Ra_ls[4:], eta_ls[1:], t1=t1_grid[1:,4:], end=end_grid[1:,4:], keys=['delta_rh', 'dT_rh'], data_path=data_path,
#                        fig_path=fig_path, load=load, save=True, fname='delta_rh-chaotic-eff', xlim=(0.7e6, 4e7),
#                        ylim=[(None) , (None)], labelsize=14, title='Chaotic time-dependence',
#                        ylabels=[r'$\delta_{rh}$', r'$\Delta T_{rh}$'], psuffixes=['_T'], legloc=['lower left', 'upper left'],
#                        postprocess_functions=[sc.T_parameters_at_sol], Ra_i='eff',
#                       compare_label='', compare_pub=sc.moresi95,
#                        fig_fmt=fig_fmt, cmap='winter', fit=True)
#
# sc.subplots_Ra_scaling(Ra_ls[:3], eta_ls, t1=t1_grid[:,:3], end=end_grid[:,:3], keys=['delta_rh', 'dT_rh'], data_path=data_path,
#                        fig_path=fig_path, load=load, save=True, fname='delta_rh-steady-eff', xlim=(0.3e5, 1e6),
#                        ylim=[(None) , (None)], labelsize=14, title='Steady-state', legloc=['lower left','upper left'],
#                        ylabels=[r'$\delta_{rh}$', r'$\Delta T_{rh}$'], psuffixes=['_T'],
#                        postprocess_functions=[sc.T_parameters_at_sol], Ra_i='eff',
#                       compare_label='', compare_pub=sc.moresi95,
#                        fig_fmt=fig_fmt, cmap='winter', fit=True)
#
#
#
# # plot scalings of chaotic time-dependence T parameters - uncorrected Ra_i
#
# sc.subplots_Ra_scaling(Ra_ls[4:], eta_ls[1:], t1=t1_grid[1:,4:], end=end_grid[1:,4:], keys=['delta_rh', 'dT_rh'], data_path=data_path,
#                        fig_path=fig_path, load=load, save=True, fname='delta_rh-chaotic', xlim=(1e7, 2e8),
#                        ylim=[(None) , (None)], labelsize=14, title='Chaotic time-dependence',
#                        ylabels=[r'$\delta_{rh}$', r'$\Delta T_{rh}$'], psuffixes=['_T'],
#                        postprocess_functions=[sc.T_parameters_at_sol], Ra_i=True, legloc=['lower left', 'upper left'],
#                       compare_label='', compare_pub=sc.moresi95,
#                        fig_fmt=fig_fmt, cmap='winter', fit=True)
#
# sc.subplots_Ra_scaling(Ra_ls[:3], eta_ls, t1=t1_grid[:,:3], end=end_grid[:,:3], keys=['delta_rh', 'dT_rh'], data_path=data_path,
#                        fig_path=fig_path, load=load, save=True, fname='delta_rh-steady', xlim=(0.1e6, 0.5e7),
#                        ylim=[(None) , (None)], labelsize=14, title='Steady-state',
#                        ylabels=[r'$\delta_{rh}$', r'$\Delta T_{rh}$'], psuffixes=['_T'],
#                        postprocess_functions=[sc.T_parameters_at_sol], Ra_i=True, legloc=['lower left', 'upper left'],
#                       compare_label='', compare_pub=sc.moresi95,
#                        fig_fmt=fig_fmt, cmap='winter', fit=True)
