""" ASPECT runs: plot scalings of various parameters vs. Ra """

import sys
sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, c_regimes_td, fig_fmt, \
    regime_grid_td, regime_names_td, load_grid, p_Earth, postprocess_kwargs  # noqa: E402
from exotop.postaspect import aspect_scalings as sc  # noqa: E402


load = True # load_grid

# compare to Nu^1/3

# sc.subplots_Ra_scaling(Ra_ls, eta_ls, t1_grid=t1_grid, end_grid=end_grid, keys=['Nu'], data_path=data_path,
#                        fig_path=fig_path, load_grid=load, Ra_i=True, compare_exponent=[1 / 3], save=True,
#                        fname='Nu-compare', labelsize=14, ylabels=['Nu'], psuffixes=['_Nu', '_T'],
#                        postprocess_functions=[sc.Nu_at_ts, sc.T_parameters_at_sol], xlim=(1e5, 2e8), ylim=[(None)],
#                        cmap='winter', compare_pub=sc.moresi95, compare_label='', fig_fmt=fig_fmt, fit=True,
#                        postprocess_kwargs=postprocess_kwargs)


# # compare scalings of other output parameters with Ra
#
# sc.subplots_Ra_scaling(Ra_ls, eta_ls, t1_grid=t1_grid, end_grid=end_grid, keys=['Nu', 'delta_0', 'T_i'],
#                        data_path=data_path, fig_path=fig_path, load_grid=load, Ra_i=True, save=True,
#                        fname='delta-Nu-Ti', labelsize=14, ylabels=['Nu', r'$\delta_0$', r'$T_i$'],
#                        psuffixes=['_T', '_Nu'], postprocess_functions=[sc.T_parameters_at_sol, sc.Nu_at_ts],
#                        xlim=(1e5, 2e8), ylim=[(None), (None), (0.8, 1)], cmap='winter', compare_pub=sc.moresi95,
#                        compare_label='Moresi & Solomatov 1995', fig_fmt=fig_fmt, fit=True,
#                        postprocess_kwargs=postprocess_kwargs)



# plot scalings of chaotic time-dependence T parameters - effective Ra_i

sc.subplots_Ra_scaling(Ra_ls[4:], eta_ls[1:], t1_grid=t1_grid[1:, 4:], end_grid=end_grid[1:, 4:],
                       keys=['delta_rh', 'dT_rh', 'Nu'], data_path=data_path, fig_path=fig_path,
                       load_grid=load,# load_grid[1:, 4:],
                       Ra_i='eff', save=True, fname='delta_rh-Nu-chaotic-eff',
                       labelsize=14, ylabels=[r'$\delta_{rh}$', r'$\Delta T_{rh}$', 'Nu'], psuffixes=['_T', '_Nu'],
                       title='Chaotic time-dependence', postprocess_functions=[sc.T_parameters_at_sol, sc.Nu_at_ts],
                       xlim=(0.7e6, 4e7), ylim=[(None), (None), (None)],
                       legloc=['lower left', 'upper left', 'upper left'], cmap='winter', compare_pub=sc.moresi95,
                       compare_label='', fig_fmt=fig_fmt, fit=True, postprocess_kwargs=postprocess_kwargs)

sc.subplots_Ra_scaling(Ra_ls[:3], eta_ls, t1_grid=t1_grid[:, :3], end_grid=end_grid[:, :3], keys=['delta_rh', 'dT_rh'],
                       data_path=data_path, fig_path=fig_path, load_grid=load, Ra_i='eff', save=True,
                       fname='delta_rh-steady-eff', labelsize=14, ylabels=[r'$\delta_{rh}$', r'$\Delta T_{rh}$'],
                       psuffixes=['_T'], title='Steady-state', postprocess_functions=[sc.T_parameters_at_sol],
                       xlim=(0.3e5, 1e6), ylim=[(None), (None)], legloc=['lower left', 'upper left'], cmap='winter',
                       compare_pub=sc.moresi95, compare_label='', fig_fmt=fig_fmt, fit=True,
                       postprocess_kwargs=postprocess_kwargs)



# plot scalings of chaotic time-dependence T parameters - uncorrected Ra_i

sc.subplots_Ra_scaling(Ra_ls[4:], eta_ls[1:], t1_grid=t1_grid[1:, 4:], end_grid=end_grid[1:, 4:],
                       keys=['delta_0', 'delta_rh', 'dT_rh', 'Nu'], data_path=data_path, fig_path=fig_path, load_grid=load,
                       Ra_i=True, save=True, fname='delta_rh-Nu-chaotic', compare_label='Moresi & Solomatov 1995', labelsize=14,
                       ylabels=[r'$\delta_{0}$', r'$\delta_{rh}$', r'$\Delta T_{rh}$', 'Nu'], psuffixes=['_T', '_Nu'],
                       title='Chaotic time-dependence', postprocess_functions=[sc.T_parameters_at_sol, sc.Nu_at_ts],
                       xlim=(1e7, 2e8), ylim=[(None), (None), (None)],
                       legloc=['lower left', 'upper left', 'lower left'], cmap='winter', compare_pub=sc.moresi95,
                       fig_fmt=fig_fmt, fit=True, postprocess_kwargs=postprocess_kwargs)

# sc.subplots_Ra_scaling(Ra_ls[:3], eta_ls, t1_grid=t1_grid[:, :3], end_grid=end_grid[:, :3], keys=['delta_rh', 'dT_rh'],
#                        data_path=data_path, fig_path=fig_path, load_grid=load, Ra_i=True, save=True,
#                        fname='delta_rh-steady', labelsize=14, ylabels=[r'$\delta_{rh}$', r'$\Delta T_{rh}$'],
#                        psuffixes=['_T'], title='Steady-state', postprocess_functions=[sc.T_parameters_at_sol],
#                        xlim=(0.1e6, 0.5e7), ylim=[(None), (None)], legloc=['lower left', 'upper left'], cmap='winter',
#                        compare_pub=sc.moresi95, compare_label='', fig_fmt=fig_fmt, fit=True,
#                        postprocess_kwargs=postprocess_kwargs)
