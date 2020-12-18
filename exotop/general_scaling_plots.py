import sys

sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, c_regimes_td, fig_fmt, \
    regime_grid_td, regime_names_td, load_grid, alpha_m, p_Earth
from exotop import aspect_scalings as sc

load = True

# compare to Nu^1/3

sc.subplots_Ra_scaling(Ra_ls, eta_ls, t1=t1_grid, end=end_grid, keys=['Nu'], data_path=data_path,
                       fig_path=fig_path, load=load, save=True, fname='Nu-compare', xlim=(1e5, 2e8),
                       ylim=[(None)], labelsize=14,
                       ylabels=['Nu'], psuffixes=['_Nu'],
                       postprocess_functions=[sc.Nu_at_ts], Ra_i=True,
                       compare_label='', compare_pub=sc.moresi95, compare_exponent=[1/3],
                       # compare_data=[([10**10.0102, 10**11.9857, 10**14.0175], [2.54284, 2.10300, 1.73935]), # d_eta 1e5
                       #               (), (), ()],
                       fig_fmt=fig_fmt, cmap='winter', fit=True)


# # compare scalings of other output parameters with Ra

sc.subplots_Ra_scaling(Ra_ls, eta_ls, t1=t1_grid, end=end_grid, keys=['Nu', 'delta_0', 'T_i'], data_path=data_path,
                       fig_path=fig_path, load=load, save=True, fname='delta-Nu-Ti', xlim=(1e5, 2e8),
                       ylim=[(None),(None),(0.8,1)], labelsize=14,
                       ylabels=['Nu', r'$\delta_0$', r'$T_i$'], psuffixes=['_T', '_Nu'],
                       postprocess_functions=[sc.T_parameters_at_sol, sc.Nu_at_ts], Ra_i=True,
                       compare_label='Moresi & Solomatov 1995', compare_pub=sc.moresi95,
                       fig_fmt=fig_fmt, cmap='winter', fit=True)



# # plot scalings of chaotic time-dependence T parameters - effective Ra_i

sc.subplots_Ra_scaling(Ra_ls[4:], eta_ls[1:], t1=t1_grid[1:,4:], end=end_grid[1:,4:], keys=['delta_rh', 'dT_rh', 'Nu'], data_path=data_path,
                       fig_path=fig_path, load=load_grid[1:,4:], save=True, fname='delta_rh-Nu-chaotic-eff', xlim=(0.7e6, 4e7),
                       ylim=[(None) , (None), (None)], labelsize=14, title='Chaotic time-dependence',
                       ylabels=[r'$\delta_{rh}$', r'$\Delta T_{rh}$', 'Nu'], psuffixes=['_T', '_Nu'], legloc=['lower left', 'upper left', 'upper left'],
                       postprocess_functions=[sc.T_parameters_at_sol, sc.Nu_at_ts], Ra_i='eff',
                       compare_label='', compare_pub=sc.moresi95,
                       fig_fmt=fig_fmt, cmap='winter', fit=True)

sc.subplots_Ra_scaling(Ra_ls[:3], eta_ls, t1=t1_grid[:,:3], end=end_grid[:,:3], keys=['delta_rh', 'dT_rh'], data_path=data_path,
                       fig_path=fig_path, load=load, save=True, fname='delta_rh-steady-eff', xlim=(0.3e5, 1e6),
                       ylim=[(None) , (None)], labelsize=14, title='Steady-state', legloc=['lower left','upper left'],
                       ylabels=[r'$\delta_{rh}$', r'$\Delta T_{rh}$'], psuffixes=['_T'],
                       postprocess_functions=[sc.T_parameters_at_sol], Ra_i='eff',
                      compare_label='', compare_pub=sc.moresi95,
                       fig_fmt=fig_fmt, cmap='winter', fit=True)



# # plot scalings of chaotic time-dependence T parameters - uncorrected Ra_i
#
sc.subplots_Ra_scaling(Ra_ls[4:], eta_ls[1:], t1=t1_grid[1:,4:], end=end_grid[1:,4:], keys=['delta_rh', 'dT_rh', 'Nu'], data_path=data_path,
                       fig_path=fig_path, load=load, save=True, fname='delta_rh-Nu-chaotic', xlim=(1e7, 2e8),
                       ylim=[(None) , (None), (None)], labelsize=14, title='Chaotic time-dependence',
                       ylabels=[r'$\delta_{rh}$', r'$\Delta T_{rh}$', 'Nu'], psuffixes=['_T', '_Nu'],
                       postprocess_functions=[sc.T_parameters_at_sol, sc.Nu_at_ts], Ra_i=True, legloc=['lower left', 'upper left', 'lower left'],
                      compare_label='', compare_pub=sc.moresi95,
                       fig_fmt=fig_fmt, cmap='winter', fit=True)

sc.subplots_Ra_scaling(Ra_ls[:3], eta_ls, t1=t1_grid[:,:3], end=end_grid[:,:3], keys=['delta_rh', 'dT_rh'], data_path=data_path,
                       fig_path=fig_path, load=load, save=True, fname='delta_rh-steady', xlim=(0.1e6, 0.5e7),
                       ylim=[(None) , (None)], labelsize=14, title='Steady-state',
                       ylabels=[r'$\delta_{rh}$', r'$\Delta T_{rh}$'], psuffixes=['_T'],
                       postprocess_functions=[sc.T_parameters_at_sol], Ra_i=True, legloc=['lower left', 'upper left'],
                      compare_label='', compare_pub=sc.moresi95,
                       fig_fmt=fig_fmt, cmap='winter', fit=True)
