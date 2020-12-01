import sys

sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, c_regimes_td, fig_fmt, \
    regime_grid_td, regime_names_td, load_grid, alpha_m
from exotop import aspect_scalings as sc


load = True # load_grid

# ### plot h scalings with Ra

sc.subplots_topo_regimes(Ra_ls, eta_ls, regime_grid_td, regime_names_td, c_regimes=c_regimes_td, save=True, t1=t1_grid, load=load,
                         fig_path=fig_path, fname='h_Ra_all', fig_fmt=fig_fmt, end=end_grid, labelsize=14, xlabel='Ra',
                         ylabel='dynamic topography (km)', xlabelpad=8, ylabelpad=-2, fit=True,
                         xlim=(0.6e6, 5e8), ylim=(1, 12), logx=True, logy=True, hscale=2e-5 * 2700 * 2890,
                         regimes_title='Time-dependence', leftleg_bbox=(-0.01, 0.95), data_path=data_path)
load = True

## plot h scalings - with dT_m*delta*alpha

sc.subplots_topo_regimes(Ra_ls, eta_ls, regime_grid_td, regime_names_td, c_regimes=c_regimes_td, save=True, t1=t1_grid,
                         T_components=True,
                         load=load, fig_path=fig_path, fname='h_T_all', fig_fmt=fig_fmt, end=end_grid,
                         labelsize=14, xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography',
                         xlabelpad=8, ylabelpad=-2, fit=True, alpha_m=alpha_m, showallscatter=True,
                         xlim=(1e-8, 0.9e-6), ylim=(6e-3, 10e-2), logx=True, logy=True,
                         regimes_title='Time-dependence',leftleg_bbox=(-0.01, 0.95), data_path=data_path)

## plot scalings of other output parameters with Ra

# sc.subplots_Ra_scaling(Ra_ls, eta_ls, t1=t1_grid, end=end_grid, keys=['Nu', 'delta_0', 'T_i'], data_path=data_path,
#                        fig_path=fig_path, load=load, save=True, fname='delta-Nu-Ti', xlim=(1e5, 5e8),
#                        ylim=[(None),(None),(0.8,1)], labelsize=14,
#                        ylabels=['Nu', r'$\delta_0$', r'$T_i$'], psuffixes=['_T', '_Nu'],
#                        postprocess_functions=[sc.T_parameters_at_sol, sc.Nu_at_ts], Ra_i=True,
#                        compare_label='Moresi & Solomatov 1995', compare_pub=sc.moresi95,
#                        fig_fmt=fig_fmt, cmap='winter', fit=True)

