import sys

sys.path.insert(0, '/home/cmg76/Works/exo-top/')
# import numpy as np
from exotop import aspect_scalings as sc
from exotop.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, c_regimes, fig_fmt, \
    regime_grid, regime_names, load_grid, alpha_m

load = 'auto'  # True

# ### plot h scalings with Ra
#
# sc.subplots_topo_regimes(Ra_ls, eta_ls, regime_grid, regime_names, c_regimes=c_regimes, save=True, t1=t1_grid, load=load,
#                          fig_path=fig_path, fname='h_Ra_all', fig_fmt=fig_fmt, end=end_grid, labelsize=14, xlabel='Ra',
#                          ylabel='dynamic topography (km)', xlabelpad=8, ylabelpad=10, fit=True, xlim=(0.7e6, 3.05e8),
#                          ylim=(1, 10), logx=True, logy=True, hscale=2e-5 * 2700 * 2890, data_path=data_path)

### plot h scalings - with dT_m*delta*alpha
#
# sc.subplots_topo_regimes(Ra_ls, eta_ls, regime_grid, regime_names, c_regimes=c_regimes, save=True, t1=t1_grid, T_components=True,
#                          load=load, fig_path=fig_path, fname='h_T_all', fig_fmt=fig_fmt, end=end_grid,
#                          labelsize=14, xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography',
#                          xlabelpad=8, ylabelpad=14, fit=True, alpha_m=alpha_m, showallscatter=True, xlim=(3e-8, 7e-7),
#                          logx=True, logy=True, data_path=data_path)

### plot scalings of other output parameters with Ra

sc.subplots_Ra_scaling(Ra_ls, eta_ls, t1=t1_grid, end=end_grid, keys=['Nu', 'delta_0', 'T_i'], data_path=data_path,
                       fig_path=fig_path, load=load, save=True, fname='delta-Nu-Ti',
                       ylabels=['Nu', r'$\delta$', r'$T_i$'], psuffixes=['_T', '_Nu'],
                       postprocess_functions=[sc.T_parameters_at_sol, sc.Nu_at_ts], compare_pub=None, fig_fmt=fig_fmt,
                       fit=True)

