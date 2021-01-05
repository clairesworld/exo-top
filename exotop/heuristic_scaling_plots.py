""" ASPECT runs: get heuristic scalings and plot model vs. data (fit error) """

import sys
sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    c_regimes_td, fig_fmt, \
    regime_grid_td, regime_names_td, load_grid, alpha_m, p_Earth  # noqa: E402
from exotop import aspect_scalings as sc  # noqa: E402


load = True #load_grid

# heuristic scalings master

const, expon = sc.plot_heuristic_scalings(Ra_ls, eta_ls, regime_grid=regime_grid_td, t1=t1_grid, load=load, end=end_grid,
                           literature_file=None, legend=True, cbar='eta', alpha_m=alpha_m, intercept=True,
                            c='k', which_h='rms', data_path=data_path, include_regimes=['chaotic'],
                            save=True, fname='model-data-chaotic', ylim=[6e-3, 4e-2])
print('fit parameters:', const, expon)

sc.plot_heuristic_scalings(Ra_ls, eta_ls, regime_grid=regime_grid_td, t1=t1_grid, load=load, end=end_grid,
                           literature_file=None, legend=True, cbar='regime', clist=c_regimes_td, intercept=True,
                           regime_names=regime_names_td, alpha_m=alpha_m,
                            c='k', which_h='rms', data_path=data_path, #include_regimes=['chaotic'],
                            save=True, fname='model-data-regimes', ylim=[6e-3, 4e-2])