""" ASPECT runs: plot scalings of h vs. Ra or T-heuristic, with subplots by eta """

import sys
sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    c_regimes_td, fig_fmt, \
    regime_grid_td, regime_names_td, load_grid, alpha_m, p_Earth  # noqa: E402
from exotop import aspect_scalings as sc  # noqa: E402

load = load_grid

# plot h scalings - with dT_m*delta*alpha

sc.subplots_topo_regimes(Ra_ls, eta_ls, regime_grid_td, regime_names_td, c_regimes=c_regimes_td, save=True, t1=t1_grid,
                         T_components=True, averagefirst=True, legloc='upper right',
                         load=load, fig_path=fig_path, fname='h_T_all', fig_fmt=fig_fmt, end=end_grid,
                         labelsize=14, xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography, $h^\prime$',
                         xlabelpad=8, ylabelpad=-2, fit=True, alpha_m=alpha_m, showallscatter=False,
                         xlim=(1e-8, 0.9e-6), ylim=(6e-3, 10e-2), logx=True, logy=True,
                         regimes_title='Stationarity', leftleg_bbox=(-0.01, 0.95), data_path=data_path)
load = True

# same but just looking at each time point

sc.subplots_topo_regimes(Ra_ls, eta_ls, regime_grid_td, regime_names_td, c_regimes=c_regimes_td, save=True, t1=t1_grid,
                         T_components=True, averagefirst=False, legloc='upper right',
                         load=load,
                         fig_path=fig_path, fname='h_T_all_scatter', fig_fmt=fig_fmt, end=end_grid,
                         labelsize=14, xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography $h^\prime$',
                         xlabelpad=8, ylabelpad=-2, fit=True, alpha_m=alpha_m, showallscatter=True,
                         xlim=(1e-8, 0.9e-6), ylim=(6e-3, 10e-2), logx=True, logy=True,
                         regimes_title='Stationarity', leftleg_bbox=(-0.01, 0.95), data_path=data_path)
load = True

# plot h scalings with Ra

sc.subplots_topo_regimes(Ra_ls, eta_ls, regime_grid_td, regime_names_td, c_regimes=c_regimes_td, save=True, t1=t1_grid,
                         load=load, show_isoviscous=True, averagefirst=True, alpha_m=alpha_m,
                         fig_path=fig_path, fname='h_Ra_all', fig_fmt=fig_fmt, end=end_grid, labelsize=14,
                         xlabel='Ra_i', Ra_i=True,
                         ylabel='dynamic topography $h^\prime$', y2label='dynamic topography $h$ (km)',
                         xlabelpad=8, ylabelpad=5, fit=True, showallscatter=False,
                         xlim=(0.6e6, 5e8), #ylim=(1, 12),
                         logx=True, logy=True,
                         p_dimensionals=p_Earth,
                         regimes_title='Stationarity', leftleg_bbox=(-0.01, 0.95), data_path=data_path)
