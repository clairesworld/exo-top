import sys

sys.path.insert(0, '/home/cmg76/Works/exo-top/')
# import numpy as np
from exotop import aspect_scalings as sc
from exotop.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, c_regimes, \
    fig_fmt, \
    regime_grid, regime_names, load_grid, alpha_m

# plot convection regimes in parameter space
# example of transitional is Ra3e7 eta1e6 - still has regular cells

# sc.plot_parameter_grid(Ra_ls, eta_ls, function=sc.regime_to_digital, regime_grid=regime_grid, regime_names=regime_names,
#                        data_path=data_path, fig_path=fig_path, load=True, vmin=1, vmax=3,
#                        save=True, fname='grid', labelsize=16, fig_fmt=fig_fmt, t1=t1_grid, end=end_grid,
#                        cticklabels=['steady', 'transitional', 'chaotic'], cticks=[1.5, 2, 2.5],
#                        overplot_h=False, nlevels_contour=14, cmap='jet', clist=c_regimes, cmap_contours='autumn')

# lid mobility

sc.plot_parameter_grid(Ra_ls, eta_ls, function=sc.sfc_mobility_at_sol, regime_grid=regime_grid,
                       data_path=data_path, fig_path=fig_path, load=load_grid, vmax=1, log=True, clabel='log surface mobility',
                       save=True, fname='surface-mobility', labelsize=16, fig_fmt=fig_fmt, t1=t1_grid, end=end_grid,
                       overplot_h=False, nlevels_contour=14, cmap='Wistia', set_over='xkcd:navy blue')
