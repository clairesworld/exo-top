import sys

sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, c_regimes_td, \
    fig_fmt, regime_grid_td, regime_names_td, load_grid, alpha_m
from exotop import aspect_scalings as sc


load = True  # load_grid

# plot convection regimes in parameter space
# example of transitional is Ra3e7 eta1e6 - still has regular cells
sc.plot_parameter_grid(Ra_ls, eta_ls, function=sc.regime_to_digital, regime_grid=regime_grid_td, regime_names=regime_names_td,
                       data_path=data_path, fig_path=fig_path, load=True, vmin=1, vmax=4, discrete=True, log=False,
                       title='Styles of time-dependence',
                       save=True, fname='time-dependence', labelsize=16, fig_fmt=fig_fmt, t1=t1_grid, end=end_grid,
                       cticklabels=['steady', 'transitional', 'chaotic', 'not stagnant lid'], cticks=[1.5, 2, 3, 3.5],
                       overplot_h=False, nlevels_contour=14, cmap='jet', clist=c_regimes_td, cmap_contours='autumn',
                       )

# lid mobility
# sc.plot_parameter_grid(Ra_ls, eta_ls, function=sc.surf_mobility_at_sol, regime_grid=regime_grid_td,
#                        data_path=data_path, fig_path=fig_path, load=True, vmax=1, log=True, clabel='log surface mobility',
#                        title='Stagnant lid regime',
#                        save=True, fname='surface-mobility', labelsize=16, fig_fmt=fig_fmt, t1=t1_grid, end=end_grid,
#                        overplot_h=False, nlevels_contour=14, discrete=False, cmap='Wistia', set_over='xkcd:navy blue')
