""" ASPECT runs: get heuristic scalings and plot model vs. data (fit error) """

import sys

sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    c_regimes_td, fig_fmt, postprocess_kwargs, \
    regime_grid_td, regime_names_td, load_grid, p_Earth  # noqa: E402
from exotop.postaspect import aspect_scalings as sc  # noqa: E402

load = True #load_grid

sc.plot_fit_parameter_grid(Ra_ls, eta_ls,  data_path=data_path, fig_path=fig_path, load_grid=load,
                       averagescheme='timefirst', which_x=('Ra', 'eta'), regime_grid=regime_grid_td, include_regimes=['chaotic'],
                        save=True, fname='fit-grid', labelsize=16, fig_fmt=fig_fmt, t1_grid=t1_grid, end_grid=end_grid,
                        # cticklabels=None, cticks=None,
                        title='', lognorm=False, log=False, clabel=r'$h_{rms}$', postprocess_kwargs=postprocess_kwargs,
                         nlevels_contour=10, cmap='Greys_r', cmap_contours='spring')

const, expon = sc.plot_model_data(Ra_ls, eta_ls, regime_grid=regime_grid_td, t1_grid=t1_grid, load_grid=load,
                                  end_grid=end_grid, literature_file=None, legend=True, cmap='winter',
                                  postprocess_kwargs=postprocess_kwargs, c='k', averagescheme='timefirst',
                                  ylim=[4e-3, 4e-2], which_x=('Ra_i_eff', 'eta'), which_h='rms', data_path=data_path, save=True,
                                  fname='fit-power-Raeff-chaotic_timeavg', cbar='eta', include_regimes=['chaotic'],
                                  intercept=False, fig_fmt=fig_fmt)

const, expon = sc.plot_model_data(Ra_ls, eta_ls, regime_grid=regime_grid_td, t1_grid=t1_grid, load_grid=load,
                                  end_grid=end_grid, literature_file=None, legend=True, cmap='winter',
                                  postprocess_kwargs=postprocess_kwargs, c='k', averagescheme='timefirst',
                                  ylim=[4e-3, 4e-2], which_x=('Ra_i', 'eta'), which_h='rms', data_path=data_path, save=True,
                                  fname='fit-power-Rai-chaotic_timeavg', cbar='eta', include_regimes=['chaotic'],
                                  intercept=False, fig_fmt=fig_fmt)

const, expon = sc.plot_model_data(Ra_ls, eta_ls, regime_grid=regime_grid_td, t1_grid=t1_grid, load_grid=load,
                                  end_grid=end_grid, literature_file=None, legend=True, cmap='winter',
                                  postprocess_kwargs=postprocess_kwargs, c='k', averagescheme='timefirst',
                                  ylim=[4e-3, 4e-2], which_x=('Ra', 'eta'), which_h='rms', data_path=data_path, save=True,
                                  fname='fit-power-Ra-chaotic_timeavg', cbar='eta', include_regimes=['chaotic'],
                                  intercept=False, fig_fmt=fig_fmt)
print('fit parameters:', const, expon)

# heuristic scalings master

# const, expon = sc.plot_model_data(Ra_ls, eta_ls, regime_grid=regime_grid_td, t1_grid=t1_grid, load_grid=load,
#                                   end_grid=end_grid, literature_file=None, legend=True,
#                                   postprocess_kwargs=postprocess_kwargs, c='k', averagescheme='timefirst',
#                                   ylim=[4e-3, 4e-2], which_h='rms', data_path=data_path, save=True,
#                                   fname='model-data-chaotic_timeavg', cbar='eta', include_regimes=['chaotic'],
#                                   intercept=True, fig_fmt=fig_fmt)
# print('fit parameters:', const, expon)

# sc.plot_heuristic_scalings(Ra_ls, eta_ls, regime_grid=regime_grid_td, t1_grid=t1_grid, load_grid=load,
#                            end_grid=end_grid, literature_file=None, legend=True, c='k', ylim=[2e-3, 4e-2],
#                            which_h='rms', data_path=data_path, save=True, fname='model-data-regimes',
#                            regime_names=regime_names_td, clist=c_regimes_td, cbar='regime', intercept=True,
#                            postprocess_kwargs=postprocess_kwargs, fig_fmt=fig_fmt)
