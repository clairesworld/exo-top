""" Call functions to plot relevant 1D thermal model results which rely on processing of ASPECT runs """

import sys
sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    c_regimes_td, fig_fmt, \
    regime_grid_td, regime_names_td, load_grid, alpha_m, p_Earth  # noqa: E402
from exotop import aspect_scalings as sc  # noqa: E402
from exotop.model_1D import the_results as plottop  # noqa: E402
benchmark_path = '../benchmarks/'

# how h varies across key input parameters

fig, axes = plottop.plot_h_relative_multi(age=4, alpha=1, wspace=0.15, legsize=23.5, ticksize=20, labelsize=28,
                                          lw=4,  # legtitle=r'\textbf{\textit{Model}}',
                                          labelpad=10, legendtop=True, tickwidth=2,
                                          plots_save=True, fname='h_aspect_scaling', fig_path=fig_path,
                                          models=['dyn_top_rms'],
                                          labels=[r'$\Delta h \sim \alpha \Delta T_{rh} \delta^\zeta$'],
                                          x_vars=['age', 'M_p', 'CMF', 'H0', 'Ea'],
                                          c=['#d88868', '#749af3'], ftype='png')

# scale ocean capacity with h_peak - fractal model would be more realistic than scaling spherical harmonic spectrum

fig, axes = plottop.plot_ocean_capacity_relative(n_stats=1000, legsize=23.5, ticksize=20, labelsize=28, wspace=0.15,
                                                 titlesize=32,
                                                 fig_path=fig_path, savefig=False,
                                                 defaults='Venusbaseline', title='Ocean volume to submerge land',
                                                 spectrum_fpath=benchmark_path+'wei_Venus/',
                                                 spectrum_fname='model_power_m2_b.csv',
                                                 c='#81f79f', alpha=1, lw=4, ymin=0.3, ymax=1.8, labelpad=10,
                                                 set_ylim=True,
                                                 fname='ocean-vol', ftype='png')
