""" Plot new 1D thermal model topography after processing of ASPECT runs """

import sys

sys.path.insert(0, '/home/cmg76/Works/exo-top/')
sys.path.insert(0, '/home/claire/Works/exo-top/')
from exotop.postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, c_rms, c_peak, \
    c_regimes_td, fig_fmt, \
    regime_grid_td, regime_names_td, load_grid, alpha_m, p_Earth  # noqa: E402
from exotop.model_1D.parameters import M_E  # noqa: E402
from exotop.postaspect import aspect_combined as combo  # noqa: E402

fig_path_home = '/home/claire/Works/exo-top/exotop/figs_scratch/'

combo.plot_h_heuristic_variations(default='Earthbaseline', x_param='M_p', x_min=0.1 * M_E, x_max=5 * M_E, y_params=None,
                                  xlabel=r'$M_p$ ($M_E$)', save=True, fig_path=fig_path_home, fig_fmt=fig_fmt,
                                  legend=True,
                                  fname='heuristic_mass_dependence', nondimensional=True,
                                  c='xkcd:drab green', lw=3, labelsize=30, ticksize=12,
                                  overplot_aspect_cases=['Ra3e8-eta1e8-wide-ascii', 'Ra3e8-eta1e7-wide-ascii',
                                                         'Ra3e8-eta1e6-wide-ascii', 'Ra1e8-eta1e8-wide-ascii',
                                                         'Ra1e8-eta1e7-wide', 'Ra1e8-eta1e6-wide'])
