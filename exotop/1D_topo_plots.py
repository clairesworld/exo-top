""" Functions to plot relevant 1D thermal model results which rely on processing of ASPECT runs """

import sys
sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    c_regimes_td, fig_fmt, \
    regime_grid_td, regime_names_td, load_grid, alpha_m, p_Earth  # noqa: E402
from exotop import aspect_scalings as sc  # noqa: E402
from exotop.model_1D import plot_top as pltop  # noqa: E402

