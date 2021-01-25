import sys
import os
sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    fig_fmt, regime_grid_td, postprocess_kwargs, regime_names_td, \
    load_grid    # noqa: E402
from exotop.postaspect import aspect_scalings as sc  # noqa: E402
from exotop.postaspect import aspect_postprocessing2 as ap
from exotop.useful_and_bespoke import hide_log_ticklabels, not_iterable, dark_background
import numpy as np
import matplotlib.pyplot as plt

ticksize = 20
axissize = 30
fig, ax = sc.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=True, data_path=data_path,
                 fig_path=fig_path, averagescheme='timefirst', p_dimensionals=None, which_x='Ra_i_eff',
                 sigma=1,
                 include_regimes=['chaotic'], save=False, fname='h_Raieff_chaotic_timeavg', labelsize=axissize,
                       legend=False,
                 xlabel=r'Ra$_{i,eff}$', ylabel='dynamic topography',
                 title=r'fit to CRa$_{i,eff}^n$', fiterror=False,
                 c_peak='k', c_rms='xkcd:lilac', fit=True, logx=True, logy=True, hscale=1,
                 show_isoviscous=False, ylim=None, xlim=None, postprocess_kwargs=postprocess_kwargs,
                 regime_grid=regime_grid_td)

ax.tick_params(axis='x', labelsize=ticksize, pad=15)
ax.tick_params(axis='y', labelsize=ticksize, pad=15)
fig, ax = dark_background(fig, ax)
sc.plot_save(fig, fname='h_Ra_chaotic', fig_path=fig_path+'slides/', facecolor=fig.get_facecolor())