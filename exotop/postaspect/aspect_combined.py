""" Functions to plot relevant 1D thermal model results which rely on processing of ASPECT runs """

import sys
# sys.path.insert(0, '/home/cmg76/Works/exo-top/')
sys.path.insert(0, '/home/claire/Works/exo-top/')
from exotop.postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    c_regimes_td, fig_fmt, \
    regime_grid_td, regime_names_td, load_grid, alpha_m, p_Earth  # noqa: E402
from exotop.postaspect.aspect_scalings import plot_save  # noqa: E402
from exotop.model_1D.parameters import M_E  # noqa: E402
from exotop.model_1D.the_results import bulk_planets  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402


def plot_h_heuristic_variations(default='Earthbaseline', x_param='M_p', x_min=0.1*M_E, x_max=5*M_E, y_params=None,
                                xlabel=r'$M_p$ ($M_E$)', save=True, age=4.5, fname='heuristic_mass_dependence',
                                ylabels=[r'$\Delta h_{rms}$ (km)', r'$\Delta T_{rh}$ (K)', r'$\delta_u$ (km)',
                                         r'$\alpha \Delta T_{rh} \delta_u$ (km)'], legend=True,
                                fig=None, axes=None, c='xkcd:drab green', lw=2, labelsize=20, legsize=14, **kwargs):

    pl_mass = bulk_planets(n=20, name=x_param, mini=x_min, maxi=x_max, like=default, random=False,
                           initial_kwargs={'tf': age}, **kwargs)
    if y_params is None:
        y_params = ['dyn_top_rms', 'T_rh', 'TBL_u', 'heuristic_h']
    if axes is None:
        fig, axes = plt.subplots(len(y_params), 1, figsize=(4, 4 * len(y_params)), sharex=True)

    for iax, ax in enumerate(axes):

        x = [vars(pl)[x_param] for pl in pl_mass]
        y = [vars(pl)[y_params[iax]][-1] for pl in pl_mass]
        print('x', x)
        print('y', y)
        ax.plot(x, y, c=c, lw=lw, label='1D model')

        if iax + 1 == len(axes):
            ax.set_xlabel(xlabel, fontsize=labelsize)
        ax.set_ylabel(ylabels[iax], fontsize=labelsize)

    if legend:
        axes[0].legend(fontsize=legsize, frameon=False)
    if save:
        plot_save(fig, fname, **kwargs)

