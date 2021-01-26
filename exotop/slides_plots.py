import sys
sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    fig_fmt, regime_grid_td, postprocess_kwargs, regime_names_td, load_grid, p_Earth    # noqa: E402
from exotop.postaspect import aspect_scalings as sc  # noqa: E402
from exotop.useful_and_bespoke import dark_background
import numpy as np
import matplotlib.lines as mlines
import matplotlib.ticker as ticker

ticksize = 22
axissize = 40
c_fit = 'xkcd:off white'
c_rms = ['xkcd:lime green', 'xkcd:lilac', 'xkcd:orange', 'xkcd:yellow']
lw = 5

regimes = ['chaotic']
for regime in regimes:
    if regime == 'all':
        include_regimes = ['steady', 'trans.', 'chaotic']
        ylim = [6e-3, 4e-2]
        xlim = [0.7e5, 3e7]
        yticks = [6e-3, 1e-2, 4e-2]
        xticks = [1e5, 1e6, 1e7]
        fitlabel = r'$\Delta h = 0.345$ Ra$^{-0.212}$'
        handles = [mlines.Line2D([], [], color=c_fit, marker='*', ls='--',
                                     markersize=0, lw=lw, label=fitlabel),
        mlines.Line2D([], [], color=c_rms[0], marker='o', ls='--',
                                     markersize=20, lw=0, label=r'$\Delta \eta = 10^{5}$'),
        mlines.Line2D([], [], color=c_rms[1], marker='o', ls='--',
                                     markersize=20, lw=0, label=r'$\Delta \eta = 10^{6}$'),
        mlines.Line2D([], [], color=c_rms[2], marker='o', ls='--',
                                     markersize=20, lw=0, label=r'$\Delta \eta = 10^{7}$'),
        mlines.Line2D([], [], color=c_rms[3], marker='o', ls='--',
                                     markersize=20, lw=0, label=r'$\Delta \eta = 10^{8}$')]
    elif regime == 'chaotic':
        include_regimes = ['chaotic']
        ylim = [6e-3, 1.4e-2]
        xlim = [1.999e6, 3.001e7]
        yticks = [6e-3, 7e-3, 8e-3, 9e-3, 1e-2]
        xticks = [1.999e6, 2e6, 1e7, 3e7, 3.001e7]
        fitlabel = r'$\Delta h = 0.094$ Ra$^{-0.151}$'
        handles = [mlines.Line2D([], [], color=c_fit, marker='*', ls='--',
                                 markersize=0, lw=lw, label=fitlabel),
                   mlines.Line2D([], [], color=c_rms[1], marker='o', ls='--',
                                 markersize=20, lw=0, label=r'$\Delta \eta = 10^{6}$'),
                   mlines.Line2D([], [], color=c_rms[2], marker='o', ls='--',
                                 markersize=20, lw=0, label=r'$\Delta \eta = 10^{7}$'),
                   mlines.Line2D([], [], color=c_rms[3], marker='o', ls='--',
                                 markersize=20, lw=0, label=r'$\Delta \eta = 10^{8}$')]

    fig, ax = sc.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=True, data_path=data_path,
                     fig_path=fig_path, averagescheme='timefirst', which_x='Ra_i_eff',
                     sigma=1,
                     include_regimes=include_regimes, save=False, fname='h_Raieff_chaotic_timeavg', labelsize=axissize,
                           legend=False, figsize=(16,9), showpeak=False, lw=lw, ms=25, elw=2, ecapsize=8,
                     xlabel=r'Ra$_{i,eff}$', ylabel=r'dynamic topography $\Delta h^\prime$', ylabelpad=20, xlabelpad=13,
                     title='', fiterror=False, c_fit=c_fit,
                     c_peak='k', c_rms=c_rms,
                           fit=True, logx=True, logy=True, hscale=1,
                     show_isoviscous=False, ylim=ylim, xlim=xlim, postprocess_kwargs=postprocess_kwargs,
                     regime_grid=regime_grid_td, p_dimensionals=None)

    ax.tick_params(axis='x', labelsize=ticksize, pad=15)
    ax.tick_params(axis='y', labelsize=ticksize, pad=15)
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
    ax.ticklabel_format(style='plain', axis='y', useOffset=False)
    # if yticks is not None:
    #     ax.set_yticks(yticks)
    if xticks is not None:
        ax.set_xticks(xticks)

    ax.legend(handles=handles, frameon=False, fontsize=25, ncol=1, bbox_to_anchor=(1.01, 1), loc='upper left')

    fig, ax = dark_background(fig, ax)
    sc.plot_save(fig, fname='h_Ra_'+regime, fig_path=fig_path+'slides/', fig_fmt=fig_fmt, facecolor=fig.get_facecolor())

    """ model vs data """
    fig, ax = sc.plot_model_data(Ra_ls, eta_ls, regime_grid=regime_grid_td, t1_grid=t1_grid, load_grid=True,
                                      end_grid=end_grid, literature_file=None, legend=False, cmap='winter',
                                      postprocess_kwargs=postprocess_kwargs, c='xkcd:off white', averagescheme='timefirst',
                                      ylim=[6e-3, 2e-2], which_x='Ra_i_eff', which_h='rms', data_path=data_path,
                                      save=False, cbar='eta', include_regimes=['chaotic'], axissize=axissize,
                                      intercept=False, fig_fmt=fig_fmt)
    ax.tick_params(axis='x', labelsize=ticksize, pad=15)
    ax.tick_params(axis='y', labelsize=ticksize, pad=15)

    fig, ax = dark_background(fig, ax)
    sc.plot_save(fig, fname='model_data_'+regime, fig_path=fig_path+'slides/', fig_fmt=fig_fmt, facecolor=fig.get_facecolor())