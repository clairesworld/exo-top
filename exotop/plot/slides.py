""" ASPECT runs: good copy plots for some slides  """
from matplotlib import ticker

from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, \
     fig_fmt, regime_grid_td, postprocess_kwargs, regime_names_td, load_grid, p_Venus  # noqa: E402
from postaspect import plt_aspect as plat  # noqa: E402
from useful_and_bespoke import dark_background, cmap_from_list
import matplotlib.lines as mlines

ticksize = 22
axissize = 40
c_fit = 'xkcd:off white'
c_rms = ['r', 'xkcd:lime green', 'xkcd:lilac', 'xkcd:orange', 'xkcd:yellow']
lw = 5
ms = 25
elw = 2
ecapsize = 8

regimes = ['chaotic']
for regime in regimes:
    if regime == 'all':
        include_regimes = ['steady', 'trans.', 'chaotic']
        c_rms = c_rms[:]
        ylim = [6e-3, 4e-2]
        xlim = [0.7e5, 3e7]
        yticks = [6e-3, 1e-2, 2e-2, 4e-2]
        xticks = [1e5, 1e6, 1e7]
        fitlabel = r'$\Delta h = 0.345$ Ra$_{i, eff}^{-0.212}$'
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
        # c_rms = c_rms[1:]
        ylim = [6e-3, 1.1e-2]
        xlim = [2e6, 8e7]
        yticks = [6e-3, 8e-3, 1e-2]
        xticks = [3e6, 1e7]
        fitlabel = r'$\Delta h = 0.11$ Ra$_{i, eff}^{-0.16}$'
        handles = [mlines.Line2D([], [], color=c_fit, marker='*', ls='--',
                                 markersize=0, lw=lw, label=fitlabel),
                   mlines.Line2D([], [], color=c_rms[1], marker='o', ls='--',
                                 markersize=20, lw=0, label=r'$\Delta \eta = 10^{6}$'),
                   mlines.Line2D([], [], color=c_rms[2], marker='o', ls='--',
                                 markersize=20, lw=0, label=r'$\Delta \eta = 10^{7}$'),
                   mlines.Line2D([], [], color=c_rms[3], marker='o', ls='--',
                                 markersize=20, lw=0, label=r'$\Delta \eta = 10^{8}$'),
                   mlines.Line2D([], [], color=c_rms[4], marker='o', ls='--',
                                 markersize=20, lw=0, label=r'$\Delta \eta = 10^{9}$')
                   ]

    # fig, ax, ax2 = plat.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=True, data_path=data_path,
    #                  fig_path=fig_path, averagescheme='timefirst', which_x='Ra_i_eff', ticksize=ticksize,
    #                  sigma=1,
    #                  include_regimes=include_regimes, save=False,  labelsize=axissize,
    #                        legend=False, figsize=(16,9), showpeak=False, lw=lw, ms=ms, elw=elw, ecapsize=ecapsize,
    #                  xlabel=r'Ra$_{i,eff}$', ylabel=r'dynamic topography $\Delta h_{rms}^\prime$', ylabelpad=20, xlabelpad=13,
    #                  title='', fiterror=False, c_fit=c_fit,
    #                  c_peak='k', c_rms=c_rms,
    #                        fit=True, logx=True, logy=True, hscale=1,
    #                  show_isoviscous=False, ylim=ylim, xlim=xlim, postprocess_kwargs=postprocess_kwargs,
    #                  regime_grid=regime_grid_td,
    #                               # p_dimensionals=p_Venus, y2label='             Venus\n'+r'             $\Delta h_{rms}$ (m)',
    #                               )

    fig, ax = plat.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=True, data_path=data_path,
                   fig_path=fig_path, averagescheme='timefirst', p_dimensionals=None, which_x='Ra_i_eff', ticksize=ticksize,
                   beta0=[0.1, -0.15], sigma=1, fiterror=False, legend=False, ylabelpad=20,
                   include_regimes=['chaotic'], save=True, fname='h_Raieff_chaotic_timeavg', labelsize=axissize,
                   xlabel=r'Ra$_{i,eff}$', ylabel=r'topography, $\Delta h^\prime_{rms}$', #legsize=legsize, cleglabels=cleglabels,
                   title=r'Fit to $C$ Ra$_{i,eff}^p$', showpeak=False,
                   cmap=None, c_rms=c_rms, c_fit=c_fit,
                   fit=True, logx=True, logy=True, ms=ms, elw=elw,
                   show_isoviscous=False, ylim=None, xlim=None, postprocess_kwargs=postprocess_kwargs,
                   regime_grid=regime_grid_td, figsize=(13, 9), errortype='standard', cbar=False)

    ax.tick_params(axis='x', labelsize=ticksize, pad=15)
    ax.tick_params(axis='y', labelsize=ticksize, pad=15)
    if yticks is not None:
        ax.set_yticks(yticks)
    if xticks is not None:
        ax.set_xticks(xticks)
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())
    # ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    # ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
    # ax.ticklabel_format(style='plain', axis='y', useOffset=False)
    ax.legend(handles=handles, frameon=False, fontsize=25, ncol=1, bbox_to_anchor=(1.01, 1), loc='upper left')

    # ax2.tick_params(axis='y', labelsize=ticksize, pad=15)
    # ax2.set_yticks([1000, 1600, 2280])
    # ax2.yaxis.set_major_formatter(ticker.ScalarFormatter())
    # ax2.yaxis.set_minor_formatter(ticker.ScalarFormatter())

    # fig, _, _ = dark_background(fig, (ax, ax2))
    fig, _ = dark_background(fig, ax)
    plat.plot_save(fig, fname='h_Ra_'+regime, fig_path=fig_path+'slides/', fig_fmt=fig_fmt, facecolor=fig.get_facecolor())




    """ model vs data, rms """
    # c_contours = 'xkcd:off white'
    # fc = 'k'
    # fig, ax = plat.plot_model_data_errorbars(Ra_ls, eta_ls, errortype='standard',
    #                                          regime_grid=regime_grid_td, t1_grid=t1_grid, load_grid=True,
    #                                   end_grid=end_grid, literature_file=None, legend=False, cmap=None, ms=ms,
    #                                   postprocess_kwargs=postprocess_kwargs, c_contours=c_contours, fc=fc, averagescheme='timefirst',
    #                                   ylim=ylim, which_x='Ra_i_eff', which_h='rms', data_path=data_path,
    #                                   clist=c_rms, figsize=(10, 10), labelsize=axissize, ylabelpad=20, xlabelpad=13,
    #                                   z_name='eta', elw=elw, ecapsize=ecapsize,
    #                                   save=False, include_regimes=include_regimes, errorsize=20, errs=[0.5, 0.2, 0.1, 0.05],
    #                                   intercept=False, fig_fmt=fig_fmt, vmin=1, vmax=3, show_cbar=False,
    #                                        ylabel=r'Model $\Delta h_{rms}^\prime$', xlabel=r'Data $\Delta h_{rms}^\prime$')
    # ax.tick_params(axis='x', labelsize=ticksize, pad=15)
    # ax.tick_params(axis='y', labelsize=ticksize, pad=15)
    # if yticks is not None:
    #     ax.set_yticks(yticks)
    #     ax.set_xticks(yticks)
    # ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    # ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())
    # ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    # ax.xaxis.set_minor_formatter(ticker.ScalarFormatter())
    #
    # fig, ax = dark_background(fig, ax)
    # plat.plot_save(fig, fname='model_data_'+regime, fig_path=fig_path+'slides/', fig_fmt=fig_fmt, facecolor=fig.get_facecolor())





    """ model vs data, peak """
    # # yticks = [1e-2, 2e-2, 3e-2]
    # # fig, ax = plat.plot_model_data_errorbars(Ra_ls, eta_ls, regime_grid=regime_grid_td, t1_grid=t1_grid, load_grid=True,
    # #                                          end_grid=end_grid, literature_file=None, legend=False, cmap=None, ms=ms,
    # #                                          postprocess_kwargs=postprocess_kwargs, c_contours=c_contours, fc=fc,
    # #                                          averagescheme='timefirst',
    # #                                          ylim=[1e-2, 3e-2], which_x='Ra_i_eff', which_h='peak', data_path=data_path,
    # #                                          clist=c_rms, figsize=(10, 10), labelsize=axissize, ylabelpad=20,
    # #                                          xlabelpad=13,
    # #                                          z_name='eta', elw=elw, ecapsize=ecapsize,
    # #                                          save=False, include_regimes=include_regimes, errorsize=20,
    # #                                          errs=[0.5, 0.2, 0.1, 0.05],
    # #                                          intercept=False, fig_fmt=fig_fmt, vmin=1, vmax=3, show_cbar=False,
    # #                                          ylabel=r'Model $\Delta h_{peak}^\prime$',
    # #                                          xlabel=r'Data $\Delta h_{peak}^\prime$')
    # # ax.tick_params(axis='x', labelsize=ticksize, pad=15)
    # # ax.tick_params(axis='y', labelsize=ticksize, pad=15)
    # # if yticks is not None:
    # #     ax.set_yticks(yticks)
    # #     ax.set_xticks(yticks)
    # # ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    # # ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())
    # # ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    # # ax.xaxis.set_minor_formatter(ticker.ScalarFormatter())
    # #
    # # fig, ax = dark_background(fig, ax)
    # # plat.plot_save(fig, fname='model_data_peak_' + regime, fig_path=fig_path + 'slides/', fig_fmt=fig_fmt,
    # #                facecolor=fig.get_facecolor())
