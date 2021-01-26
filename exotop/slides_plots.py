import sys
sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    fig_fmt, regime_grid_td, postprocess_kwargs, regime_names_td, load_grid, p_Earth    # noqa: E402
from exotop.postaspect import aspect_scalings as sc  # noqa: E402
from exotop.useful_and_bespoke import dark_background
import numpy as np
import matplotlib.lines as mlines

ticksize = 22
axissize = 40
c_fit = 'xkcd:off white'
lw = 5

regimes = 'all'
if regimes == 'all':
    include_regimes = ['steady', 'trans.', 'chaotic']
    ylim = [6e-3, 4e-2]
    xlim = [0.7e5, 3e7]
    yticks = [6e-3, 1e-2, 2e-2]
    xticks = [1e6, 1e7]
    c_rms = ['xkcd:lime green', 'xkcd:lilac', 'xkcd:orange', 'xkcd:yellow']
elif regimes == 'chaotic':
    include_regimes = ['steady', 'trans.', 'chaotic']
    ylim = ylim_ch
    xlim = xlim_ch
    yticks = yticks_ch
    xticks = xticks_ch
    c_rms = ['xkcd:lilac', 'xkcd:orange', 'xkcd:yellow']

fig, ax, ax2 = sc.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=True, data_path=data_path,
                 fig_path=fig_path, averagescheme='timefirst', which_x='Ra_i_eff',
                 sigma=1,
                 include_regimes=regimes, save=False, fname='h_Raieff_chaotic_timeavg', labelsize=axissize,
                       legend=False, figsize=(16,9), showpeak=False, lw=lw, ms=25, elw=2, ecapsize=8,
                 xlabel=r'Ra$_{i,eff}$', ylabel=r'dynamic topography $\Delta h^\prime$', ylabelpad=20,
                 title='', fiterror=False, c_fit=c_fit,
                 c_peak='k', c_rms=c_rms,
                       fit=True, logx=True, logy=True, hscale=1,
                 show_isoviscous=False, ylim=ylim, xlim=xlim, postprocess_kwargs=postprocess_kwargs,
                 regime_grid=regime_grid_td, p_dimensionals=p_Earth)

ax.tick_params(axis='x', labelsize=ticksize, pad=15)
ax2.tick_params(axis='y', labelsize=ticksize, pad=15)
ax.tick_params(axis='y', labelsize=ticksize, pad=15)
if yticks is not None:
    ax.set_yticks(yticks)
if xticks is not None:
    ax.set_xticks(xticks)

handles = []
handles.append(mlines.Line2D([], [], color=c_fit, marker='*', ls='--',
                          markersize=0, lw=lw, label=r'$\Delta h = 0.094$ Ra$^{-0.151}$'))
handles.append(mlines.Line2D([], [], color=c_rms[0], marker='o', ls='--',
                          markersize=20, lw=0, label=r'$\Delta \eta = 10^{5}$'))
handles.append(mlines.Line2D([], [], color=c_rms[1], marker='o', ls='--',
                          markersize=20, lw=0, label=r'$\Delta \eta = 10^{6}$'))
handles.append(mlines.Line2D([], [], color=c_rms[2], marker='o', ls='--',
                          markersize=20, lw=0, label=r'$\Delta \eta = 10^{7}$'))
handles.append(mlines.Line2D([], [], color=c_rms[3], marker='o', ls='--',
                          markersize=20, lw=0, label=r'$\Delta \eta = 10^{8}$'))

ax.legend(handles=handles, frameon=False, fontsize=25, ncol=1, bbox_to_anchor=(1.05, 1), loc='upper left')

fig, ax, ax2 = dark_background(fig, ax, ax2)
sc.plot_save(fig, fname='h_Ra_'+regimes, fig_path=fig_path+'slides/', fig_fmt=fig_fmt, facecolor=fig.get_facecolor())