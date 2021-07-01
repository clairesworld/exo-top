""" ASPECT runs: plot scalings of h vs. Ra or T-heuristic, with subplots by eta """

import matplotlib.pyplot as plt
from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    c_regimes_td, fig_fmt, regime_grid_td, regime_names_td, load_grid, p_Earth, postprocess_kwargs
from postaspect import plt_aspect as plat
from useful_and_bespoke import colourised_legend
import matplotlib.ticker as ticker
from matplotlib import rc
from datetime import date

"""setup"""

rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})
today = date.today().strftime("%b-%d-%Y")
load = True  # load_grid
labelsize = 16
ms = 10
legsize = 12
ylabelpad = 7
xlabelpad = 10
mark = 'o'
vmin = None  # 1
vmax = None  # 3
cmap = None  # 'rainbow'  # 'art-nouveau-03'
c_rms = ['k', 'xkcd:seaweed', 'xkcd:lilac', 'xkcd:orange',
         'xkcd:yellow']  # update: must correspond to entire z range i.e. eta_ls
cleglabels = [r'$\Delta \eta = 10^{5}$', r'$\Delta \eta = 10^{6}$', r'$\Delta \eta = 10^{7}$',
              r'$\Delta \eta = 10^{8}$', r'$\Delta \eta = 10^{9}$']  # these correspond to entire z range as well

include_regimes = ['chaotic']  # , 'not ready'
averagescheme = 'timefirst'
which_x = 'Ra_i_eff'
hlim = (6e-3, 1.2e-2)

fig, axes = plt.subplots(1, 2, figsize=(10, 5))

"""all eta on one axis"""

fig, ax0 = plat.plot_h_vs(Ra=Ra_ls, eta=eta_ls, t1_grid=t1_grid, end_grid=end_grid, load_grid=load, data_path=data_path,
                          fig_path=fig_path, averagescheme=averagescheme, p_dimensionals=None, which_x=which_x,
                          beta0=[0.1, -0.1], sigma=1, fiterror=False, legend=False, fig=fig, ax=axes[0],
                          include_regimes=include_regimes, save=False, labelsize=labelsize,
                          xlabel=r'Ra$_{i,{\rm eff}}$', ylabel=r'Dynamic topography $h^\prime_{\rm rms}$',
                          legsize=legsize,
                          cleglabels=cleglabels,
                          title='', showpeak=False, vmin=vmin, vmax=vmax,
                          cmap=cmap, c_rms=c_rms,
                          fit=True, logx=True, logy=True, ms=ms, xlabelpad=xlabelpad, ylabelpad=ylabelpad,
                          show_isoviscous=False, ylim=hlim, xlim=None, postprocess_kwargs=postprocess_kwargs,
                          regime_grid=regime_grid_td, errortype='standard', cbar=False)

"""model versus data"""

fig, ax1 = plat.plot_model_data_errorbars(Ra_ls, eta_ls, regime_grid=regime_grid_td, t1_grid=t1_grid, load_grid=load,
                                          end_grid=end_grid, legend=False, ms=ms, fig=fig, ax=axes[1],
                                          postprocess_kwargs=postprocess_kwargs, averagescheme=averagescheme,
                                          ylim=hlim, which_x=which_x, which_h='rms', data_path=data_path,
                                          clist=c_rms, cmap=cmap, z_name='eta',
                                          save=False, include_regimes=include_regimes, errs=[0.2, 0.1, 0.05],
                                          fig_fmt=fig_fmt, vmin=vmin, vmax=vmax, xlabelpad=xlabelpad, ylabelpad=3,
                                          show_cbar=False, errortype='standard',
                                          title='',
                                          ylabel=r'Model',
                                          xlabel=r'Data')

for ax in axes:
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())
axes[1].xaxis.set_minor_formatter(ticker.ScalarFormatter())
axes[1].xaxis.set_major_formatter(ticker.ScalarFormatter())

# show colours outside
ax1 = colourised_legend(axes[1], clist=c_rms[1:], cleglabels=cleglabels[1:], lw=0, marker=mark, markersize=ms,
                        legsize=legsize, ncol=1)

fig.suptitle(r'Fit to $C$ Ra$_{i,eff}^p$', fontsize=labelsize)
plt.tight_layout()
plt.subplots_adjust(wspace=0.35)
fig.savefig(fig_path + 'h_Ra_scalings' + today + '.png', bbox_inches='tight')
