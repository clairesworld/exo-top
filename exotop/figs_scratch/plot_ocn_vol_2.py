# GOOD COPY OF PLOT FOR SLIDES
import sh_things as sh
import matplotlib.pyplot as plt
from useful_and_bespoke import dark_background, cmap_from_ascii, get_continuous_cmap, minmaxnorm, cornertext
import model_1D.the_results as results
import model_1D.evolve as evol
import model_1D.parameters as p
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
import numpy as np
from datetime import date
from matplotlib import rc

rc('text', usetex=True)  # turn off for running over ssh
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'CMU Serif'
today = date.today().strftime("%b-%d-%Y")
fig_path = '/home/claire/Works/exo-top/exotop/figs_ms/'
# fig_path = '/home/cmg76/Works/exo-top/exotop/figs_ms/'
data_path = '/home/cmg76/Works/aspect/runs/model-output/'
case = 'Ra1e8-eta1e7-wide'
# d, dT, alpha = 1, 1, 1
d, dT, alpha = 2890, 3000, 3e-5  # Hoggard AGU Monograph dim factors


""" money plot """
spec_path = '/home/claire/Works/exo-top/exotop/top_spectra/'
# spec_path = '/home/cmg76/Works/exo-top/exotop/top_spectra/'
cmap_path = '/plot/cmaps/'
# cmap_name = 'c3t3a'
# cmap = cmap_from_ascii(cmap_name, path=cmap_path, end='.txt').reversed()
# hex_list = ['#ffd8e2', '#dce5e9', '#d6fcdd', '#bbfde1', '#a4fde1', '#49fffc', '#4cd5e8', '#56c3e0', '#72b2d3',
#             '#3b5e71']
hex_list = ['#ffd8e2', '#dce5e9', '#d6fcdd', '#bbfde1', '#a4fde1', '#49fffc', '#4cd5e8', '#56c3e0',
            '#3b5e71']
float_list = [1e-5, 5e-5, 6e-5, 8e-5, 1e-4, 2e-4, 2.3e-4, 3e-4]
fln = None  # minmaxnorm(float_list, a=0, b=1)

slides = False
xlabel = 'Planet mass ' + r'($M_{\oplus}$)'
labelsize = 30  # 33.5  # 30
legsize = 22
ticksize = 26  # 20
clabelpad = 50

if slides:
    textc = 'xkcd:off white'
    c_dt = 'xkcd:orange red'
    alpha_w = 0.8  # 0.6
    alpha_dist = 0.25
    cmap = 'gist_earth_r'
else:
    textc = 'k'
    c_dt = 'xkcd:bordeaux'
    alpha_w = 0.5 #0.3
    alpha_dist = 0.15
    cmap = get_continuous_cmap(hex_list, N=14, float_list=fln)
    cmap.set_under(hex_list[0])
    cmap.set_over(hex_list[-1])

mini_mc = [200e3, 2.63e10, #240e3, 1.6e10
           # 0.3
           ]
maxi_mc = [300e3, 5.32e13, #340e3, 2.6e12,
           # 3
           ]
names_mc = ['Ea', 'eta_pre',
            # 'x_Eu'
            ]
baseline = 'baseline'
planet_kwargs = {
    # 'T_s':730
}
run_kwargs = {
    # 'rms_type': 'Ra_i_eff'
}

nplanets = 7  # 32
n_stats = 2  # 500
dist_res = 1000  #500 # 1000
n_sigma = 1

pickledir = '/home/claire/Works/exo-top/exotop/figs_scratch/pickle/'

fig, axes = plt.subplots(1, 3, figsize=(30, 10))
rad_vals = [0.3, 1, 3]
peak_ratios = [3.5,  3.9]
# c_spec = [c_dt, 'xkcd:squash', 'xkcd:reddish orange']
c_spec = [c_dt, 'xkcd:reddish orange']
ls_rad = ['--', '-', '-.']
ls_spec = ['--', ':']
labels_spec = ['Pure dynamic topography',  'Red noise topography']
labels_cols = ['Low U and Th, cold', 'Solar U and Th', 'High U and Th, hot']
for ii, spec in enumerate(['base_spectrum_l1.pkl',  'spectrum_-2.pkl']):
    print('spectrum', ii+1, '/ 3')
    for jj, rad in enumerate(rad_vals):
        print('   rad budget', jj + 1, '/ 3')
        if jj == len(rad_vals) - 1 and (ii == len(c_spec) - 1):
            show_cbar = True
        else:
            show_cbar = False
        planet_kwargs.update({'x_Eu': rad})

        picklefile = pickledir + 'ocnplot-' + str(ii) + '-' + str(jj)
        fig, ax, planets_axes = results.plot_ocean_capacity(fig=fig, axes=axes[jj], M0=1,
                                              pickleto=picklefile,
                                              peak_ratio=peak_ratios[ii],
                                              mass_frac_sfcwater=np.logspace(-6, np.log10(3e-3), num=120),
                                              # mass_frac_sfcwater=np.logspace(-6, np.log10(2e-3), num=60), #vmin=1e-6, vmax=1e-2,
                                              textc=textc, titlesize=32,
                                              save=False, spectrum_fname=spec, c=c_spec[ii], ls=ls_spec[ii],
                                              spectrum_fpath=spec_path, extra_def=True,
                                              ticksize=ticksize, labelsize=labelsize, clabelpad=clabelpad,
                                              relative=True,
                                              vol_0='Earth', simple_scaling=False, leg_bbox=(0, 1.01), log=True,
                                              x_range=[(0.1 * p.M_E, 5 * p.M_E)], nplanets=nplanets, dist_res=dist_res,
                                              cmap=cmap, version=0, alpha_w=alpha_w, alpha_dist=alpha_dist, alpha=1,
                                              names_mc=names_mc, mini_mc=mini_mc, maxi_mc=maxi_mc, defaults=baseline,
                                              show_contours=None, ensemble=True, n_stats=n_stats, verbose=False,
                                              update_kwargs=planet_kwargs, run_kwargs=run_kwargs, wspace=0.15,
                                              legend=False, lw=4, labelpad=10, show_cbar=show_cbar,
                                              x_vars=['M_p'], units=['$M_E$'], xscales=[p.M_E ** -1],
                                              xlabels=[xlabel], n_sigma=n_sigma, legsize=legsize,
                                              )

# format
for iax, ax in enumerate(axes):
    # if slides:
    #     ax.axhline(y=1, c='xkcd:off white', alpha=0.5, zorder=0)
    # else:
    #     ax.axhline(y=1, c='k', alpha=0.5, zorder=0, lw=0.5)
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    # ax.set_ylim((2e-1, 1.8))  # old Ra i eff scaling
    ax.set_ylim((0.02, 1))  # low b scaling
    # ax.set_yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1, 1.4])  # old Ra i eff scaling
    # ax.set_yticks([0.01, 0.1, 0.2])  # low b scaling
    ax.set_yscale('log')  # ensure
    if iax > 0:  # remove tick labels but keep ticks
        tk = [item.get_text() for item in ax.get_yticklabels()]
        empty_string_labels = [''] * len(tk)
        ax.set_yticklabels(empty_string_labels)
    ax.set_xlabel(xlabel, fontsize=labelsize, c=textc, labelpad=20)
    ax.set_xlim((0.1, 5))
    ax.set_xticks([0.1, 1, 2, 3, 4, 5])
    if iax == 0:
        ax.set_ylabel('Basin capacity (Earth oceans)', fontsize=labelsize, c=textc,
                      labelpad=20)
    else:
        ax.set_ylabel('')
    ax.set_title(labels_cols[iax], fontsize=labelsize, c=textc)

axes[0] = cornertext(axes[0], 'WATER\nPLANETS', pos='top left', size=labelsize, pad=0.03)
axes[-1] = cornertext(axes[-1], 'LAND\nPLANETS', pos='bottom right', size=labelsize, pad=0.03)


# handles = [mlines.Line2D([], [], color=c_dt, ls='-', lw=3,
#                          label='Pure dynamic topography'),
#            mlines.Line2D([], [], color='xkcd:squash', ls='--', lw=3,
#                          label='Venus-like topography'),
#            mlines.Line2D([], [], color='xkcd:reddish orange', ls='-.', lw=3,
#                          label='Red noise topography'),
#            # mlines.Line2D([], [], color='g', ls='--', lw=3,
#            #               label='Simple scaling')
#            ]
handles = []
names = ['Pure dynamic topography', 'Red noise']  # str(rad*100) + r'$\%$')
for jj, n in enumerate(names):
    handles.append(mlines.Line2D([], [], color=c_spec[jj], ls=ls_spec[jj], lw=3, label=n))
axes[0].legend(handles=handles, bbox_to_anchor=(0, 1.05, 1, 0.2), loc="lower left", frameon=False, fontsize=legsize,
               title=r'\textbf{Topography spectral model}', title_fontsize=legsize, ncol=3)

# ax = axes[1]
# ax.set_xlabel('U and Th abundance\n($\%$ relative to solar)', fontsize=labelsize, c=textc,
#               labelpad=20)
# ax.set_xlim((30, 300))
# ax.set_xticks([30, 100, 200, 300])


if slides:
    fig, *axes = dark_background(fig, axes)

plt.subplots_adjust(wspace=0.1)
fig.savefig(fig_path + 'ocn_vol2_fast_1000_' + today + '.pdf', bbox_inches='tight')
# plt.show()



