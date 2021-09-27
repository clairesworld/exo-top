# GOOD COPY OF PLOT FOR SLIDES
import sh_things as sh
import matplotlib.pyplot as plt
from useful_and_bespoke import dark_background, cmap_from_ascii, get_continuous_cmap, minmaxnorm, cornertext
# import postaspect.plt_aspect as plat
import model_1D.the_results as results
import model_1D.evolve as evol
import model_1D.parameters as p
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
import numpy as np
from postaspect.plt_aspect import plot_error_contours
from datetime import date

today = date.today().strftime("%b-%d-%Y")
fig_path = '/home/claire/Works/exo-top/exotop/figs_ms/'
data_path = '/home/claire/Works/aspect/runs/model-output/'
case = 'Ra1e8-eta1e7-wide'
# d, dT, alpha = 1, 1, 1
d, dT, alpha = 2890, 3000, 3e-5  # Hoggard AGU Monograph dim factors

""" test single planet """
# from model_1D import evolve as evol
# from model_1D import oceans
# from model_1D import parameters
# pl0 = \
# evol.bulk_planets(n=1, name='M_p', mini=parameters.M_E, maxi=parameters.M_E, like='Venusbaseline',  # verbose=True,
#                   t_eval=None, random=False, postprocessors=['topography'],)[0]
#
# spectrum_fname='base_spectrum.pkl'
# spectrum_fpath='/home/claire/Works/exo-top/exotop/figs_scratch/'
# degree, phi0 = sh.load_model_spectrum_pkl(fname=spectrum_fname, path=spectrum_fpath)
# print('rms of original phi0', sh.parseval_rms(phi0, sh.l_to_k(degree, 2)))
# print('rms from aspect model', pl0.dyn_top_aspect_prime[-1])
# print('dimensional rms', pl0.dyn_top_rms[-1])
#
# pl0 = oceans.max_ocean(pl0, at_age=4.5, name_rms='dyn_top_aspect_prime', phi0=phi0, n_stats=1)
# vol_0 = pl0.max_ocean
# print('basin capacity in M_E', vol_0*1000/parameters.M_E)


""" money plot """
spec_path = '/home/claire/Works/exo-top/exotop/top_spectra/'
cmap_path = '/plot/cmaps/'
# cmap_name = 'c3t3a'
# cmap = cmap_from_ascii(cmap_name, path=cmap_path, end='.txt').reversed()
hex_list = ['#ffd8e2', '#dce5e9', '#d6fcdd', '#bbfde1', '#a4fde1', '#49fffc', '#4cd5e8', '#56c3e0', '#72b2d3',
            '#3b5e71']
float_list = [1e-5, 5e-5, 6e-5, 8e-5, 1e-4, 1.1e-4, 1.3e-4, 2e-4, 2.3e-4, 3e-4]
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
    alpha_w = 0.3
    alpha_dist = 0.15
    cmap = get_continuous_cmap(hex_list, N=12, float_list=fln)
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

nplanets = 9
n_stats = 100
dist_res = 500
n_sigma = 1


fig, axes = plt.subplots(1, 3, figsize=(30, 10))
rad_vals = [0.3, 1, 3]
# c_spec = [c_dt, 'xkcd:squash', 'xkcd:reddish orange']
c_spec = [c_dt, 'xkcd:squash', 'xkcd:reddish orange']
ls_rad = ['--', '-', '-.']
labels_spec = ['Pure dynamic topography', 'Spectrally Venus-like topography', 'Red noise topography']
for ii, spec in enumerate(['base_spectrum_l1.pkl', 'Venus_spectrum_l1.pkl', 'spectrum_-2.pkl']):
    for jj, rad in enumerate(rad_vals):
        if jj == 2 and ii == 2:
            show_cbar = True
        else:
            show_cbar = False
        planet_kwargs.update({'x_Eu': rad})
        fig, ax, planets_l = results.plot_ocean_capacity(fig=fig, axes=axes[ii], M0=1,
                                              mass_frac_sfcwater=np.logspace(-6, -3, num=60), #vmin=1e-6, vmax=1e-2,
                                              textc=textc, titlesize=32,
                                              save=False, spectrum_fname=spec, c=c_spec[ii], ls=ls_rad[jj],
                                              spectrum_fpath=spec_path,
                                              ticksize=ticksize, labelsize=labelsize, clabelpad=clabelpad,
                                              relative=True,
                                              vol_0='Earth', simple_scaling=False, leg_bbox=(0, 1.01), log=True,
                                              x_range=[(0.1 * p.M_E, 5 * p.M_E)], nplanets=nplanets, dist_res=dist_res,
                                              cmap=cmap, version=0, alpha_w=alpha_w, alpha_dist=alpha_dist, alpha=1,
                                              names_mc=names_mc, mini_mc=mini_mc, maxi_mc=maxi_mc, defaults=baseline,
                                              show_contours=None, ensemble=True, n_stats=n_stats, verbose=True,
                                              update_kwargs=planet_kwargs, run_kwargs=run_kwargs, wspace=0.15,
                                              legend=False, lw=4, labelpad=10, show_cbar=show_cbar,
                                              x_vars=['M_p'], units=['$M_E$'], xscales=[p.M_E ** -1],
                                              xlabels=[xlabel], n_sigma=n_sigma, legsize=legsize,
                                              )

        # hacky but will run faster - need to rearrange order of rad loop and spec loop
        # planets = planets_l[0]  # because it's a list over axes
        # degree, phi0 = sh.load_model_spectrum_pkl(fname=spec, path=spec_path)
        # for pl in planets:
        #     # recalculate ocn
        #     pl = evol.postprocess_planet(pl, postprocessors=['ocean_capacity'], phi0=phi0, n_stats=n_stats)

        # show single test case
        # fig, axes = results.plot_ocean_capacity_relative(fig=fig, axes=axes, n_stats=n_stats, relative=True, nplanets=nplanets, version=0,
        #                                                  update_kwargs=planet_kwargs, run_kwargs=run_kwargs,
        #                                                  legsize=legsize, ticksize=ticksize, labelsize=labelsize, wspace=0.15,
        #                                                  titlesize=32, fig_path=fig_path, save=False, log=True, alpha_w=alpha_w,
        #                                                  vol_0='Earth', simple_scaling=False, M0=1, legend=False,
        #                                                  defaults=baseline, textc=textc,
        #                                                  # title='Water volume to submerge land',
        #                                                  spectrum_fpath='/home/claire/Works/exo-top/exotop/figs_scratch/',
        #                                                  # benchmark_path+'wei_Venus/',
        #                                                  spectrum_fname='base_spectrum_l1.pkl',
        #                                                  #                                                  c='#81f79f',
        #                                                  c=c_dt, cmap=cmap,
        #                                                  alpha=1, lw=4, ymin=0.3, ymax=1.8, labelpad=10,
        #                                                  set_ylim=True, x_vars=['M_p'], units=['$M_E$'],
        #                                                  x_range=[(0.1 * p.M_E, 5 * p.M_E)], xscales=[p.M_E ** -1],
        #                                                  xlabels=[xlabel],
        #                                                  leg_bbox=(0, 1.01), clabelpad=clabelpad,
        #                                                  fname='ocean-vol', ytitle=1.05,  # vmax=3e-3,
        #                                                  mass_frac_sfcwater=np.logspace(-5, -3, num=30),
        #                                                  ensemble=False,
        #                                                  # [1e-5, 3e-5, 1e-4, 3e-4, 1e-3]
        #                                                  show_contours=None  # [3e-5, 1e-4, 3e-4]
        #                                                  )

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
    ax.set_title(labels_spec[iax], fontsize=labelsize, c=textc)

axes[0] = cornertext(axes[0], 'WATER\nPLANETS', pos='top left', size=labelsize, pad=0.03)
axes[-1] = cornertext(axes[-1], 'LAND\nPLANETS', pos='bottom right', size=labelsize, pad=0.03)

# ax.text(0.03, 0.97, '4.5 Gyr\n300 kJ mol$^{-1}$\n0.3 CMF\n4.6 pW kg$^{-1}$', fontsize=legsize,
#         horizontalalignment='left', c=textc,
#         verticalalignment='top',
#         transform=ax.transAxes)

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
for jj, rad in enumerate(rad_vals):
    handles.append(mlines.Line2D([], [], color='k', ls=ls_rad[jj], lw=3, label=str(rad*100) + r'$\%$'))
axes[0].legend(handles=handles, bbox_to_anchor=(0, 1.05, 1, 0.2), loc="lower left", frameon=False, fontsize=legsize,
               title=r'\textbf{U, Th budget relative to solar}', title_fontsize=legsize, ncol=3)

# ax = axes[1]
# ax.set_xlabel('U and Th abundance\n($\%$ relative to solar)', fontsize=labelsize, c=textc,
#               labelpad=20)
# ax.set_xlim((30, 300))
# ax.set_xticks([30, 100, 200, 300])


if slides:
    fig, *axes = dark_background(fig, axes)

plt.subplots_adjust(wspace=0.1)
fig.savefig(fig_path + 'ocn_vol_ensemble_soft_' + today + '.png', bbox_inches='tight')
plt.show()
