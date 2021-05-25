# GOOD COPY OF PLOT FOR SLIDES
import sh_things as sh
import matplotlib.pyplot as plt
from useful_and_bespoke import dark_background, cmap_from_ascii
# import postaspect.plt_aspect as plat
import model_1D.the_results as results
import model_1D.parameters as p
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
import numpy as np

cmap_path = '/home/claire/Works/exo-top/exotop/plot/cmaps/'
cmap_name = 'c3t3a'
cmap = cmap_from_ascii(cmap_name, path=cmap_path, end='.txt').reversed()
fig_path = '/home/claire/Works/exo-top/exotop/figs_scratch/'
data_path = '/home/claire/Works/aspect/runs/model-output/'
case = 'Ra1e8-eta1e7-wide'
# d, dT, alpha = 1, 1, 1
d, dT, alpha = 2890, 3000, 3e-5  # Hoggard AGU Monograph dim factors
labelsize = 20

# only do this once
# sh.make_model_spectrum(case, R=2, data_path=data_path, fig_path='', newfname='base_spectrum',
#                         bl_fudge=2*np.pi, max_dscale=1, plot=False, verbose=False)
# l, S = sh.load_model_spectrum_pkl(path='')
#
# # original
# h_ratio = 1
# clm = sh.random_harms_from_psd(S, l, R=2, h_ratio=h_ratio, plot=False, verbose=False)
#
# h_rms, h_peak = sh.coeffs_to_grid(clm, plot_grid=False, plot_spectrum=True, d=d, alpha_m=alpha, dT=dT, R=2*d,
#                                   verbose=False, cbar=False, labelsize=labelsize, cmap='nipy_spectral', cline='xkcd:off white')
#
# fig = plt.gcf()
# ax = plt.gca()
# ax.set_xlabel("Spherical harmonic degree", fontsize=labelsize, labelpad=16)
# ax.set_ylabel("Power (km$^2$ km$^2$)", fontsize=labelsize, labelpad=16)
# ax.set_xticks([2, 5, 10, 20, 50])
# # ax.xaxis.set_minor_formatter(ticker.ScalarFormatter())
# ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
# fig, ax = dark_background(fig, ax)
# fig.savefig('psd_example.png', bbox_inches='tight', transparent=True)
# # plt.show()

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
slides = False
nplanets = 2
n_stats = 1000

if slides:
    textc = 'xkcd:off white'
else:
    textc = 'k'

xlabel = 'Planet mass ' + r'($M_{\oplus}$)'
labelsize = 36 #33.5  # 30
legsize = 30 # 24  # 20
ticksize = 22  # 20
clabelpad = 35
print('\nfirst call')
fig, axes = results.plot_ocean_capacity_relative(n_stats=n_stats, relative=True, nplanets=nplanets, version=0,
                                                 legsize=legsize, ticksize=ticksize, labelsize=labelsize, wspace=0.15,
                                                 titlesize=32, fig_path=fig_path, save=False, log=True, alpha_w=0.3,
                                                 vol_0='Earth', simple_scaling=False,
                                                 defaults='Venusbaseline', textc=textc,
                                                 # title='Water volume to submerge land',
                                                 spectrum_fpath='/home/claire/Works/exo-top/exotop/figs_scratch/',
                                                 # benchmark_path+'wei_Venus/',
                                                 spectrum_fname='base_spectrum_l1.pkl',
                                                 #                                                  c='#81f79f',
                                                 c='xkcd:bordeaux', cmap=cmap,
                                                 alpha=1, lw=4, ymin=0.3, ymax=1.8, labelpad=10,
                                                 set_ylim=True, x_vars=['M_p'], units=['$M_E$'],
                                                 x_range=[(0.1 * p.M_E, 5 * p.M_E)], xscales=[p.M_E ** -1],
                                                 xlabels=[xlabel],
                                                 leg_bbox=(0, 1.01), clabelpad=clabelpad,
                                                 fname='ocean-vol', ytitle=1.05,  # vmax=3e-3,
                                                 mass_frac_sfcwater=np.logspace(-5, -3.4, num=30),
                                                 # [1e-5, 3e-5, 1e-4, 3e-4, 1e-3]
                                                 )
print('\nsecond call')
fig, axes = results.plot_ocean_capacity_relative(n_stats=n_stats, relative=True, nplanets=nplanets,
                                                 fig=fig, axes=axes, vol_0='Earth',  # 8.468613612559923e+17,
                                                 legsize=legsize, ticksize=ticksize, labelsize=labelsize, wspace=0.15,
                                                 titlesize=32, fig_path=fig_path, save=False,
                                                 simple_scaling=False, log=True,
                                                 defaults='Venusbaseline',
                                                 # title='Water volume to submerge land',
                                                 spectrum_fpath='/home/claire/Works/exo-top/exotop/figs_scratch/',
                                                 spectrum_fname='Venus_spectrum_l1.pkl',
                                                 #                                                  c='#81f79f',
                                                 c='xkcd:squash', ls='--',
                                                 alpha=1, lw=4, ymin=0.3, ymax=1.8, labelpad=10,
                                                 set_ylim=True, x_vars=['M_p'], units=['$M_E$'],
                                                 x_range=[(0.1 * p.M_E, 5 * p.M_E)], xscales=[p.M_E ** -1],
                                                 xlabels=[xlabel],
                                                 leg_bbox=(0, 1.01), clabelpad=clabelpad, ytitle=1.05, vmax=3e-3,
                                                 mass_frac_sfcwater=None)

print('\nthird call')
fig, axes = results.plot_ocean_capacity_relative(n_stats=n_stats, relative=True, nplanets=nplanets,
                                                 fig=fig, axes=axes, vol_0='Earth',  # 8.468613612559923e+17,
                                                 legsize=legsize, ticksize=ticksize, labelsize=labelsize, wspace=0.15,
                                                 titlesize=32, fig_path=fig_path, save=False,
                                                 simple_scaling=False, log=True,
                                                 defaults='Venusbaseline',
                                                 spectrum_fpath='/home/claire/Works/exo-top/exotop/figs_scratch/',
                                                 spectrum_fname='spectrum_-2.pkl',
                                                 c='xkcd:reddish orange', ls='-.',
                                                 alpha=1, lw=4, ymin=0.3, ymax=1.8, labelpad=10,
                                                 set_ylim=True, x_vars=['M_p'], units=['$M_E$'],
                                                 x_range=[(0.1 * p.M_E, 5 * p.M_E)], xscales=[p.M_E ** -1],
                                                 xlabels=[xlabel],
                                                 leg_bbox=(0, 1.01), clabelpad=clabelpad, ytitle=1.05, vmax=3e-3,
                                                 mass_frac_sfcwater=None)

ax = axes[0]
ax.axhline(y=1, c='xkcd:off white', alpha=0.5, zorder=0)
# for ax in axes:
#     ax.set_xscale('log')
#     ax.set_yscale('log')
ax.set_ylabel('Basin capacity (Earth oceans)', fontsize=labelsize, c=textc,
              labelpad=20)
ax.set_xlabel(xlabel, fontsize=labelsize, c=textc,
              labelpad=20)

ax.set_xlim((0.1, 5))
ax.set_ylim((1e-1, 1e0))
# ax.set_ylim((1e-2, 4e0))

ax.text(0.03, 0.97, '4.5 Ga\n300 kJ mol$^{-1}$\n0.3 CMF\n4.6 pW kg$^{-1}$', fontsize=legsize,
        horizontalalignment='left', c=textc,
        verticalalignment='top',
        transform=ax.transAxes)

ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
ax.xaxis.set_minor_formatter(ticker.NullFormatter())
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
ax.yaxis.set_minor_formatter(ticker.NullFormatter())
ax.set_xticks([0.1, 1, 2, 3, 4, 5])
ax.set_yticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1])
# ax.set_yticks([0.2, 0.3, 1, 2])

handles = [mlines.Line2D([], [], color='xkcd:bordeaux', ls='-', lw=3,
                         label='Pure dynamic topography'),
           mlines.Line2D([], [], color='xkcd:squash', ls='--', lw=3,
                         label='Venus-like topography'),
           mlines.Line2D([], [], color='xkcd:reddish orange', ls='-.', lw=3,
                         label='Red noise topography'),
           # mlines.Line2D([], [], color='g', ls='--', lw=3,
           #               label='Simple scaling')
           ]
ax.legend(handles=handles, bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left", frameon=False, fontsize=22, ncol=1)

if slides:
    fig, *axes = dark_background(fig, axes)
fig.savefig(fig_path + 'ocn_vol_v3.png', bbox_inches='tight')
plt.show()
