# GOOD COPY OF PLOT FOR SLIDES
import sh_things as sh
import matplotlib.pyplot as plt
from useful_and_bespoke import dark_background
# import postaspect.plt_aspect as plat
import model_1D.the_results as results
import model_1D.parameters as p
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
import numpy as np

fig_path = '/home/claire/Works/exo-top/exotop/figs_scratch/'
data_path = '/home/claire/Works/aspect/runs/model-output/'
case = 'Ra1e8-eta1e7-wide'
# d, dT, alpha = 1, 1, 1
d, dT, alpha = 2890, 3000, 3e-5  # Hoggard AGU Monograph
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
print('first call')
labelsize = 40
fig, axes = results.plot_ocean_capacity_relative(n_stats=10, relative=True, nplanets=8,
                                                 legsize=20, ticksize=25, labelsize=labelsize, wspace=0.15,
                                                 titlesize=32, fig_path=fig_path, save=False,
                                                 showwaterscale=True, log=True,
                                                 vol_0='Earth', simple_scaling=False,
                                                 defaults='Venusbaseline', textc='xkcd:off white',
                                                 # title='Water volume to submerge land',
                                                 spectrum_fpath='/home/claire/Works/exo-top/exotop/figs_scratch/',
                                                 # benchmark_path+'wei_Venus/',
                                                 spectrum_fname='base_spectrum_l1.pkl',
                                                 #                                                  c='#81f79f',
                                                 c='xkcd:light red',
                                                 alpha=1, lw=4, ymin=0.3, ymax=1.8, labelpad=10,
                                                 set_ylim=True, x_vars=['M_p'], units=['$M_E$'],
                                                 x_range=[(0.1 * p.M_E, 6 * p.M_E)], xscales=[p.M_E ** -1],
                                                 xlabels=['Planet mass\n($M_E$)'],
                                                 leg_bbox=(0, 1.01), clabelpad=70,
                                                 fname='ocean-vol', ytitle=1.05, vmax=3e-3,
                                                 mass_frac_sfcwater=[1e-5, 3e-5, 1e-4, 3e-4, 1e-3])
print('second call')
fig, axes = results.plot_ocean_capacity_relative(n_stats=10, relative=True, nplanets=8,
                                                 fig=fig, axes=axes, vol_0='Earth',  # 8.468613612559923e+17,
                                                 legsize=20, ticksize=25, labelsize=labelsize, wspace=0.15,
                                                 titlesize=32, fig_path=fig_path, save=False,
                                                 simple_scaling=True,
                                                 showwaterscale=True, log=True,
                                                 defaults='Venusbaseline', textc='xkcd:off white',
                                                 # title='Water volume to submerge land',
                                                 spectrum_fpath='/home/claire/Works/exo-top/exotop/figs_scratch/',
                                                 spectrum_fname='Venus_spectrum_l1.pkl',
                                                 #                                                  c='#81f79f',
                                                 c='xkcd:hot pink', ls='--',
                                                 alpha=1, lw=4, ymin=0.3, ymax=1.8, labelpad=10,
                                                 set_ylim=True, x_vars=['M_p'], units=['$M_E$'],
                                                 x_range=[(0.1 * p.M_E, 6 * p.M_E)], xscales=[p.M_E ** -1],
                                                 xlabels=['Planet mass\n($M_E$)'],
                                                 leg_bbox=(0, 1.01), clabelpad=70, ytitle=1.05, vmax=3e-3,
                                                 mass_frac_sfcwater=None)

ax = axes[0]
ax.axhline(y=1, c='xkcd:off white', alpha=0.5, zorder=0)
# for ax in axes:
#     ax.set_xscale('log')
#     ax.set_yscale('log')
ax.set_ylabel('Relative basin capacity', fontsize=labelsize,  # c='xkcd:off white',
              labelpad=20)
ax.set_xlabel('Planet mass ($M_E$)', fontsize=labelsize,  # c='xkcd:off white',
              labelpad=20)

ax.set_xlim((0.1, 6))
ax.set_ylim((2e-1, 4e0))
# ax.set_ylim((1e-2, 4e0))

ax.text(0.05, 0.95, '4.5 Ga\n300 kJ mol$^{-1}$\n0.3 CMF\n4.6 pW kg$^{-1}$', fontsize=20,
        horizontalalignment='left', c='xkcd:off white',
        verticalalignment='top',
        transform=ax.transAxes)

ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
ax.xaxis.set_minor_formatter(ticker.NullFormatter())
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
ax.yaxis.set_minor_formatter(ticker.NullFormatter())
ax.set_xticks([0.1, 1, 2, 3, 4, 5, 6])
ax.set_yticks([0.3, 1, 3])

handles = [mlines.Line2D([], [], color='xkcd:light red', ls='-', lw=3,
                         label='Pure dynamic topography'),
           mlines.Line2D([], [], color='xkcd:hot pink', ls='--', lw=3,
                         label='Venus-like topography'),
           # mlines.Line2D([], [], color='g', ls='--', lw=3,
           #               label='Simple scaling')
           ]
ax.legend(handles=handles, bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left", frameon=False, fontsize=20, ncol=2)

fig, *axes = dark_background(fig, axes)
fig.savefig(fig_path + 'ocn_vol_test.png', bbox_inches='tight')
