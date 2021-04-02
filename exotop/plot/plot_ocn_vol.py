# GOOD COPY OF PLOT FOR SLIDES
import sh_things as sh
import matplotlib.pyplot as plt
from useful_and_bespoke import dark_background
import postaspect.plt_aspect as plat
import matplotlib.ticker as ticker
import numpy as np

fig_path = ''
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


""" money plot """

labelsize = 40
fig, axes = plat.plot_ocean_capacity_relative(n_stats=10, relative=True, nplanets=25,
                                                 legsize=20, ticksize=25, labelsize=labelsize, wspace=0.15,
                                                 titlesize=32, fig_path=fig_path, save=False,
                                                 showwaterscale=True, log=True,
                                                 defaults='Venusbaseline',  # textc='xkcd:off white',
                                                 title='Ocean volume to submerge land',
                                                 spectrum_fpath='/home/claire/Works/exo-top/exotop/figs_scratch/',
                                                 # benchmark_path+'wei_Venus/',
                                                 spectrum_fname='base_spectrum.pkl',  # 'model_power_m2_b.csv',
                                                 #                                                  c='#81f79f',
                                                 c='xkcd:light red',
                                                 alpha=1, lw=4, ymin=0.3, ymax=1.8, labelpad=10,
                                                 set_ylim=True, x_vars=['M_p'], units=['$M_E$'],
                                                 x_range=[(0.1 * p.M_E, 6 * p.M_E)], xscales=[p.M_E ** -1],
                                                 xlabels=['Planet mass\n($M_E$)'],
                                                 leg_bbox=(0, 1.01), clabelpad=70,
                                                 fname='ocean-vol', ytitle=1.05,
                                                 mass_frac_sfcwater=[1e-5, 3e-5, 1e-4, 3e-4, 1e-3])
ax = axes[0]
# for ax in axes:
#     ax.set_xscale('log')
#     ax.set_yscale('log')
ax.set_ylabel('$V_{\mathrm{max}}/V_{\mathrm{max},0}$', fontsize=labelsize,  # c='xkcd:off white',
              labelpad=20)
ax.set_xlabel('Planet mass ($M_E$)', fontsize=labelsize,  # c='xkcd:off white',
              labelpad=20)

ax.set_xlim((0.1, 6))
ax.set_ylim((3e-1, 1e1))

ax.text(0.05, 0.95, '4.5 Ga\n300 kJ mol$^{-1}$\n0.3 CMF\n4.6 pW kg$^{-1}$', fontsize=20,
        horizontalalignment='left',  # c='xkcd:off white',
        verticalalignment='top',
        transform=ax.transAxes)

ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
ax.xaxis.set_minor_formatter(ticker.NullFormatter())
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
ax.yaxis.set_minor_formatter(ticker.NullFormatter())
ax.set_xticks([0.1, 1, 2, 3, 4, 5, 6])

fig, *axes = dark_background(fig, axes)
fig.savefig(fig_path + 'ocn_vol-light.png', bbox_inches='tight')
