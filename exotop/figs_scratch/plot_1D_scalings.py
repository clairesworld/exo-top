import numpy as np
import matplotlib.pyplot as plt

from model_1D import the_results as plottop
from model_1D import inputs as ins
from model_1D import parameters as p
from model_1D.astroenvironment import radius_zeng, grav
from useful_and_bespoke import dark_background, imscatter
import matplotlib.ticker as ticker
import matplotlib.lines as mlines

# set paths
fig_path = ''  # laptop
fig_format = '.png'
benchmark_path = '../benchmarks/'
planet_icon_path = '/home/claire/Pictures/science-graphics/planet_png/'
labelsize = 28
legsize = 20
xlabelpad = 20
ticksize = 20

# how h varies across key input parameters
x_vars = [#'t',
          'M_p', 'H_f', 'CMF']
units = [#'Gyr',
         '$M_E$', 'pW kg$^{-1}$', 'CMF']
log = [#False,
       True, False, False]
x_range = [#(2, 4.5),
           (0.1 * p.M_E, 6 * p.M_E), (2e-12, 10e-12), (0.1, 0.7)]
xticks = [#[2, 4.5],
          [0.1, 1, 6], [2, 5, 10], [0.1, 0.3, 0.7]]
xscales = [#p.sec2Gyr,
           p.M_E ** -1, 1e12, 1]
xlabels = [#'Age\n(Gyr)',
           'Planet mass\n($M_E$)', 'Rad. heating, $t_f$\n(pW kg$^{-1}$)',  'Core mass fraction']


fig, axes = plottop.plot_change_with_observeables_ensemble(age=4.5, dist_res=500, x_res=8,
                                                           defaults='baseline',
                                                           ticksize=ticksize, labelsize=labelsize, fig_height=6,
                                                           legend=True, lw=4, ylabel='$\Delta h_{rms}$ (m)',
                                                           labelpad=20, legendtop=True, tickwidth=2,
                                                           save=False, fname='relative_h_slides', fig_path=fig_path,
                                                           update_kwargs={'visc_type': 'KW'},
                                                           models=['dyn_top_rms'], labels=[''],
                                                           x_vars=x_vars, units=units, log=log, x_range=x_range,
                                                           xscales=xscales, xlabels=xlabels,
                                                           linec='#d88868', textc='xkcd:off white',
                                                           alpha=1, legsize=legsize)

for i, ax in enumerate(axes):
    ax.set_xlim([x*xscales[i] for x in x_range[i]])  # mass
    ax.set_xticks(xticks[i])
    ax.set_yscale('log')
    ax.set_ylim((400, 1200))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
axes[0].yaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
axes[0].yaxis.set_minor_formatter(ticker.NullFormatter())
axes[0].set_yticks([500, 600, 700, 800, 900, 1000])

# VENUS: 850 m
h_Venus = 865.4906656355711
M_Venus = 0.815
h_Mars = 6688.627942023225
M_Mars = 0.107

ax = axes[0]
# imscatter(M_Venus, h_Venus, planet_icon_path + 'Venus.png', zoom=0.04, ax=ax)
# imscatter(M_Mars, h_Mars, planet_icon_path + 'Mars.png', zoom=0.08, ax=ax)
ax.scatter(M_Venus, h_Venus, marker='$V$', c='xkcd:goldenrod', s=200, zorder=100)

# # Huang cases 1-13, 15
# h_Huang = np.array([200.15279436132423 , 688.2014927583677 , 673.7880493468331 , 402.07565967751117 , 695.2136989391211 , 672.4561163950626 , 214.12066607342535 , 488.4601789919337 , 878.5607285545191 , 292.43829959982384 , 311.3352436867767 , 339.3664129742059 , 640.1361418805931 , 430.1894190342128 ])
# for h in h_Huang:
#     ax.scatter(M_Venus, h, marker='^', s=70, alpha=0.5, c='xkcd:orchid', label=r'Huang+ (2013)', zorder=1)


handles = [mlines.Line2D([], [], color='#d88868', ls='-',
                         markersize=0, lw=4, label='$\Delta h^\prime = 0.12$ Ra$_{i, eff}^{-0.17}$'),
           #          mlines.Line2D([], [], color='#749af3', ls='-',
           #                                  markersize=0, lw=4, label='h = f($\delta_{rh}, \Delta T_{rh}$)'),
           mlines.Line2D([], [], color='xkcd:goldenrod', marker='$V$',
                         markersize=15, lw=0, label=r'Venus'),
           #                    mlines.Line2D([], [], color='xkcd:orchid', marker='^',
           #                                  markersize=15, lw=0, label=r'Huang+ (2013) 3D model'),
           ]

legend = axes[0].legend(handles=handles, frameon=False, fontsize=legsize,
                        borderaxespad=0,
                        loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=3, )

fig, *axes = dark_background(fig, axes)
plt.subplots_adjust(wspace=0.2)
# plt.tight_layout()
# plt.show()
fig.savefig(fig_path+'h_parameters_slides'+fig_format, bbox_inches='tight')
