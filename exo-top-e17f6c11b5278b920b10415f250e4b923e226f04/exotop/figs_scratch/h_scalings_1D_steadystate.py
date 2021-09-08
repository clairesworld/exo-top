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
linec = 'xkcd:ocean green'  # 'xkcd:british racing green'  # '#d88868'

# how h varies across key input parameters
x_vars = [#'t',
          'M_p', 'CMF', 'H_f']
units = [#'Gyr',
         '$M_E$', 'CMF', 'pW kg$^{-1}$']
log = [False,
       True, False, False]
x_range = [#(2, 5),
           (0.1 * p.M_E, 5 * p.M_E),  (0.1, 0.4), (2e-12, 10e-12)]
xticks = [#[2, 3, 4, 5],
          [0.1, 1, 5], [0.1, 0.2, 0.3, 0.4], [2, 5, 10]]
xscales = [#p.sec2Gyr,
           p.M_E ** -1, 1, 1e12]
xlabels = [#'Age\n(Gyr)',
           'Planet mass\n($M_E$)', 'Core mass fraction', 'Present rad. heating\n(pW kg$^{-1}$)']

fig, axes = plottop.plot_change_with_observeables(age=4.5, dist_res=1000, x_res=7, relative=False,
                                                           defaults='baseline',
                                                           ticksize=ticksize, labelsize=labelsize, fig_height=6,
                                                           legend=True, lw=4, ylabel='$h_{rms}$ (m)',
                                                           labelpad=20, legendtop=True, tickwidth=2,
                                                           save=False, fname='relative_h_slides', fig_path=fig_path,
                                                           # update_kwargs={'visc_type': 'KW'},
                                                           models=['dyn_top_rms'], labels=[''],
                                                           x_vars=x_vars, units=units, log=log, x_range=x_range,
                                                           xscales=xscales, xlabels=xlabels,
                                                           linec=linec,
                                                           textc='k', #'xkcd:off white',
                                                           alpha=0.3, legsize=legsize,
                                                           solve_ODE=False, steady=True,
                                                           )

for i, ax in enumerate(axes):
    ax.set_xlim([x*xscales[i] for x in x_range[i]])
    ax.set_xticks(xticks[i])
    ax.set_yscale('log')
    ax.set_ylim((350, 1100))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
axes[0].yaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
axes[0].yaxis.set_minor_formatter(ticker.NullFormatter())
axes[0].set_yticks([400, 500, 600, 700, 800, 900, 1000])



handles = [mlines.Line2D([], [], color=linec, ls='-',
                         markersize=0, lw=4, label='$\Delta h^\prime = 0.14$ Ra$_{i, eff}^{-0.18}$'),
           #          mlines.Line2D([], [], color='#749af3', ls='-',
           #                                  markersize=0, lw=4, label='h = f($\delta_{rh}, \Delta T_{rh}$)'),
           # mlines.Line2D([], [], color='xkcd:goldenrod', marker='$V$',
           #               markersize=15, lw=0, label=r'Venus'),
           #                    mlines.Line2D([], [], color='xkcd:orchid', marker='^',
           #                                  markersize=15, lw=0, label=r'Huang+ (2013) 3D model'),
           ]

# legend?
# legend = axes[0].legend(handles=handles, frameon=False, fontsize=legsize,
#                         borderaxespad=0,
#                         loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=3, )

# fig, *axes = dark_background(fig, axes)
plt.subplots_adjust(wspace=0.2)
# plt.tight_layout()
# plt.show()
fig.savefig(fig_path+'h_parameters_ss'+fig_format, bbox_inches='tight')
