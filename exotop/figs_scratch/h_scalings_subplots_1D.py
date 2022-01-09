import numpy as np
import matplotlib.pyplot as plt

from model_1D import the_results as plottop
from model_1D import inputs as ins
from model_1D import parameters as p
from model_1D.astroenvironment import radius_zeng, grav
from useful_and_bespoke import dark_background, imscatter
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
from datetime import date
import pickle as pkl

today = date.today().strftime("%b-%d-%Y")

# set paths
fig_path = '/home/claire/Works/exo-top/exotop/figs_ms/'
fig_format = '.pdf'  # '.png'
benchmark_path = '../benchmarks/'
planet_icon_path = '/home/claire/Pictures/science-graphics/planet_png/'
labelsize = 30
legsize = 22  # 30
xlabelpad = 50
ticksize = 25
linec = 'xkcd:ocean green'  # 'xkcd:british racing green'  # '#d88868'
alpha = 0.3  # 0.3

planet_kwargs = {
    # 'visc_type': 'KW'
}
run_kwargs = {
    # 'rms_type': 'Ra_i_eff'
}

names = ['Ea', 'eta_pre',
                   ]
mini_mc = [200e3, 2.63e10, #240e3, 1.6e10
           # 0.3
           ]
maxi_mc = [300e3, 5.32e13, #340e3, 2.6e12,
           # 3
           ]
# mini = [300e3, 1.6e11,
#                    ]
# maxi = [300e3, 1.6e11,
#                    ]

# how h varies across key input parameters
x_vars = ['t',
          'M_p', 'CMF', 'x_Eu']  # 'H_f'
units = [' Gyr',
         r' $M_{\oplus}$', ' CMF', r'$\%$ U, Th']  # pW kg$^{-1}$
log = [False,
       True, False, False]
x_range = [(2, 4.5),
           (0.1 * p.M_E, 5 * p.M_E), (0.1, 0.4), (0.3, 3)]  # (2e-12, 10e-12)
xticks = [[2, 3, 4],
          [0.1, 1, 5], [0.1, 0.2, 0.3, 0.4], [30, 100, 300]]  # [2, 5, 10]
xscales = [p.sec2Gyr,
           p.M_E ** -1, 1, 1e2]  # 1e12
xlabels = ['Age\n(Gyr)',
           'Planet mass\n' + r'($M_{\oplus}$)', 'Core Mass Fraction',
           'U and Th budget\n($\%$ relative to solar)']  # 'Radiogenic heating\n(pW kg$^{-1}$)'

dist_res = 5
x_res = 7
n_sigma = 1

fig, axes = plottop.plot_change_with_observeables_ensemble(dist_res=dist_res, x_res=x_res, n_sigma=n_sigma,
                                                           names=names, mini=mini_mc, maxi=maxi_mc,
                                                           defaults='baseline', age=4.5,
                                                           ticksize=ticksize, labelsize=labelsize, fig_height=6,
                                                           legend=True, lw=4, ylabel=r'RMS topography (m)',
                                                           labelpad=20, legendtop=True, tickwidth=1,
                                                           save=False,
                                                           update_kwargs=planet_kwargs, run_kwargs=run_kwargs,
                                                           model_param='dyn_top_rms', labels=[''],
                                                           x_vars=x_vars, units=units, log=log, x_range=x_range,
                                                           xscales=xscales, xlabels=xlabels,
                                                           linec=linec, leg_loc='upper right',
                                                           textc='k',  # 'xkcd:off white',
                                                           alpha=alpha, legsize=legsize, xlabelpad=xlabelpad,
                                                           # verbose=True
                                                           )

# # add second scaling relationship
# fig, axes = plottop.plot_change_with_observeables_ensemble(dist_res=10, x_res=5,
#                                                            defaults='baseline', age=4.5,
#                                                            ticksize=ticksize, labelsize=labelsize, fig_height=6,
#                                                            legend=False, lw=4, ylabel=r'RMS topography (m)',
#                                                            labelpad=20, legendtop=True, tickwidth=2,
#                                                            fig=fig, axes=axes, save=False,
#                                                            update_kwargs=planet_kwargs, run_kwargs=run_kwargs,
#                                                            model_param='dyn_top_rms_1param', labels=[''],
#                                                            x_vars=x_vars, units=units, log=log, x_range=x_range,
#                                                            xscales=xscales, xlabels=xlabels,
#                                                            linec='#749af3', leg_loc='upper right',
#                                                            textc='k',  # 'xkcd:off white',
#                                                            alpha=alpha, legsize=legsize, xlabelpad=xlabelpad,
#                                                            )

for i, ax in enumerate(axes):
    ax.set_xlim([x * xscales[i] for x in x_range[i]])
    ax.set_xticks(xticks[i])
    ax.set_yscale('log')
    # ax.set_ylim((10, 30000))
    # ax.set_yticks((200, 800, 1000))
    ax.set_ylim((10, 800))
    ax.set_yticks((10, 100, 800))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
axes[0].yaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
axes[0].yaxis.set_minor_formatter(ticker.NullFormatter())
# axes[0].set_yticks([600, 700, 800, 900, 1000, 1200,  1500])

axes[0].set_xlim(x_range[0])  # ensure time axis lims

# VENUS: 850 m
# h_Venus = 865.4906656355711
# M_Venus = 0.815
# h_Mars = 6688.627942023225
# M_Mars = 0.107

# ax = axes[0]
# imscatter(M_Venus, h_Venus, planet_icon_path + 'Venus.png', zoom=0.04, ax=ax)
# imscatter(M_Mars, h_Mars, planet_icon_path + 'Mars.png', zoom=0.08, ax=ax)
# ax.scatter(M_Venus, h_Venus, marker='$V$', c='xkcd:goldenrod', s=200, zorder=100)

# # Huang cases 1-13, 15
# h_Huang = np.array([200.15279436132423 , 688.2014927583677 , 673.7880493468331 , 402.07565967751117 , 695.2136989391211 , 672.4561163950626 , 214.12066607342535 , 488.4601789919337 , 878.5607285545191 , 292.43829959982384 , 311.3352436867767 , 339.3664129742059 , 640.1361418805931 , 430.1894190342128 ])
# for h in h_Huang:
#     ax.scatter(M_Venus, h, marker='^', s=70, alpha=0.5, c='xkcd:orchid', label=r'Huang+ (2013)', zorder=1)


handles = [
    mlines.Line2D([], [], color=linec, ls='-',
                  markersize=0, lw=4, label=r'$h^\prime_{\rm rms}$ = $f$(Ra$_i$, $b$)'),
    mlines.Line2D([], [], color='#749af3', ls='-',
                  markersize=0, lw=4, label=r'$h^\prime_{\rm rms}$ = $f$(Ra$_{i, {\rm eff}}$)'),
    # mlines.Line2D([], [], color='xkcd:goldenrod', marker='$V$',
    #               markersize=15, lw=0, label=r'Venus'),
    #                    mlines.Line2D([], [], color='xkcd:orchid', marker='^',
    #                                  markersize=15, lw=0, label=r'Huang+ (2013) 3D model'),
]

# legend?
# legend = axes[0].legend(handles=handles, frameon=False, fontsize=legsize,
#                         borderaxespad=0,
#                         loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=3, )


# rock strength scaling
def h_peak_rock(M=None, Y=100e6, rho_c=2700, C=1 / 2, **kwargs):
    # C is 1/3 to 1/2 - min stress difference supported elastically underneath load (Jeffreys' theorem)
    R = radius_zeng(M, CMF=0.3)  # in km
    g = grav(M, R)
    return (C ** -1 * Y) / (rho_c * g)

# m = np.linspace(0.1, 5, 10)
# h_max = h_peak_rock(m*p.M_E)
# axes[1].plot(m, h_max, c='k', alpha=0.6)
# print('rock strength', h_max, 'm')


# fig, *axes = dark_background(fig, axes)
plt.subplots_adjust(wspace=0.2)
# plt.tight_layout()

with open(fig_path + r"h_parameters_fig.pkl", "wb") as f:
    pkl.dump(fig, f)

# plt.suptitle(r'Ra$_{i, {\rm eff}}$ scaling')
fig.savefig(fig_path + 'h_parameters' + today + fig_format, bbox_inches='tight')
# plt.show()
