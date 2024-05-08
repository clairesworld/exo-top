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
from matplotlib import rc


rc('text', usetex=True)  # turn off for running over ssh
today = date.today().strftime("%b-%d-%Y")
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'CMU Serif'

# set paths
fig_path = '/home/claire/Works/exo-top/exotop/figs_ms/'
fig_format = '.png'
benchmark_path = '../benchmarks/'
labelsize = 14
legsize = 14
xlabelpad = 20
ticksize = 12
linec = 'xkcd:greyish green' # 'xkcd:pale turquoise'  # 'xkcd:british racing green'  # '#d88868'

names_mc = ['Ea', 'eta_pre',
            # 'x_Eu'
            ]
mini_mc = [200e3, 2.63e10, #240e3, 1.6e10
           # 0.3
           ]
maxi_mc = [300e3, 5.32e13, #340e3, 2.6e12,
           # 3
           ]

yticks = None
xlabel = 'Time (Gyr)'

# good version
yvars = ['T_m',
         # 'h_rad_m',
         'eta_m',
         'b',
         'Ra_i',
         'h_dim_factor',
         'dyn_top_rms',
         # 'dyn_top_rms_1param',
         # 'urey'
         ]
ylabels = [r'$T_m$ (K)',
           # r'$h_{rad}$ (pW/kg)',
           r'$\eta_m$ (Pa s)',
           r'$b$',
           r'Ra$_i$',
           r'$d \Delta T \alpha_m$' + '\n' + r'($\times 10^5$ m)',
           r'$h_{\rm rms}$ (m)', # + '\n' + r'$f$(Ra$_i$, $b$)',
           # r'$h_{\rm rms}$ (m)' + '\n' + r'$f$(Ra$_{i, {\rm eff}}$)',
           # 'Ur'
           ]
ylims = [(1600, 3500),
         # None,
         (1e17, 1e20),
         (8, 16),
         (1e7, 1e13),
         (1e-1, 4e0),
         (1e0, 1e3),
         # (1e2, 1e3),
         # None
         ]
log = [False,
       # False,
       True,
       False,
       True,
       False,
       True,
       # True,
       # False
       ]
yscales = [1,
           # 1e12,
           1,
           1,
           1,
           1e-5,
           1,
           # 1,
           # 1
           ]

num_dist = 100
line_alpha = 0.08 #0.11
masses = np.array([0.1, 1, 5]) * p.M_E
verbose = False

# # testing version
# yvars = ['T_m',
#          'dT_rh',
#          'eta_m',
#          'Ra_i'
#          ]
# ylabels = yvars
# ylims = [None]*len(yvars)
# log = [False, False, True, True]
# yscales = [1, 1, 1, 1]
# num_dist = 10
# masses = np.array([0.1, 1, 5]) * p.M_E
# verbose = True
# labelsize = 12

planet_kwargs = {
    'beta_u': 0.333
}
run_kwargs = {
    # 'rms_type': 'Ra_i_eff'
}

fig, axes = plt.subplots(len(yvars), len(masses), figsize=(5, 9.5))
if len(np.shape(axes)) == 1:
    axes = axes[:, np.newaxis]
for ii, mass in enumerate(masses):
    if ii > 0:
        ylabels = None
        # yticks = [[]] * len(yvars)
    print('\n', mass / p.M_E, 'M_E')
    planet_kwargs.update({'M_p': mass})
    if mass >= (1 * p.M_E):
        legtext = str(int(mass / p.M_E)) + r' $M_{\oplus}$'
    else:
        legtext = str(mass / p.M_E) + r' $M_{\oplus}$'
    fig, ax_col = plottop.plot_distribution(yvars, default='baseline', tf=4.5, propagate_fit_err=False,
                                            update_kwargs=planet_kwargs, run_kwargs=run_kwargs,
                                            num=num_dist, names=names_mc, mini=mini_mc, maxi=maxi_mc,
                                            xlabelpad=None, ylabelpad=10, n_sigma=1, ylims=ylims, tickpad=5,
                                            fig=fig, axes=axes[:, ii], c='k', lw=0.5, alpha=line_alpha, c_mean=linec,
                                            xticks=None, yticks=yticks, log=log,
                                            xlabel='', ylabels=ylabels, yscales=yscales, labelsize=labelsize,
                                            ticksize=ticksize, legsize=legsize, save=False,
                                            fig_path=fig_path, legtext=legtext, verbose=verbose)

    # # highlight baseline case
    # fig, ax_col = plottop.plot_distribution(yvars, default='baseline',
    #                                         update_kwargs=planet_kwargs, run_kwargs=run_kwargs,
    #                                         num=1, names=['Ea', 'eta_pre'], mini=[300e3, 1.6e11], maxi=[300e3, 1.6e11],
    #                                         xlabelpad=None, ylabelpad=10, n_sigma=0, ylims=ylims, tickpad=5,
    #                                         fig=fig, axes=axes[:, ii], c='k', lw=0.5, alpha=0.05, c_mean='r',
    #                                         xticks=None, yticks=yticks, log=log,
    #                                         xlabel='', ylabels=ylabels, yscales=yscales, labelsize=labelsize,
    #                                         ticksize=ticksize, legsize=legsize, save=False,
    #                                         fig_path=fig_path, legtext='',
    #                                         verbose=verbose)
    for ax in ax_col:
        # postprocessing matplotlib
        if ii > 0:
            tk = [item.get_text() for item in ax.get_yticklabels()]
            empty_string_labels = [''] * len(tk)
            ax.set_yticklabels(empty_string_labels)

# plt.minorticks_off()  # TODO
fig.supxlabel(xlabel, fontsize=labelsize, y=0.02)
plt.tight_layout()
plt.subplots_adjust(wspace=0.1)
plt.savefig(fig_path + 'evol_dist' + today + '.pdf', bbox_inches='tight')
# plt.show()
