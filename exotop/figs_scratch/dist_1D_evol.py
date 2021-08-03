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
legsize = 16
xlabelpad = 20
ticksize = 12
linec = 'xkcd:pale turquoise'  # 'xkcd:british racing green'  # '#d88868'

names_mc = ['Ea', 'eta_pre',
            # 'x_Eu'
            ]
mini_mc = [240e3, 1.6e10,
           # 0.3
           ]
maxi_mc = [300e3, 2.6e12,
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
         'dyn_top_rms'
         ]
ylabels = [r'$T_m$ (K)',
           # r'$h_{rad}$ (pW/kg)',
           r'$\eta_m$ (Pa s)',
           r'$b$',
           r'Ra$_i$',
           r'$d_m \Delta T_m \alpha_m$'+'\n'+r'($\times 10^5$ m)',
           r'$h_{\rm rms}$ (m)']
ylims = [(1600, 3000),
         # None,
         (1e17, 1e20),
         (8, 20),
         (1e7, 1e13),
         (1e-1, 5e0),
         (1e0, 1e3)]
log = [False,
       # False,
       True,
       False,
       True,
       False,
       True]
yscales = [1,
           # 1e12,
           1,
           1,
           1,
           1e-5,
           1]

num_dist = 1000
masses = np.array([0.1, 1, 5])*p.M_E
verbose = False
labelsize = 16

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

fig, axes = plt.subplots(len(yvars), len(masses), figsize=(5, 9.5))
if len(np.shape(axes)) == 1:
    axes = axes[:, np.newaxis]
for ii, mass in enumerate(masses):
    if ii > 0:
        ylabels = None
        # yticks = [[]] * len(yvars)
    print('\n', mass / p.M_E, 'M_E')
    fig, ax_col = plottop.plot_distribution(yvars, default='baseline', update_kwargs={'M_p': mass},
                                       num=num_dist, names=names_mc, mini=mini_mc, maxi=maxi_mc,
                                       xlabelpad=None, ylabelpad=10, n_sigma=1, ylims=ylims, tickpad=5,
                                       fig=fig, axes=axes[:, ii], c='k', lw=0.5, alpha=0.05, c_mean=linec,
                                       xticks=None, yticks=yticks, log=log,
                                       xlabel='', ylabels=ylabels, yscales=yscales, labelsize=labelsize,
                                       ticksize=ticksize, legsize=legsize, save=False,
                                       fig_path=fig_path, legtext=str(mass / p.M_E) + r' $M_E$',
                                       verbose=verbose)
    for ax in ax_col:
        # postprocessing
        if ii > 0:
            tk = [item.get_text() for item in ax.get_yticklabels()]
            empty_string_labels = [''] * len(tk)
            ax.set_yticklabels(empty_string_labels)

# plt.minorticks_off()  # TODO
fig.supxlabel(xlabel, fontsize=labelsize, y=0.02)
plt.tight_layout()
plt.subplots_adjust(wspace=0.1)
plt.savefig('evol_dist.png', bbox_inches='tight')
plt.show()
