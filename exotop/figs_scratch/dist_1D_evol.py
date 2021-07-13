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
labelsize = 20  # 16
legsize = 16
xlabelpad = 20
ticksize = 12
linec = 'xkcd:pale turquoise'  # 'xkcd:british racing green'  # '#d88868'

names_mc = ['Ea', 'eta_pre']
mini_mc = [240e3, 1.5e10]
maxi_mc = [300e3, 2.5e12]

yticks = None
xlabel = 'Time (Gyr)'

# good version
# yvars = ['h_rad_m', 'T_m', 'eta_m', 'Ra_i_eff', 'h_dim_factor', 'dyn_top_rms']
# ylabels = [r'$h_{rad}$ (pW/kg)', r'$T_m$ (K)', r'$\eta_m$ (Pa s)', r'Ra$_{i, {\rm eff}}$',
#            r'$d_m \Delta T_m \alpha_m$ (km)', r'$h_{\rm rms}$ (km)']
# ylims = [None, (1450, 2250), (1e18, 1e23), (1e4, 1e12), (3e1, 3e2), (1e-1, 3e1)]
# log = [False, False, True, True, True, True]
# yscales = [1e12, 1, 1, 1, 1e-3, 1e-3]

# testing version
yvars = ['T_m', 'T_c', 'D_l']
ylabels = ['$T_m$ (K)', '$T_c$ (K)', '$D_l$ (km)']
ylims = [None, None, None]
log = [False, False, False]
yscales = [1, 1, 1e-3]

# masses = np.array([0.1, 1, 3, 5])*p.M_E
masses = np.array([1]) * p.M_E

num_dist = 5
fig, axes = plt.subplots(len(yvars), len(masses), figsize=(7, 9.5))
if len(np.shape(axes)) == 1:
    axes = axes[:, np.newaxis]
for ii, mass in enumerate(masses):
    if ii > 0:
        ylabels = None
        yticks = [[]] * len(yvars)
    print('\n', mass / p.M_E, 'M_E')
    fig, _ = plottop.plot_distribution(yvars, default='baseline', update_kwargs={'M_p': mass},
                                       num=num_dist, names=names_mc, mini=mini_mc, maxi=maxi_mc,
                                       xlabelpad=None, ylabelpad=10, n_sigma=1, ylims=ylims, tickpad=5,
                                       fig=fig, axes=axes[:, ii], c='k', lw=0.5, alpha=0.05, c_mean=linec,
                                       xticks=None, yticks=yticks, log=log,
                                       xlabel='', ylabels=ylabels, yscales=yscales, labelsize=labelsize,
                                       ticksize=ticksize, legsize=legsize, save=False,
                                       fname='evol_dist', fig_path=fig_path, legtext=str(mass / p.M_E) + r' $M_E$',
                                       verbose=True)

# plt.minorticks_off()  # TODO
fig.supxlabel(xlabel, fontsize=labelsize, y=0.02)
plt.tight_layout()
plt.subplots_adjust(wspace=0.1)

plt.savefig('evol_dist.png', bbox_inches='tight')
plt.show()
