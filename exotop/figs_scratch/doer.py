from model_1D import the_results as plt1
import numpy as np
from model_1D.parameters import M_E

fig_path = ''
labelsize = 16
ticksize = 12
legsize = 16

yvars = ['T_m', 'T_c', 'delta_rh', 'd_m', 'eta_m', 'Ra_i_eff']
ylabels = [r'$T_m$ (K)', r'$T_c$ (K)', r'$\delta_u$ (km)', r'$d_m$ (km)', r'$\eta_m$ (Pa s)', r'Ra$_{i, eff}$']
yscales = [1, 1, 1e-3, 1e-3, 1, 1]
log = [False, False, False, False, True, True]
ylims = [(1600, 2250), (1800, 2400), (0, 60), None, (1e17, 1e21), (0.8e8, 1e12)]
run_kwargs = {}

update_kwargs = {'M_p': 1*M_E}
fig, ax = plt1.plot_distribution(yvars, default='baseline',
                      num=200, update_kwargs=update_kwargs, run_kwargs=run_kwargs,
                      names=['Ea', 'eta_pre', 'T_m0', 'T_c0', 'D_l0'],
                      mini=[240e3, 1.5e10, 1600, 2000, 100e3],
                      maxi=[300e3, 2.5e12, 2000, 2500, 300e3],
                      xlabelpad=None, ylabelpad=None, n_sigma=1, ylims=ylims,
                      fig=None, ax=None, c='k', lw=0.5, alpha=0.2, c_mean='xkcd:sea green', log=log, yticks=None,
                      xlabel='Time (Gyr)', ylabels=ylabels, yscales=yscales, labelsize=labelsize, ticksize=ticksize,
                                 legsize=legsize, save=True,
                      fname='evol_dist_1ME', fig_path='', legtext=r'1 $M_E$', verbose=False)

update_kwargs = {'M_p': 5*M_E}
fig, ax = plt1.plot_distribution(yvars, default='baseline',
                      num=200, update_kwargs=update_kwargs, run_kwargs=run_kwargs,
                      names=['Ea', 'T_m0', 'T_c0', 'D_l0'],
                      mini=[240e3, 1600, 2000, 100e3],
                      maxi=[300e3, 2000, 2500, 300e3],
                      xlabelpad=None, ylabelpad=None, n_sigma=1, ylims=ylims,
                      fig=None, ax=None, c='k', lw=0.5, alpha=0.2, c_mean='xkcd:sea green', log=log, yticks=None,
                      xlabel='Time (Gyr)', ylabels=ylabels, yscales=yscales, labelsize=labelsize, ticksize=ticksize,
                                 legsize=legsize, save=True,
                      fname='evol_dist_5ME', fig_path='', legtext=r'5 $M_E$', verbose=False)
