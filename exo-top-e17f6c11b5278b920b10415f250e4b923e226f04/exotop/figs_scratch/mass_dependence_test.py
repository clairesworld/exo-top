from model_1D import the_results as plottop
from model_1D import parameters as p
import matplotlib.pyplot as plt

fig_path = ''  # laptop
fig_format = '.png'
labelsize = 14
legsize = 12
xlabelpad = 20
ticksize = 12
linec = 'xkcd:ocean green'  # 'xkcd:british racing green'  # '#d88868'

# how h varies across mass
x_vars = ['M_p', 'CMF', 'H_f']
units = ['$M_E$', 'CMF', 'pW kg$^{-1}$']
log = [False, False, False]
x_range = [(1 * p.M_E, 2 * p.M_E),  (0.2, 0.4), (4e-12, 8e-12)]
xticks = None # [[0.1, 1, 5], [0.1, 0.2, 0.3, 0.4], [2, 5, 10]]
xscales = [p.M_E ** -1, 1, 1e12]
xlabels = ['$M_p$', 'CMF', '$H_f$\n(pW kg$^{-1}$)']

y_vars = ['dT_m', 'd_m', 'eta_m', 'Ra_i_eff', 'dyn_top_rms']
ylabels = [r'$\Delta T_m$', r'$d$', r'$\eta_m$', r'Ra$_{i, eff}$', r'h$_{rms}$']
ylims = [(200, 350), (1.5e6, 4.5e6), (3e18, 2e20), (1e4, 1e9), (50, 900)]

fig, axes = plt.subplots(len(y_vars), len(x_vars))
for yy, yvar in enumerate(y_vars):
    print('\n---------', yvar)
    if yy < len(y_vars) - 1:
        xlabels1 = ['']*len(y_vars)
    else:
        xlabels1 = xlabels
    fig, _ = plottop.plot_change_with_observeables_ensemble(age=4.5, dist_res=500, x_res=5, verbose=True,
                                                               defaults='baseline', fig=fig, axes=axes[yy, :],
                                                               ticksize=ticksize, labelsize=labelsize, fig_height=6,
                                                               legend=False, lw=4, ylabel=ylabels[yy],
                                                               labelpad=20, legendtop=True, tickwidth=2,
                                                               save=False,  fig_path=fig_path,
                                                               update_kwargs={'visc_type': 'KW'},
                                                               model_param=yvar, labels=[''],
                                                               x_vars=x_vars, units=units, log=log, x_range=x_range,
                                                               xscales=xscales, xlabels=xlabels1, ylim=ylims[yy],
                                                               linec=linec, textc='k',  alpha=0.3, legsize=legsize)

plt.show()
