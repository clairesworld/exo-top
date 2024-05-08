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
from matplotlib import rc
from matplotlib.pyplot import rcParams
from matplotlib import use as mpluse
import matplotlib.collections

mpluse('QTAgg')
#mpluse('Agg')
rc('text', usetex=True)  # turn off for running over ssh
# rcParams['font.family'] = 'serif'
# rcParams['font.serif'] = 'CMU Serif'
today = date.today().strftime("%b-%d-%Y")

# set paths
fig_path='/home/claire/Desktop/rw-poster/'
fig_format = '.png'
benchmark_path = '../benchmarks/'
# planet_icon_path = '/home/claire/Pictures/science-graphics/planet_png/'
labelsize = 49
legsize = 26  # 30
xlabelpad = 15
ylabelpad = 15
ticksize = 30
linec = 'xkcd:greyish green'  #'xkcd:ocean green'  # 'xkcd:british racing green'  # '#d88868'
linec2 = None  #'xkcd:pale turquoise'
alpha = 0.3  # 0.3
lw = 3  # 4

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
x_vars = [
          'M_p'
          ]  # 'H_f'
units = [

         'Planet mass (Earth masses)']  # pW kg$^{-1}$
log = [
       True]
x_range = [
           (0.1 * p.M_E, 5 * p.M_E)]  # (2e-12, 10e-12)
xticks = [
          ]  # [2, 5, 10]
xscales = [
           p.M_E ** -1]  # 1e12
xlabels = [
           'Planet mass\n' + r'(Earth masses)']  # 'Radiogenic heating\n(pW kg$^{-1}$)'

dist_res = 300  # works with 10000 - production run
x_res = 7  # works with 7 - production run
n_sigma = 1

# fig, axes = plottop.plot_change_with_observeables_ensemble(dist_res=dist_res, x_res=x_res, n_sigma=n_sigma,
#                                                            names=names, mini=mini_mc, maxi=maxi_mc,
#                                                            defaults='baseline', age=4.5,
#                                                            ticksize=ticksize, labelsize=labelsize, fig_height=6,
#                                                            legend=True, lw=lw, ylabel=r'RMS topography (m)',
#                                                            labelpad=20, legendtop=True, tickwidth=1,
#                                                            save=False,
#                                                            update_kwargs=planet_kwargs, run_kwargs=run_kwargs,
#                                                            model_param='dyn_top_rms', labels=[''],
#                                                            x_vars=x_vars, units=units, log=log, x_range=x_range,
#                                                            xscales=xscales, xlabels=xlabels,
#                                                            linec=linec,
#                                                            leg_loc='upper right',
#                                                            textc='k',  # 'xkcd:off white',
#                                                            alpha=alpha, legsize=legsize, xlabelpad=xlabelpad,
#                                                            ylabelpad=ylabelpad,
#                                                            extra_def=True,
#                                                            # verbose=True
#                                                            picklefrom=fig_path + 'h_parameters_fig.pkl',
#                                                            )
fig = pkl.load(open(fig_path + 'h_parameters_fig.pkl', "rb"))
axes = fig.get_axes()
axes[0].set_xlabel('Planet mass $(M_\oplus)$', labelpad=20, fontsize=labelsize)
axes[0].set_ylabel('RMS topography (m)', labelpad=18, fontsize=labelsize)
axes[0].tick_params(axis='both', labelsize=ticksize)

for i, ax in enumerate(axes):
    ax.set_xlim([x * xscales[i] for x in x_range[i]])
    ax.set_xticks([0.1, 1, 5])
    ax.set_yscale('log')
    # ax.set_ylim((10, 30000))
    # ax.set_yticks((200, 800, 1000))
    ax.set_ylim((10, 800))
    ax.set_yticks((10, 100, 800))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    print(ax.get_xlim())
    ax.tick_params(axis='both', which='major', pad=10)
axes[0].yaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
axes[0].yaxis.set_minor_formatter(ticker.NullFormatter())
# axes[0].set_yticks([600, 700, 800, 900, 1000, 1200,  1500])

axes[0].set_title('1D model', color='xkcd:off white', fontsize=labelsize)#

print(fig.texts)
print(axes[0].texts)
for txt in axes[0].texts:
    txt.set_visible(False)
axes[0].annotate(text='', xy=(2.01,31), xytext=(2.01,67), arrowprops=dict(arrowstyle='<->, head_width=1, head_length=1',
                                                                          lw=5, fc='xkcd:off white', ec='xkcd:off white',
                                                                          ))
axes[0].text(1.5,35.3,'rheological\nuncertainty', c='xkcd:off white', fontsize=35, ha='right', va='top')

fig, *axes = dark_background(fig, axes)
plt.subplots_adjust(wspace=0.1)
# plt.tight_layout()

# with open(fig_path + r"h_parameters_fig.pkl", "wb") as f:
#     pkl.dump(fig, f)

# change all spines
for axis in ['top', 'bottom', 'left', 'right']:
    axes[0].spines[axis].set_linewidth(5)

# increase tick width
axes[0].tick_params(width=5)
plt.gca().get_lines()[0].set_color("#35afbfff")
plt.gca().get_lines()[0].set_linewidth(5)
handles = plt.gca().get_children()  #plt.gca().get_legend_handles_labels()
for h in handles:
    print('h', h)
    if isinstance(h, matplotlib.collections.PolyCollection):
        h.set_color("#35afbf")

# plt.suptitle(r'Ra$_{i, {\rm eff}}$ scaling')
fig.savefig(fig_path + 'h_parameters.png', bbox_inches='tight', transparent=True, dpi=600)
plt.show()
