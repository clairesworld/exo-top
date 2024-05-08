import matplotlib.pyplot as plt
import numpy as np
import model_1D.parameters as p
import model_1D.the_results as plottop

fig_format = '.png'
fig_path = '/home/claire/Desktop/'

def evols_with_Mp(pl_baseline='Earthbaseline', v_eval=[0.1, 1, 2, 3, 4, 5, 6], save=True, fig_path=fig_path,
                  names=None, log=False, nrows=None, ncols=None, pl_update_args={}, model_update_args={}, age=4.5,
                  xlim=None, vscale=p.M_E, vunits='M$_E$', vname='M_p', **kwargs):
    if names is None:
        names = {'Ra_i': ('log(Ra)', 1),
                 'dyn_top_rms': ('$\Delta h_{rms}$ (m)', 1),
                 }
    if nrows is None:
        nrows = 1
    if ncols is None:
        ncols = len(names)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))
    # col = ['#4cc470','#68b894','#a4a0c0', '#e084f0', '#a86494', '#5c5034', '#506080', '#3868b8']
    col = ['#4cc470', '#7cb0a0', '#a4a0c0', '#c890dc', '#ec80f8', '#b468a4', '#7c5454']
    # col_rgb = [(tuple(int(h.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))) for h in col]
    col_rgb = [[76, 196, 112], [124, 176, 160], [164, 160, 192], [200, 144, 220], [236, 128, 248],
               [180, 104, 164], [124, 84, 84]]

    # print(col_rgb)
    # cm = mpl.colors.LinearSegmentedColormap.from_list('test',np.array(col_rgb)/255.0)

    if xlim is None:
        xlim = (0, age)
    model_update_args.update({'tf': age})
    for ii, v in enumerate(reversed(v_eval)):
        pl_update_args.update({vname: v * vscale})

        fig, axes = plottop.benchmark_thermal_plots(ident='Earthbaseline', names=names,
                                                    pl_update_args=pl_update_args,
                                                    model_update_args=model_update_args,
                                                    # compare_dir=['benchmarks/thiriet_Mars1', 'benchmarks/breuer_Mars'],
                                                    fig_path=fig_path, fig=fig, axes=axes,
                                                    plots_save=False, print_tf=True,
                                                    legsize=14, ncols=ncols, labelsize=34,
                                                    ticksize=25,
                                                    show_qsfc_error=False, title='', legax=None,
                                                    fontname='Ubuntu Mono',
                                                    label=str(v) + ' ' + vunits,
                                                    line_args={'c': col[ii], 'lw': 4, 'ls': '-'},
                                                    xlim=xlim, **kwargs
                                                    )

    legend = axes.flatten()[-1].legend(frameon=False, fontsize=25, loc='upper left',
                                       bbox_to_anchor=(1.05, 0.95),
                                       borderaxespad=0, ncol=1)
    # sc = axes[0].scatter([0,1], [0,1], cmap=cm)
    # fig.subplots_adjust(right=0.8)
    # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    # fig.colorbar(sc, cax=cbar_ax)
    # sc.set_visible=False

    for ax in axes.flatten():
        ax.xaxis.set_tick_params(width=2)
        ax.yaxis.set_tick_params(width=2)
        ax.set_facecolor('#f7f7f7')
        if log:
            ax.set_yscale('log')

    fig.subplots_adjust(wspace=.49)
    if save:
        fig.savefig(fig_path + 'therm2' + fig_format, bbox_inches='tight', dpi=200)
    return fig, axes


# evols_with_Mp(pl_baseline='Earthbaseline', save=True)


names = {'Ra_i':('Ra_i',1),
                  'H_rad_m':('$H_rad$ (TW)',1e-12),
                  'h_rad_m':('$h_rad$ (mW/m2)',1e3),
         # 'dyn_top_rms':('$\Delta h_{rms}$ (m)',1),
         'urey':('Ur',1),
#          'dyn_top_heuristic':('$\Delta h_{rms}$ heuristic (m)',1),
        'q_sfc':('$q_{sfc}$ (mW/m2)',1e3),
         'T_c':('$T_c$ (K)',1),
         'Q_sfc':('$Q_{sfc}$ (TW)',1e-12),
            'T_m':('$T_m$ (K)',1),
        }

pl_update_args = {}
# model_update_args = {'T_m0':1600, 'T_c0':2200}
# model_update_args = {'age':15}
model_update_args = {}

fig, axes = evols_with_Mp(pl_baseline='Earthbaseline', save=True, names=names, log=False, nrows=3, ncols=3,
             pl_update_args=pl_update_args, model_update_args=model_update_args, xlim=(1, 4.5))

plt.show()