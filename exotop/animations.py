import sys
import os
sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    fig_fmt, regime_grid_td, postprocess_kwargs, regime_names_td, \
    load_grid    # noqa: E402
from exotop.postaspect import aspect_scalings as sc  # noqa: E402
from exotop.postaspect import aspect_postprocessing2 as ap
from exotop.useful_and_bespoke import hide_log_ticklabels, not_iterable, dark_background
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm, colors, rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd

def animate_T(case, data_path=data_path, fig_path=fig_path, labelsize=30, ticksize=16, cmap='coolwarm', fps=15):
    # preload data
    dat = ap.Aspect_Data(directory=data_path + 'output-' + case + '/')
    df = sc.pickleio(case=case, load=True, suffix='_T', postprocess_functions=None, data_path=data_path)
    n = df.sol.to_numpy()

    if len(n) > 100:
        n = n[::2]

    T_n = []
    for nn in n:
        x, y, _, T = dat.read_temperature(nn, verbose=False)
        T_n.append(T)
    print('loaded', len(n), 'T fields')

    # set up initial image
    x = dat.x
    y = dat.y

    fig, ax = plt.subplots(figsize=(20, 10))
    # gs = gridspec.GridSpec(1, 9)
    # ax = plt.subplot(gs[:-1])

    ax.set_xlabel('x', fontsize=labelsize, labelpad=20)
    ax.set_ylabel('y', fontsize=labelsize, labelpad=20)
    ax.set_xticks([])
    ax.set_yticks([])
    # ax.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.axis('equal')

    divider = make_axes_locatable(ax)
    cax = divider.new_vertical(size='10%', pad=0.3)
    fig.add_axes(cax)
    cax.set_xlabel('Temperature', fontsize=18)
    cax.tick_params(axis='x', which='major', labelsize=ticksize)
    ax.set_title('Nondimensional temperature', fontsize=labelsize, pad=90, color='xkcd:off white')
    im = ax.pcolormesh(x, y, ap.reduce_dims(T_n[0]), cmap=cmap)
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cax.xaxis.tick_top()
    cax.xaxis.set_label_position('top')
    cax.xaxis.set_ticks_position('top')
    fig, ax, cax = dark_background(fig, [ax, cax])

    # set frame rate according to model time:
    if fps is None:
        t = df.time.to_numpy()
        dt = t[-1] - t[0]
        print('model time span:', dt)
        factor = 1 / 68
        fps = len(n) / (dt) * factor
        print('equivalent fps:', fps)

    # do animation
    def init():
        im.set_array([])
        return im

    def animate(i, T_n):
        T = ap.reduce_dims(T_n[i])
        im.set_array(T.ravel())
        return im

    try:
        ani = animation.FuncAnimation(fig, animate, frames=len(n),
                                      #                               init_func=init,
                                      fargs=(T_n,), blit=False, repeat=True,
                                      interval=200, )  # interval: delay between frames in ms
        # HTML(ani.to_jshtml())
        ani.save(fig_path+case + '-T.gif', writer='imagemagick', fps=fps, savefig_kwargs={'facecolor': fig.get_facecolor()})
        # print('Finished!!')
    except:
        # out of cache?
        ani = animation.FuncAnimation(fig, animate, frames=range(0, len(n), 2),
                                      #                               init_func=init,
                                      fargs=(T_n,), blit=False, repeat=True,
                                      interval=200, )  # interval: delay between frames in ms
        ani.save(fig_path+case + '-T.gif', writer='imagemagick', fps=fps, savefig_kwargs={'facecolor': fig.get_facecolor()})



    # T profile
    fig2, ax = plt.subplots(figsize=(2, 4))
    ax.set_xlabel('', fontsize=labelsize, labelpad=20)
    ax.set_ylabel('', fontsize=labelsize, labelpad=20)
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.set_xticks([0, 0.5, 1])
    ax.set_yticks([])

    df_n = df.iloc[0]

    delta_rh_n = np.array(df_n['delta_rh'])
    D_l_n = np.array(df_n['y_L'])
    T_f = np.array(df_n['T_av'].tolist())
    y_f = np.array(df_n['y'].tolist())

    line, = ax.plot(T_f, y_f, c='xkcd:off white', lw=3)
    ax.fill_between(T_f, [D_l_n - delta_rh_n] * len(T_f), [D_l_n] * len(T_f), fc='xkcd:tangerine', alpha=0.9,
                    label='Thermal bdy layer')
    #                 label='Thermal bdy layer\n' + r'$\delta$ = ' + '{:04.3f}'.format(delta_rh_n))
    ax.legend(frameon=False, fontsize=20, ncol=1, bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    fig2, ax = dark_background(fig2, ax)

    def init2():
        line.set_xdata(([np.nan] * len(y_f)))
        return line,

    def animate2(i, df):
        ax.collections.clear()
        df_n = df.iloc[i]
        T_f = np.array(df_n['T_av'].tolist())
        line.set_xdata(T_f)  # update the data.

        D_l_n = np.array(df_n['y_L'])
        delta_rh_n = np.array(df_n['delta_rh'])
        ax.fill_between(T_f, [D_l_n - delta_rh_n] * len(T_f), [D_l_n] * len(T_f), fc='xkcd:tangerine', alpha=0.9,
                        label='Thermal bdy layer')
        return line,

    ani2 = animation.FuncAnimation(fig2, animate2, frames=len(n),
                                  #                               init_func=init2,
                                  fargs=(df,), blit=False, repeat=True,
                                  interval=200, )  # interval: delay between frames in ms
    # HTML(ani.to_jshtml())
    ani2.save(fig_path+case + '-prof.gif', writer='imagemagick', fps=fps, savefig_kwargs={'facecolor': fig.get_facecolor()})


def animate_h(case, data_path=data_path, fig_path=fig_path, labelsize=30, ticksize=16, fps=15):
    # preload data
    df = sc.pickleio(case=case, load=True, suffix='_T', postprocess_functions=None, data_path=data_path)
    n = df.sol.to_numpy()
    ts = df.index.to_numpy()

    h_n, h_peak, h_rms = [], [], []
    for i in range(len(n)):
        x, h = sc.read_topo_stats(case, ts[i], data_path=data_path)
        h_norm = sc.trapznorm(h)
        peak, rms = sc.peak_and_rms(h_norm)
        h_n.append(h)
        h_peak.append(peak)
        h_rms.append(rms)
    print('loaded', len(n), 'h profiles')

    fig, ax = plt.subplots(figsize=(20, 3.5))
    ax.set_xlabel('', fontsize=labelsize, labelpad=20)
    ax.set_ylabel('', fontsize=labelsize, labelpad=20)
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.set_xticks([])
    # ax.set_yticks([])

    hprof, = ax.plot(x, h_n[0], c='xkcd:off white', lw=3)
    hmean, = ax.plot(x, [h_rms[0]]*len(x), c='xkcd:off white', lw=2, ls='--')
    ax.set_xlim([0, 8])
    ax.set_ylim([-5e-2, 5e-2])
    fig, ax = dark_background(fig, ax)

    def init():
        hprof.set_ydata(([np.nan] * len(x)))
        hmean.set_ydata(([np.nan] * len(x)))
        return hprof, hmean,

    def animate(i, h_n, h_rms):
        hprof.set_ydata((h_n[i]))  # update the data.
        hmean.set_ydata(([h_rms[i]]*len(x)))
        return hprof, hmean,

    ani = animation.FuncAnimation(fig, animate, frames=len(n), init_func=init,
                                  fargs=(h_n,h_rms), blit=True, repeat=True,
                                  interval=200, )  # interval: delay between frames in ms
    ani.save(fig_path+case + '-h.gif', writer='imagemagick', fps=fps, savefig_kwargs={'facecolor': fig.get_facecolor()})

for jj, etastr in enumerate(eta_ls):
    if jj <= 20:
        cases, cases_var = sc.get_cases_list(Ra_ls, etastr, end_grid[jj])
        for ii, case in enumerate(cases):
            if (os.path.exists(data_path + 'output-' + case)) and (ii >= 0):
                # animate_T(case, data_path=data_path, fig_path=fig_path, labelsize=30, ticksize=16, cmap='coolwarm')
                animate_h(case, data_path=data_path, fig_path=fig_path, labelsize=30, ticksize=16)
                print('finished case')
