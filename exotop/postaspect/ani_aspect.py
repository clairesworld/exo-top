import sys
import os
# sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    fig_fmt, regime_grid_td, postprocess_kwargs, regime_names_td, \
    load_grid    # noqa: E402
from postaspect import plt_aspect as sc  # noqa: E402
from postaspect import aspectdata as ap
from postaspect.aspect_post import pickleio
from useful_and_bespoke import hide_log_ticklabels, not_iterable, dark_background
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
from matplotlib import cm, colors, rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd

def animate_T_field(case, data_path=data_path, fig_path=fig_path, labelsize=30, ticksize=16, cmap='gist_heat', fps=15,
                    shading='nearest', c='k', dark=True):
    fig, ax, im, n, T_n = static_T_field(case, data_path=data_path, fig_path=fig_path, labelsize=labelsize, ticksize=ticksize,
                                     cmap=cmap, shading=shading, return_artists=True, save=False, c=c, dark=dark)

    # # set frame rate according to model time:
    # if fps is None:
    #     t = df.time.to_numpy()
    #     dt = t[-1] - t[0]
    #     print('model time span:', dt)
    #     factor = 1 / 68
    #     fps = len(n) / (dt) * factor
    #     print('equivalent fps:', fps)

    def init():
        im.set_array([])
        return im

    # do animation
    def animate(i, T_n):
        T = ap.reduce_dims(T_n[i])
        im.set_array(T.ravel())
        return im

    try:
        ani = animation.FuncAnimation(fig, animate, frames=len(n),
                                                                    init_func=init,
                                      fargs=(T_n,), blit=False, repeat=True,
                                      interval=200, )  # interval: delay between frames in ms
        # HTML(ani.to_jshtml())
        ani.save(fig_path+case + '-T.gif', writer='imagemagick', fps=fps, savefig_kwargs={'facecolor': fig.get_facecolor()})
        # print('Finished!!')
    except:
        # out of cache?
        ani = animation.FuncAnimation(fig, animate, frames=range(0, len(n), 2),
                                                                    init_func=init,
                                      fargs=(T_n,), blit=False, repeat=True,
                                      interval=200, )  # interval: delay between frames in ms
        ani.save(fig_path+case + '-T.gif', writer='imagemagick', fps=fps, savefig_kwargs={'facecolor': fig.get_facecolor()})


def animate_T_prof(case, data_path=data_path, fig_path=fig_path, labelsize=30, ticksize=16, fps=15, dark=True):

    fig2, ax, line, df = static_T_prof(case, data_path=data_path, fig_path=fig_path, labelsize=labelsize,
                                       ticksize=ticksize, dark=dark, return_artists=True, save=False)

    def init2():
        line.set_xdata(([np.nan] * len(np.array(df.iloc[0]['y'].tolist()))))
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

    ani2 = animation.FuncAnimation(fig2, animate2, frames=len(df),
                                  #                               init_func=init2,
                                  fargs=(df,), blit=False, repeat=True,
                                  interval=200, )  # interval: delay between frames in ms
    # HTML(ani.to_jshtml())
    ani2.save(fig_path+case + '-prof.gif', writer='imagemagick', fps=fps, savefig_kwargs={'facecolor': fig2.get_facecolor()})


def animate_h(case, data_path=data_path, fig_path=fig_path, labelsize=30, ticksize=16, fps=15, dark=True):

    fig, ax, hprof, hmean, h_n, h_rms = static_h(case, data_path=data_path, fig_path=fig_path, labelsize=labelsize,
                                     ticksize=ticksize, dark=dark, return_artists=True, save=False)

    def init():
        hprof.set_ydata(([np.nan] * len(h_n[0])))
        hmean.set_ydata(([np.nan] * len(h_n[0])))
        return hprof, hmean,

    def animate(i, h_n, h_rms):
        hprof.set_ydata((h_n[i]))  # update the data.
        hmean.set_ydata(([h_rms[i]]*len(h_n[i])))
        return hprof, hmean,

    ani = animation.FuncAnimation(fig, animate, frames=len(h_rms), init_func=init,
                                  fargs=(h_n,h_rms), blit=True, repeat=True,
                                  interval=200, )  # interval: delay between frames in ms
    ani.save(fig_path+case + '-h.gif', writer='imagemagick', fps=fps, savefig_kwargs={'facecolor': fig.get_facecolor()})


def static_h(case, data_path=data_path, fig_path=fig_path, labelsize=30, ticksize=16, dark=False,
             return_artists=False, c='k', save=True, i_ts=0, avg=False, fig=None, ax=None, **kwargs):
    if dark:
        foreground = 'xkcd:off white'
    else:
        foreground = c

    # preload data
    df = pickleio(case=case, load=True, suffix='_T', postprocess_functions=None, data_path=data_path)
    n = df.sol.to_numpy()
    ts = df.index.to_numpy()

    h_n, h_peak, h_rms = [], [], []
    for i in range(len(n)):
        x, h = sc.read_topo_stats(case, ts[i], data_path=data_path)
        h_norm = sc.trapznorm(h)  # normalize to 0 mean
        peak, rms = sc.peak_and_rms(h_norm)
        h_n.append(h_norm)
        h_peak.append(peak)
        h_rms.append(rms)
    print('loaded', len(n), 'h profiles')

    if fig is None and ax is None:
        fig, ax = plt.subplots(figsize=(20, 1))
    ax.set_xlabel('', fontsize=labelsize, labelpad=20)
    ax.set_ylabel('', fontsize=labelsize, labelpad=20)
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.set_xticks([])
    # ax.set_yticks([])

    hprof, = ax.plot(x, h_n[i_ts], c=foreground, lw=3)
    hmean, = ax.plot(x, [h_rms[i_ts]] * len(x), c=foreground, lw=2, ls='--')
    ax.set_xlim([0, 8])
    ax.set_ylim([-5e-2, 5e-2])
    if dark:
        fig, ax = dark_background(fig, ax)
    if save:
        sc.plot_save(fig, case+'_h_prof', fig_path=fig_path)

    if return_artists:
        return fig, ax, hprof, hmean, h_n, h_rms
    else:
        return fig, ax


def static_T_prof(case, data_path=data_path, fig_path=fig_path, labelsize=30, ticksize=16, legsize=20, dark=False,
                  return_artists=False, c='k', alpha=0.9, save=True, i_n=0, avg=False, fig=None, ax=None, **kwargs):

    if dark:
        foreground = 'xkcd:off white'
    else:
        foreground = c

    df = pickleio(case=case, load=True, suffix='_T', postprocess_functions=None, data_path=data_path)
    n = df.sol.to_numpy()

    # T profile
    if fig is None and ax is None:
        fig, ax = plt.subplots(figsize=(2, 4))
    ax.set_xlabel('', fontsize=labelsize, labelpad=20)
    ax.set_ylabel('', fontsize=labelsize, labelpad=20)
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.set_xticks([0, 0.5, 1])
    ax.set_yticks([])

    if avg:
        # load time-averages
        T_f, y_f = ap.time_averaged_profile_from_df(df, 'T_av')
        uv_mag_av, y = ap.time_averaged_profile_from_df(df, 'uv_mag_av')
        dic_av = ap.T_parameters_at_sol(case, n=None, T_av=T_f, uv_mag_av=uv_mag_av, y=y_f,
                                         **postprocess_kwargs, **kwargs)  # actually a dict
        for k in ['T_av', 'uv_mag_av', 'y']:
            dic_av.pop(k, None)
        df_n = pd.DataFrame({key: value for (key, value) in dic_av.items()}, index=[0])
    else:
        df_n = df.iloc[i_n]
        T_f = np.array(df_n['T_av'].tolist())
        y_f = np.array(df_n['y'].tolist())

    delta_rh_n = np.array(df_n['delta_rh'])
    D_l_n = np.array(df_n['y_L'])

    line, = ax.plot(T_f, y_f, c=foreground, lw=3)
    ax.fill_between(T_f, [D_l_n - delta_rh_n] * len(T_f), [D_l_n] * len(T_f), fc='xkcd:tangerine', alpha=alpha,
                    label='Thermal bdy layer')
    #                 label='Thermal bdy layer\n' + r'$\delta$ = ' + '{:04.3f}'.format(delta_rh_n))
    ax.legend(frameon=False, fontsize=legsize, ncol=1, bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])

    if dark:
        fig, ax = dark_background(fig, ax)
    if save:
        sc.plot_save(fig, case + '_T_prof', fig_path=fig_path)

    if return_artists:
        return fig, ax, line, df
    else:
        return fig, ax


def static_T_field(case, data_path=data_path, fig_path=fig_path, labelsize=30, ticksize=16, cmap='gist_heat',
                   shading='nearest', return_artists=False, save=True, i_n=0, avg=False, c='k', dark=False,
                   fig=None, ax=None):
    if dark:
        foreground = 'xkcd:off white'
    else:
        foreground = c

    # preload data
    dat = ap.Aspect_Data(directory=data_path + 'output-' + case + '/')
    df = pickleio(case=case, load=True, suffix='_T', postprocess_functions=None, data_path=data_path)
    n = df.sol.to_numpy()

    if len(n) > 100:
        n = n[::2]

    T_n = []
    for nn in n:
        x, y, _, T = dat.read_temperature(nn, verbose=False)
        T_n.append(T)
    print('loaded', len(n), 'T fields')

    if avg:
        T_n = np.array(T_n)
        print('T_n', np.shape(T_n))
        T_im = np.mean(T_n, axis=0)
        print('T_im', np.shape(T_im))
    else:
        T_im = T_n[i_n]

    # set up initial image
    x = dat.x
    y = dat.y

    if ax is None and fig is None:
        fig, ax = plt.subplots(figsize=(20, 10))

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
    im = ax.pcolormesh(x, y, ap.reduce_dims(T_im), cmap=cmap, shading=shading)
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cax.xaxis.tick_top()
    cax.xaxis.set_label_position('top')
    cax.xaxis.set_ticks_position('top')

    if dark:
        fig, ax, cax = dark_background(fig, [ax, cax])
    if save:
        sc.plot_save(fig, case + '_T_field', fig_path=fig_path)

    if return_artists:
        return fig, ax, im, n, T_n
    else:
        return fig, ax


def T_h_gridspec(case, data_path=data_path, fig_path=fig_path, labelsize=30, ticksize=16, cmap='gist_heat',
                   shading='nearest', save=True, i_n=0, avg=False, c='k', dark=True):
    # not animated

    gs = gridspec.GridSpec(1, 9)
    ax = plt.subplot(gs[:-1])
    fig, ax = static_T_field(case, data_path=data_path, )