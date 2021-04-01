""" ASPECT runs: functions plotting """

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.image as mpimg
from matplotlib.colors import LogNorm, Normalize
import matplotlib.lines as mlines
from scipy.interpolate import interp1d
from postaspect.ani_aspect import static_T_field
from postaspect import aspectdata as ad
from postaspect import aspect_post as pro
from postaspect.setup_postprocessing import data_path_bullard, fig_path_bullard, highlight_colour, \
    cmap_path  # noqa: E402
from useful_and_bespoke import colorize, iterable_not_string, cmap_from_list, not_iterable, \
    colourbar, cmap_from_ascii, not_string, minmaxnorm, mahalanobis, colourised_legend
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})


# rc('text', usetex=True)

def plot_save(fig, fname, fig_path=fig_path_bullard, fig_fmt='.png', bbox_inches='tight', tight_layout=True, **kwargs):
    path = fig_path + fname + fig_fmt
    directory = os.path.dirname(path)
    os.makedirs(directory, exist_ok=True)
    if tight_layout:
        fig.tight_layout()
    fig.savefig(path, bbox_inches=bbox_inches, **kwargs)
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  saved to ', path, '!\n')


def plot_getx(Ra, eta, case=None, which_x=None, averagescheme='timefirst', data_path=data_path_bullard,
              t1=None, load=None, postprocess_kwargs=None, return_all=False, **kwargs):
    if postprocess_kwargs is None:
        postprocess_kwargs = {}
    psuffixes = []
    if which_x in ['Ra_i', 'Ra_i_eff', 'h_components', 'Ra_F_eff']:
        psuffixes.append('_T')
    elif which_x not in ['eta', 'Ra']:
        # probably in temp thign also
        print(' WARNING: possibly not implemented x variable')
        psuffixes.append('_T')
    if which_x in ['Ra_F_eff']:
        psuffixes.append('_Nu')

    df, df1 = None, None
    if not (not psuffixes):
        df = pro.pickleio_multi(case, psuffixes=psuffixes, t1=t1, load=load,
                                data_path=data_path, postprocess_kwargs=postprocess_kwargs, **kwargs)

        if averagescheme == 'timelast':
            df1 = df.mean(axis=0)
        elif averagescheme == 'timefirst':
            if '_T' in psuffixes:
                # load time-averages
                T_av, y = pro.time_averaged_profile_from_df(df, 'T_av')
                uv_mag_av, y = pro.time_averaged_profile_from_df(df, 'uv_mag_av')
                dic_av = pro.T_parameters_at_sol(case, n=None, T_av=T_av, uv_mag_av=uv_mag_av, y=y,
                                                 **postprocess_kwargs, **kwargs)  # actually a dict
                # really hacky bit
                for k in ['T_av', 'uv_mag_av', 'y']:
                    dic_av.pop(k, None)
                    df = df.drop(k, axis=1)  # drop lists you don't need
                df_av = pd.DataFrame({key: value for (key, value) in dic_av.items()}, index=[0])
                df1 = df.mean(axis=0).to_frame().transpose()  # mean of other parameters
                df1.set_index(pd.Series([0]))
                df1.update(df_av)  # update with properly timefirst-averaged temperature params

            else:
                print('not-implemented timefirst average with this x variable')
        else:
            # use each xy point (y=h) for fitting
            df1 = df
    x = getx_fromdf(Ra, eta, df=df1, case=case, which_x=which_x, averagescheme=averagescheme,
                    data_path=data_path_bullard,
                    t1=t1, load=load, postprocess_kwargs=postprocess_kwargs, **kwargs)

    if return_all:
        x_all = getx_fromdf(Ra, eta, df=df, case=case, which_x=which_x, averagescheme=None,
                            data_path=data_path_bullard,
                            t1=t1, load=load, postprocess_kwargs=postprocess_kwargs, **kwargs)
        return x, x_all
    else:
        return x


def getx_fromdf(Ra, eta, df=None, case=None, which_x=None, averagescheme=None, data_path=data_path_bullard,
                t1=None, load=None, postprocess_kwargs=None, **kwargs):
    if df is not None:
        try:
            keys = df.columns
        except AttributeError:
            keys = df.keys()
    if ('h_components' in which_x) or ((df is not None) and (which_x in keys)):
        try:
            missing = (which_x not in df.columns) or ((which_x in df.columns) and df[which_x].isnull().values.any())
        except AttributeError:
            missing = (which_x not in df.keys()) or ((which_x in df.keys()) and np.isnan(df[which_x]).any())
        if not missing:
            try:
                x = df[which_x].to_numpy()
            except AttributeError:  # if a float
                x = df[which_x]
        else:
            if averagescheme == 'timelast':
                print('    Averaging T components calcualted at each timestep')
            elif averagescheme == 'timefirst':
                print('    Calculating T components using time-averaged profiles')
            elif missing:
                print('    Calculating T components')
            else:
                print('    Unrecognized conditions: hoping for the best')
            x = pro.T_components_of_h(case, df=df, data_path=data_path, t1=t1, load=load, **postprocess_kwargs,
                                      **kwargs)

    elif 'Ra_i_eff' in which_x:  # calculate effective Ra using time-mean of T field params
        x = pro.Ra_i_eff(Ra_1=float(Ra), d_eta=float(eta), T_i=df['T_i'].to_numpy(),
                         T_l=df['T_l'].to_numpy(), delta_L=df['delta_L'].to_numpy())
    elif 'Ra_i' in which_x:
        x = pro.Ra_interior(Ra_1=float(Ra), d_eta=float(eta), T_i=df['T_i'].to_numpy())

    elif 'Ra_F_eff' in which_x:
        x = pro.Ra_F_eff(d_eta=float(eta), T_i=df['T_i'].to_numpy(), delta_L=df['delta_L'].to_numpy(),
                         q_sfc=df['Nu'].to_numpy())

    elif 'Ra' in which_x:
        x = float(Ra)

    elif 'eta' in which_x:
        x = float(eta)

    else:
        raise Exception('Invalid variable for x-axis / not implemented: ' + which_x)

    # make sure u can do numpy things
    convert = False
    try:
        np.log10(x)
    except TypeError:
        convert = True
    if (isinstance(x, pd.Series)):
        convert = True
    if convert:
        try:
            # print('x', x)
            x = x.to_numpy()
            # print('converted', x)
        except AttributeError:
            # dunno
            pass
    if iterable_not_string(x) and (averagescheme is not None):
        x = x[0]
    return x


def plot_geth(case=None, averagescheme=None, data_path=data_path_bullard, return_all=False,
              t1=None, postprocess_kwargs=None, load=True, **kwargs):
    # get the y values, depending on averaging scheme

    if averagescheme in ['timelast', 'timefirst']:
        df = pro.pickleio_multi(case, psuffixes=['_h_all'], t1=t1, load=load,
                                data_path=data_path, postprocess_kwargs=postprocess_kwargs, **kwargs)
        h_rms = df.h_rms.mean()
        h_peak = df.h_peak.mean()
    # elif averagescheme == 'timefirst':
    #     # load time-averages
    #     df_h = pickleio_average(case, suffix='_h_mean', postprocess_fn=h_timeaverage, t1=t1, load=load,
    #                             data_path=data_path, postprocess_kwargs=postprocess_kwargs, **kwargs)
    #     h_rms = df_h.iloc[0].h_rms
    #     h_peak = df_h.iloc[0].h_peak
    else:
        # use each xy point (y=h) for fitting
        df = pro.pickleio_multi(case, psuffixes=['_h'], t1=t1, load=load,
                                data_path=data_path, postprocess_kwargs=postprocess_kwargs, **kwargs)
        h_rms = df.h_rms.to_numpy()
        h_peak = df.h_peak.to_numpy()

    if return_all:
        return h_rms, h_peak, df.h_rms.to_numpy(), df.h_peak.to_numpy()
    else:
        return h_rms, h_peak


def plot_h_vs_2component(Ra=None, eta=None, t1_grid=None, end_grid=None, load_grid='auto', data_path=data_path_bullard,
                         fig_path=fig_path_bullard, averagescheme=None, p_dimensionals=None,
                         fig_fmt='.png', which_xs=None, include_regimes=None, regime_grid=None, legend=False,
                         save=True, fname='h', clabel=None, cbar=True, clabelpad=17, sigma=1, showpeak=False,
                         labelsize=16, xlabel='', ylabel='dynamic topography', y2label='', title='', cmap='winter',
                         fit=False, logx=True, logy=True, hscale=1, show_isoviscous=False, vmin=None, vmax=None,
                         fig=None, ax=None, ylim=None, xlim=None, postprocess_kwargs={}, regime_names=None, **kwargs):
    # for fitting h to 2 component power law, which_xs is iterable of x values, the first of which is the plotting one
    Ra, eta, (t1_grid, load_grid, end_grid, regime_grid) = pro.reshape_inputs(Ra, eta,
                                                                              (t1_grid, load_grid, end_grid,
                                                                               regime_grid))
    if include_regimes is None:
        include_regimes = regime_names
    if clabel is None:
        clabel = which_xs[1]
    if ax is None:
        fig = plt.figure()
        ax = plt.gca()

    quants = dict.fromkeys(['h_rms', 'h_peak', *which_xs])
    yx_peak_all, yx_rms_all = [], []
    D_m2_all = []

    # loop over cases
    for jj, etastr in enumerate(eta):
        cases, cases_var = pro.get_cases_list(Ra, etastr, end_grid[jj])
        for ii, case in enumerate(cases):
            if regime_grid[jj][ii] in include_regimes:
                t1_ii = t1_grid[jj][ii]
                load_ii = load_grid[jj][ii]
                # dat = ad.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False, read_statistics=True)

                # load x values for plotting
                xs = [plot_getx(Ra[ii], etastr, case=case, which_x=which_x, data_path=data_path,
                                averagescheme=averagescheme, return_all=True,
                                t1=t1_ii, load=load_ii, postprocess_kwargs=postprocess_kwargs, **kwargs) for which_x in
                      which_xs]
                x = [a[0] for a in xs]  # returned value according to average scheme
                x_all = [a[1] for a in xs]  # all data for getting errorbars

                # get the y values, depending on averaging scheme
                h_rms, h_peak, h_rms_all, h_peak_all = plot_geth(case=case, t1=t1_ii, data_path=data_path,
                                                                 averagescheme=averagescheme, return_all=True,
                                                                 postprocess_kwargs=postprocess_kwargs, **kwargs)
                try:
                    len(x_all[0])
                except TypeError:  # e.g. scalar
                    x_all[0] = [x_all[0]] * len(h_rms_all)
                try:
                    len(x_all[1])
                except TypeError:  # e.g. scalar
                    x_all[1] = [x_all[1]] * len(x_all[0])  # this way all three will always be same length
                #
                # # calculate Mahalanobis distance for chi square later
                # div = int(np.ceil(len(h_rms_all) / len(x_all[0])))
                # try:
                #     data = pd.DataFrame(
                #         {'y': np.log10(h_rms_all[::div]), 'x0': np.log10(x_all[0]), 'x1': np.log10(x_all[1])})
                # except (TypeError, AttributeError) as e:
                #     pee = np.asarray(h_rms_all).astype(np.float64)[::div]
                #     poo = np.asarray(x_all[0]).astype(np.float64)
                #     poo2 = np.asarray(x_all[1]).astype(np.float64)
                #     data = pd.DataFrame({'y': np.log10(pee), 'x0': np.log10(poo), 'x1': np.log10(poo2)})
                # d_m = mahalanobis(x=data, data=data, cov=None)
                # D_m2 = np.mean(d_m ** 2)
                # D_m2_all.append(D_m2)

                # append to working
                yx_peak_all.append((h_peak, x))
                yx_rms_all.append((h_rms, x))
                qdict = pro.parameter_percentiles(case, df={'h_rms': h_rms_all, 'h_peak': h_peak_all,
                                                            **{key: value for (key, value) in zip(which_xs, x_all)}},
                                                  keys=quants.keys(), sigma=sigma, plot=False)
                for key in quants.keys():
                    try:
                        quants[key] = np.vstack((quants[key], qdict[key]))
                    except ValueError:  # haven't added anything yet
                        quants[key] = qdict[key].reshape((1, 3))  # reshape so it works if you just have one row

    # get errorbars and plot them
    err = dict.fromkeys(quants.keys())
    z_vec = quants[which_xs[1]][:, 1]
    c_list = colorize(np.log10(z_vec), cmap=cmap, vmin=vmin, vmax=vmax)[0]
    try:
        for key in quants.keys():
            err[key] = np.array([quants[key][:, 1] - quants[key][:, 0], quants[key][:, 2] - quants[key][:, 1]])
        for jj, z in enumerate(z_vec):
            # get subset of points with this z-value
            ind = np.nonzero(z_vec == z)[0]
            if showpeak:
                ax.errorbar(quants[which_xs[0]][ind, 1], quants['h_peak'][ind, 1],
                            yerr=err['h_peak'].T[ind].T,
                            xerr=err[which_xs[0]].T[ind].T,
                            elinewidth=0.5, fmt='d', mfc=c_list[jj], c=c_list[jj], capsize=5, alpha=0.5,
                            markeredgecolor=highlight_colour)
            ax.errorbar(quants[which_xs[0]][ind, 1], quants['h_rms'][ind, 1],
                        yerr=err['h_rms'].T[ind].T,
                        xerr=err[which_xs[0]].T[ind].T,
                        elinewidth=0.5, fmt='o', mfc=c_list[jj], c=c_list[jj], capsize=5)

    except TypeError:  # no cases in given regimes
        pass

    if fit:
        ax = fit_cases_on_plot(yx_rms_all, ax, dist=D_m2_all,
                               c_list=colorize(np.log10(np.unique(z_vec)), cmap=cmap, vmin=vmin, vmax=vmax)[0],
                               labelsize=labelsize, n_fitted=len(which_xs) + 1, cmap=cmap, **kwargs)

    if show_isoviscous:
        df_JFR = pro.read_JFR('2Dcart_fixed_T_stats_updated.csv', path='/raid1/cmg76/aspect/benchmarks/JFR/')
        Ra_iso = df_JFR['Ra']
        h_rms_iso = df_JFR['RMS_topo']
        ax.plot(Ra_iso, h_rms_iso, c='k', ls='--', lw=0.5)

    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])  # for fair comparison
    if xlim is not None:
        ax.set_xlim(xlim)
    if legend and showpeak:
        leg = ax.legend(handles=[mlines.Line2D([], [], color='k', marker='d', alpha=0.5,
                                               markersize=15, markeredgecolor=highlight_colour,
                                               label=r'$h_{peak}$, data'),
                                 mlines.Line2D([], [], color='k', marker='o', alpha=0.5,
                                               markersize=15, label=r'$h_{rms}$, data')])
        ax.add_artist(leg)
    if cbar:
        dum = ax.scatter(z_vec, z_vec, c=z_vec, cmap=cmap, vmin=vmin, vmax=vmax, visible=False, zorder=0,
                         norm=LogNorm())
        cb = colourbar(dum, label=clabel, labelsize=labelsize, labelpad=clabelpad)
    ax.set_ylabel(ylabel, fontsize=labelsize)
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_title(title, fontsize=labelsize)

    if p_dimensionals is not None:
        ax2 = ax.twinx()
        # ax2.set_ylabel(y2label, fontsize=labelsize)
        ymin, ymax = ax.get_ylim()
        # apply function and set transformed values to right axis limits
        ax2.set_ylim((pro.dimensionalise_h(ymin, p_dimensionals), pro.dimensionalise_h(ymax, p_dimensionals)))
        # set an invisible artist to twin axes
        # to prevent falling back to initial values on rescale events
        ax2.plot([], [])
    if save:
        plot_save(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt)
    return fig, ax


def plot_h_vs(Ra=None, eta=None, t1_grid=None, end_grid=None, load_grid='auto', data_path=data_path_bullard,
              fig_path=fig_path_bullard, averagescheme=None, p_dimensionals=None, fiterror=False,
              fig_fmt='.png', which_x=None, include_regimes=None, regime_grid=None,
              save=True, fname='h', legend=False, sigma=1, showpeak=False, cleglabels=None,
              labelsize=16, legsize=16, xlabel='', ylabel='dynamic topography', y2label='', title='', alpha=1,
              c_peak='xkcd:forest green', c_rms='xkcd:periwinkle', cmap=None, c_fit='k', ms=10, lw=1, z_name='eta',
              xlabelpad=10, ylabelpad=10, elw=1, ecapsize=5, errortype='time', ticksize=None, vmin=None, vmax=None,
              fit=False, logx=True, logy=True, hscale=1, show_isoviscous=False, figsize=(7, 7), cbar=False,
              fig=None, ax=None, ylim=None, xlim=None, postprocess_kwargs=None, regime_names=None, **kwargs):
    if postprocess_kwargs is None:
        postprocess_kwargs = {}
    Ra, eta, (t1_grid, load_grid, end_grid, regime_grid) = pro.reshape_inputs(Ra, eta,
                                                                              (t1_grid, load_grid, end_grid,
                                                                               regime_grid))
    if include_regimes is None:
        include_regimes = regime_names
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = plt.gca()
    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])  # for fair comparison
    if xlim is not None:
        ax.set_xlim(xlim)
    if z_name == 'eta':
        z_vec = [float(eee) for eee in eta]  # the colourised vector, must be numeric
    else:
        raise Exception('not implemented colouring scheme, use z_name=eta')
    colourful = False
    if iterable_not_string(c_rms) or (cmap is not None):
        colourful = True
    quants = dict.fromkeys(['h_rms', 'h_peak', which_x])
    yx_peak_all, yx_rms_all = [], []
    D_m2_all = []
    sdy_all, sdx_all = [], []
    zz_all = []  # the colourised vector

    # loop over cases
    for jj, etastr in enumerate(eta):
        cases, cases_var = pro.get_cases_list(Ra, etastr, end_grid[jj])
        for ii, case in enumerate(cases):
            if regime_grid[jj][ii] in include_regimes:
                t1_ii = t1_grid[jj][ii]
                load_ii = load_grid[jj][ii]
                # dat = ad.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False, read_statistics=True)

                # extract x values for plotting
                x, x_times = plot_getx(Ra[ii], etastr, case=case, which_x=which_x, return_all=True, data_path=data_path,
                                       t1=t1_ii, load=load_ii, postprocess_kwargs=postprocess_kwargs,
                                       averagescheme=averagescheme, **kwargs)

                # get the y values, depending on averaging scheme
                h_rms, h_peak, h_rms_times, h_peak_times = plot_geth(case=case, t1=t1_ii, return_all=True,
                                                                     data_path=data_path, averagescheme=averagescheme,
                                                                     postprocess_kwargs=postprocess_kwargs, **kwargs)

                # # calculate Mahalanobis distance for chi square later
                # div = int(np.ceil(len(h_rms_times) / len(x_times)))
                # try:
                #     data = pd.DataFrame({'y': np.log10(h_rms_times[::div]), 'x': np.log10(x_times)})
                # except (TypeError, AttributeError) as e:
                #     pee = np.asarray(h_rms_times).astype(np.float64)[::div]
                #     poo = np.asarray(x_times).astype(np.float64)
                #     data = pd.DataFrame({'y': np.log10(pee), 'x': np.log10(poo)})
                # # V = np.cov(np.array([np.log10(h_rms_times[::div]), np.log10(x_times)]).T)
                # # try:
                # #     IV = np.linalg.inv(V)
                # # except np.linalg.LinAlgError:
                # #     IV = np.linalg.pinv(V)  # pseudo-inverse
                # d_m = mahalanobis(x=data, data=data, cov=None)
                # # d_m = distance.mahalanobis(np.log10(h_rms_times[::div]), np.log10(x_times), IV)
                # D_m2 = np.mean(d_m**2)
                # # D_m2 = np.var(np.log10(h_rms_times))
                # D_m2_all.append(D_m2)

                # calculate statistics
                sdy = np.nanstd(h_rms_times)
                sdx = np.nanstd(x_times)
                d_times = {'h_rms': h_rms_times, 'h_peak': h_peak_times, which_x: x_times}
                qdict = pro.parameter_percentiles(case,
                                                  df=d_times,
                                                  keys=quants.keys(), plot=False, sigma=2)

                # append to working
                yx_peak_all.append((h_peak, x))
                yx_rms_all.append((h_rms, x))
                sdy_all.append(sdy)
                sdx_all.append(sdx)
                if z_name == 'eta':
                    zz_all.append(jj)  # the colourised vector

                for key in quants.keys():
                    if errortype is 'time':
                        try:
                            quants[key] = np.vstack((quants[key], qdict[key]))  # add to array of errors
                        except ValueError:  # haven't added anything yet
                            quants[key] = qdict[key].reshape((1, 3))  # reshape so it works if you just have one row
                    elif errortype is 'standard':
                        try:
                            SE_mean = np.std(d_times[key]) / np.sqrt(len(d_times[key]))  # todo: log!
                        except TypeError:  # key is a scalar
                            SE_mean = 0
                        avg = np.mean(d_times[key])
                        SE_vec = np.array([avg - SE_mean, avg, avg + SE_mean])
                        # print(key, 'mean:', avg, 'SE of mean:', SE_mean)
                        d_times[key + '_SE'] = SE_mean
                        try:
                            quants[key] = np.vstack((quants[key], SE_vec))  # add to array of errors
                        except ValueError:  # haven't added anything yet
                            quants[key] = SE_vec.reshape((1, 3))  # reshape so it works if you just have one row

    # get errorbars and plot them
    err = dict.fromkeys(quants.keys())
    means = dict.fromkeys(quants.keys())
    means['h_peak'] = np.asarray([np.mean(a[0]) for a in yx_peak_all])
    means['h_rms'] = np.asarray([np.mean(a[0]) for a in yx_rms_all])
    means[which_x] = np.asarray([np.mean(a[1]) for a in yx_rms_all])
    # print('means x', means[which_x])
    try:
        for key in quants.keys():
            err[key] = np.asarray([quants[key][:, 1] - quants[key][:, 0], quants[key][:, 2] - quants[key][:, 1]])
        if showpeak:
            ax.errorbar(means[which_x], means['h_peak'], yerr=sdy_all,  # err['h_peak'],
                        xerr=sdx_all,  # err[which_x],
                        elinewidth=0.5, ms=ms,
                        fmt='d', c=c_peak, alpha=alpha, capsize=5, markeredgecolor=highlight_colour)
        mark = 'o'
        if colourful:
            if vmin is None:
                vmin = np.min(zz_all)
            if vmax is None:
                vmax = np.max(zz_all)
            if (cmap is not None) and (c_rms is None):
                try:
                    c_rms = colorize(zz_all, cmap=cmap, vmin=vmin, vmax=vmax)[0]
                except Exception as e:
                    cmap = cmap_from_ascii(cmap, path=cmap_path, end='.txt', ncol=4)
                    c_rms = colorize(zz_all, cmap=cmap, vmin=vmin, vmax=vmax)[0]
            for pp in range(len(means[which_x])):
                ax.errorbar(means[which_x][pp], means['h_rms'][pp], yerr=np.asarray([err['h_rms'][:, pp]]).T,
                            xerr=np.asarray([err[which_x][:, pp]]).T, elinewidth=elw, alpha=alpha,
                            fmt=mark, c=c_rms[zz_all[pp]], capsize=ecapsize, ms=ms)
        else:
            ax.errorbar(means[which_x], means['h_rms'], yerr=err['h_rms'], xerr=err[which_x], elinewidth=elw,
                        fmt=mark, c=c_rms, capsize=ecapsize, ms=ms)
    except TypeError as e:
        if quants['h_rms'] is None:  # no cases in given regimes as quants is dict of None
            pass
        else:
            raise e

    if cbar and colourful:
        if cmap is None:
            print('cbar not implemented without cmap')
        else:
            cax = colourbar(mappable=None, ax=ax, vmin=vmin, vmax=vmax, label='', labelsize=labelsize,
                            ticksize=ticksize, ticks=[float(z) for z in z_vec],
                            ticklabels=None, labelpad=17,
                            rot=None, discrete=False, cmap=cmap, tickformatter=None, pad=0.05, log=True)
    elif legend and colourful:
        # show colours outside
        clist = []
        for zz, colour in enumerate(c_rms):  # only used z vec
            if zz in np.unique(zz_all):
                clist.append(colour)
        ax = colourised_legend(ax, clist=clist, cleglabels=cleglabels, lw=0, ls='--', marker=mark, markersize=ms,
                               legsize=legsize, ncol=1)

    if fit:
        if which_x == 'h_components':
            n_fitted = 1
        else:
            n_fitted = 2
        if fiterror and errortype is 'time':
            xerr = sdx_all
            yerr = sdy_all
        elif fiterror and errortype is 'standard':
            xerr = d_times[which_x + '_SE']
            yerr = d_times['h_rms_SE']
        else:
            xerr = 1
            yerr = 1
        ax = fit_cases_on_plot(yx_rms_all, ax, c=c_fit, labelsize=labelsize, n_fitted=n_fitted, dist=D_m2_all,
                               xerr=xerr, yerr=yerr, legend=legend, lw=lw,
                               sigma=sigma, **kwargs)

    ax.set_ylabel(ylabel, fontsize=labelsize, labelpad=ylabelpad)
    ax.set_xlabel(xlabel, fontsize=labelsize, labelpad=xlabelpad)
    ax.set_title(title, fontsize=labelsize)

    if show_isoviscous:
        df_JFR = pro.read_JFR('2Dcart_fixed_T_stats_updated.csv', path='/raid1/cmg76/aspect/benchmarks/JFR/')
        Ra_iso = df_JFR['Ra']
        h_rms_iso = df_JFR['RMS_topo']
        ax.plot(Ra_iso, h_rms_iso, c='k', ls='--', lw=0.5)
    if legend and showpeak:
        leg = ax.legend(handles=[mlines.Line2D([], [], color=c_peak, marker='d', alpha=0.5,
                                               markersize=10, markeredgecolor=highlight_colour,
                                               label=r'$h_{peak}$, data'),
                                 mlines.Line2D([], [], color=c_rms, marker='o', alpha=0.5,
                                               markersize=10, label=r'$h_{rms}$, data')], fontsize=legsize)
        ax.add_artist(leg)
    if p_dimensionals is not None:
        def h_todim(u):
            return u * p_dimensionals['alpha_m'] * p_dimensionals['dT_m'] * p_dimensionals['d_m']

        def h_tonondim(u):
            return u / (p_dimensionals['alpha_m'] * p_dimensionals['dT_m'] * p_dimensionals['d_m'])

        ax2 = ax.secondary_yaxis('right', functions=(h_todim, h_tonondim))
        ax2.set_ylabel(y2label, fontsize=labelsize, labelpad=ylabelpad + 60, rotation=270)
        ax2.tick_params(axis='y', which='major', labelsize=ticksize)
        axes = (ax, ax2)
    else:
        axes = (ax,)

    if save:
        plot_save(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt)
    return (fig, *axes)


def subplots_topo_regimes(Ra_ls, eta_ls, regime_grid, regime_names, c_regimes=None, save=True, t1_grid=None, nrows=2,
                          ncols=2, which_x='Ra', leftleg_bbox=(-0.05, 1), p_dimensionals=None,
                          load_grid='auto', fig_path=fig_path_bullard, fname='h_Ra_all', fig_fmt='.png', end_grid=None,
                          show_bounds=False, regimes_title='', Ra_i=False, show_isoviscous=False, y2label='',
                          labelsize=14, xlabel='Ra', include_regimes=None, ylabel='dynamic topography', xlabelpad=12,
                          ylabelpad=2, **kwargs):
    Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid) = pro.reshape_inputs(Ra_ls, eta_ls, (
        t1_grid, load_grid, end_grid, regime_grid))
    if include_regimes is None:
        include_regimes = regime_names
    if c_regimes is None:
        c_regimes = ['xkcd:sage green', 'xkcd:blood red', 'xkcd:dark violet']
    if show_isoviscous:
        show_isoviscous_flag = True
    else:
        show_isoviscous_flag = False

    print(r'Plotting h vs.', which_x, 'for', include_regimes, 'regimes')
    fig = plt.figure(figsize=(7, 7))
    bigax = fig.add_subplot(111)  # The big subplot
    bigax.spines['top'].set_color('none')
    bigax.spines['bottom'].set_color('none')
    bigax.spines['left'].set_color('none')
    bigax.spines['right'].set_color('none')
    bigax.tick_params(labelcolor='w', which='both', bottom=False, left=False, right=False, top=False)
    bigax.set_xlabel(xlabel, fontsize=labelsize, labelpad=xlabelpad)
    bigax.set_ylabel(ylabel, fontsize=labelsize, labelpad=ylabelpad)
    if p_dimensionals is not None:
        bigax2 = bigax.twinx()
        bigax2.spines['top'].set_color('none')
        bigax2.spines['bottom'].set_color('none')
        bigax2.spines['left'].set_color('none')
        bigax2.spines['right'].set_color('none')
        bigax2.tick_params(labelcolor='w', which='both', bottom=False, left=False, right=False, top=False)
        bigax2.set_ylabel(y2label, fontsize=labelsize)

    for ii, eta_ii in enumerate(eta_ls):
        print(r' $\Delta \eta$:', eta_ii, ' (', ii, '/', len(eta_ls) - 1, ')')
        z = int(str(nrows) + str(ncols) + str(ii + 1))
        ax = fig.add_subplot(z)
        t1_ii = t1_grid[ii]
        end_ii = end_grid[ii]
        load_ii = load_grid[ii]

        for ir, regime_name in enumerate(include_regimes):
            Ra_regime = [Ra_ls[j] for j in np.nonzero(regime_grid[ii] == regime_name)[0]]
            Ra_regime_idx = [j for j in np.nonzero(regime_grid[ii] == regime_name)[0]]

            if ir == 0 and show_isoviscous_flag:
                show_isoviscous = True
            else:
                show_isoviscous = False

            if not (not Ra_regime):  # if this regime is not empty
                fig, ax = plot_h_vs(Ra_regime, eta_ii, t1_grid=t1_ii[Ra_regime_idx], end_grid=end_ii[Ra_regime_idx],
                                    load_grid=load_ii[Ra_regime_idx], regime_grid=regime_grid[ii, Ra_regime_idx],
                                    p_dimensionals=p_dimensionals, which_x=which_x,
                                    save=False, labelsize=labelsize, xlabel='', ylabel='', y2label=y2label,
                                    c_peak=c_regimes[ir], c_rms=c_regimes[ir], Ra_i=Ra_i, regime_names=regime_names,
                                    show_isoviscous=show_isoviscous, fig=fig, ax=ax, **kwargs)
                if show_bounds:
                    ax.axvline(float(Ra_regime[-1]) * 2, c='k', lw=0.5, alpha=0.6, ls='--')
                    # ax.text(ax.get_xlim()[0], ylim[0], regime_name, fontsize=8, va='bottom', ha='left')
                # print('Plotted', len(Ra_regime), regime_name, 'case(s)')

        # ax.set_title(r'$\Delta \eta$=' + eta_ii, fontsize=labelsize-2)
        ax.text(0.01, 0.98, r'$\Delta \eta$=' + eta_ii, fontsize=labelsize - 4, ha='left', va='top',
                transform=ax.transAxes)  # label

        if ii % ncols != 0 and p_dimensionals is None:
            ax.yaxis.tick_right()
        elif ii % ncols != 0:
            ax.set_yticklabels([])

    # add legends
    ax = bigax
    handles1 = [ax.scatter([], [], label='peak', marker='d', edgecolors=highlight_colour, c='k'),
                ax.scatter([], [], label='rms', marker='o', c='k')]
    if show_isoviscous_flag:
        hplot, = ax.plot([], [], label='rms, 2D cartesian isoviscous', ls='--', c='k')
        handles1.append(hplot)
    outer_legend = ax.legend(handles=handles1,
                             borderaxespad=0., ncol=len(handles1), bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                             frameon=False,  # mode="expand"
                             )
    ax.add_artist(outer_legend)
    handles2 = []
    for ir, regime_name in enumerate(regime_names):
        handles2.append(ax.scatter([], [], label=regime_name, marker='o', c=c_regimes[ir]))
    regime_legend = ax.legend(handles=handles2,
                              borderaxespad=0., title=regimes_title, bbox_to_anchor=leftleg_bbox, loc='upper right',
                              frameon=False)
    ax.add_artist(regime_legend)

    fig.subplots_adjust(wspace=0.05, hspace=0.15, left=0.22)
    if save:
        plot_save(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt, bbox_inches=None,
                  bbox_extra_artists=(outer_legend, regime_legend), tight_layout=False)
    return fig, fig.axes


def plot_Ra_scaling(Ra_data=None, y_data=None, fig_path=fig_path_bullard,
                    save=True, fname='claire',
                    labelsize=16, ylabel='', xlabel='Ra', title='', fig_fmt='.png',
                    c_scatter='xkcd:forest green',
                    fit=False, logx=True, logy=True,
                    fig=None, ax=None, ylim=None, xlim=None, **kwargs):
    if fig is None:
        fig = plt.figure()
        ax = plt.gca()

    ax.plot(Ra_data, y_data, '-o', c=c_scatter)

    if fit:
        yx_all = [(y, r) for r, y in zip(Ra_data, y_data)]
        ax = fit_cases_on_plot(yx_all, ax, weights=None, c=c_scatter, labelsize=labelsize, **kwargs)

    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])  # for fair comparison
    if xlim is not None:
        ax.set_xlim(xlim)
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_title(title, fontsize=labelsize)
    ax.set_ylabel(ylabel, fontsize=labelsize)

    if save:
        plot_save(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt)
    return fig, ax


def subplots_Ra_scaling(Ra_ls=None, eta_ls=None, t1_grid=None, end_grid='', keys=None, data_path=data_path_bullard,
                        fig_path=fig_path_bullard, load_grid='auto', regime_grid=None, include_regimes=None, Ra_i=False,
                        compare_exponent=None,
                        save=True, fname='Ra_scalings', labelsize=16, ylabels=None, psuffixes='', title='',
                        xlim=None, ylim=None, legloc=None, averagescheme=None,
                        regime_names=None,
                        cmap='magma', compare_pub=None, compare_label=None, vmin=None, vmax=None,
                        fig=None, axes=None, fig_fmt='.png', postprocess_kwargs=None, **kwargs):
    # Ra or eta is list of strings, t1 is a list of numbers the same length
    # instead of plotting vs Ra or eta, plot vs theoretical components of scaling relationship

    if postprocess_kwargs is None:
        postprocess_kwargs = {}
    Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid) = pro.reshape_inputs(Ra_ls, eta_ls, (
        t1_grid, load_grid, end_grid, regime_grid))
    if include_regimes is None:
        include_regimes = regime_names
    if ylabels is None:
        ylabels = keys
    if iterable_not_string(keys):
        nkeys = len(keys)
    elif keys is None:
        raise Exception('No y-axis keys provided!')
    else:
        nkeys = 1
    if ylim is None:
        ylim = [None] * nkeys
    if legloc is None:
        legloc = ['lower left'] * nkeys
    if fig is None:
        fig, axes = plt.subplots(nkeys, 1, figsize=(7, nkeys * 2.5), sharex=True)
        if nkeys == 1:
            axes = np.array([axes])
    logeta_fl = [np.log10(float(a)) for a in eta_ls]
    c_list = colorize(logeta_fl, cmap=cmap, vmin=vmin, vmax=vmax)[0]

    for jj, eta_str in enumerate(eta_ls):
        cases, Ra_var = pro.get_cases_list(Ra_ls, eta_str, end_grid[jj])
        c_scatter = c_list[jj]

        plot_data = {'Ra': []}
        for key in keys:
            plot_data[key] = []

        for ii, case in enumerate(cases):
            t1_ii = t1_grid[jj][ii]
            load_ii = load_grid[jj][ii]

            if (t1_ii != 1) and (os.path.exists(data_path + 'output-' + case)):
                Ra_ii = float(Ra_var[ii])
                # load data
                if load_ii == 'auto':
                    dat = ad.Aspect_Data(directory=data_path + 'output-' + case + '/',
                                         read_statistics=False, read_parameters=False, **kwargs)
                else:
                    dat = None
                dfs = []
                for ip, suffix in enumerate(psuffixes):
                    df1 = pro.pickleio(case, suffix=suffix, t1=t1_ii,
                                       dat_new=dat, data_path=data_path, load=load_ii,
                                       postprocess_kwargs=postprocess_kwargs, **kwargs)
                    dfs.append(df1)
                try:
                    df = pd.concat(dfs, axis=1)
                    df = df.loc[:, ~df.columns.duplicated()]
                except Exception as e:
                    for dfi in dfs:
                        print(dfi)
                    raise e

                if averagescheme == 'timefirst':
                    print('    subplots_Ra_scaling(): Calculating T components using time-averaged profiles')
                    T_av, y = pro.time_averaged_profile_from_df(df, 'T_av')
                    uv_mag_av, y = pro.time_averaged_profile_from_df(df, 'uv_mag_av')
                    df_av = pro.T_parameters_at_sol(case, n=None, T_av=T_av, uv_mag_av=uv_mag_av, **postprocess_kwargs,
                                                    **kwargs)
                    if 'Nu' in keys:
                        df_av['Nu'] = np.mean(df['Nu'])  # using average flux anyways so scheme doesn't matter
                    df = df_av

                for key in keys:
                    med = np.median(df[key])
                    if np.isnan(med):
                        raise Exception('NaN in median, key:', key, '\n', df[key])
                    plot_data[key].append(np.median(df[key]))
                if Ra_i == 'eff':
                    plot_data['Ra'].append(pro.Ra_i_eff(Ra_1=Ra_ii, d_eta=float(eta_str), T_i=np.median(df['T_i']),
                                                        T_l=np.median(df['T_l']), delta_L=np.median(df['delta_L'])))
                elif Ra_i:
                    plot_data['Ra'].append(pro.Ra_interior(Ra_1=Ra_ii, d_eta=float(eta_str), T_i=np.median(df['T_i'])))
                else:
                    plot_data['Ra'].append(Ra_ii)

                if compare_pub is not None:
                    d_compare = compare_pub(Ra=Ra_ii, d_eta=float(eta_str), case=case, dat=dat, df=df,
                                            load=load_ii, **kwargs)
                    for k, key in enumerate(keys):
                        try:
                            axes[k].plot(d_compare['Ra_i'], d_compare[key], '^', c=c_scatter,
                                         markeredgecolor=highlight_colour, alpha=0.8, zorder=100)
                        except KeyError:
                            print('Key', key, 'not returned by', compare_pub)
                        except Exception as e:
                            print('d_compare[Ra_i]', d_compare['Ra_i'])
                            raise e
                if compare_exponent is not None:
                    for k, key in enumerate(keys):
                        axes[k].plot(plot_data['Ra'], np.array(plot_data['Ra']) ** compare_exponent[k],
                                     label=key + '^' + str(compare_exponent[k]), c='k', lw=1)

        for k, key in enumerate(keys):
            xlabel = ''
            if k == len(keys) - 1:
                if Ra_i == 'eff':
                    xlabel = r'Ra$_{i,eff}$'
                elif Ra_i:
                    xlabel = r'Ra$_i$'
                else:
                    xlabel = r'Ra$_1$'
            fig, axes[k] = plot_Ra_scaling(Ra_data=plot_data['Ra'], y_data=plot_data[key], xlim=xlim, ylim=ylim[k],
                                           save=False, labelsize=labelsize, ylabel=ylabels[k], c_scatter=c_scatter,
                                           fig=fig, ax=axes[k], xlabel=xlabel, legend=True, legloc=legloc[k], **kwargs)

    # if compare_pub is not None:  # add top legend
    ax = axes[0]
    if compare_pub is None:
        outer_handles = [None]
        outer_labels = [None]
        ncol = 1
    else:
        outer_handles = [ax.scatter([], [], marker='^', c=c_scatter, edgecolors=highlight_colour),
                         ax.scatter([], [], marker='o', c=c_scatter)]
        outer_labels = [compare_label, 'This work']
        ncol = len(outer_handles)
    outer_legend = ax.legend(handles=outer_handles, labels=outer_labels,
                             borderaxespad=0., ncol=ncol, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                             frameon=False,  # mode="expand"
                             )
    ax.add_artist(outer_legend)
    bbox_extra_artists = (outer_legend,)
    tight_layout = False
    bbox_inches = None
    fig.subplots_adjust(left=0.2)
    # else:
    #     bbox_extra_artists = None
    #     tight_layout = True
    #     bbox_inches = 'tight'

    # colorbar proxy artist
    scat = axes[-1].scatter(logeta_fl, logeta_fl, visible=False, c=np.array(logeta_fl), cmap=cmap,
                            vmin=vmin, vmax=vmax)  # dummy
    cbar = fig.colorbar(scat, ax=axes.ravel().tolist())
    print('axes:', axes.ravel().tolist())
    cbar.set_label(r'log($\Delta \eta$)', fontsize=labelsize, rotation=270, labelpad=22)
    plt.suptitle(title, fontsize=labelsize, y=1.02)

    if save:
        plot_save(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt, bbox_inches=bbox_inches,
                  bbox_extra_artists=bbox_extra_artists, tight_layout=tight_layout)
    return fig, axes


def subplots_cases(cases, labels=None, labelsize=16, labelpad=5, t1=None, save=True, dt_xlim=(0.0, 0.065),
                   fname='cases', data_path=data_path_bullard, fig_path=fig_path_bullard, fig_fmt='.png',
                   load='auto', includegraphic=False, c_rms='xkcd:forest green', c_peak='xkcd:periwinkle',
                   suptitle='', includepdf=True, includeTz=True, regime_grid=None, show_sols=False, **kwargs):
    # rows are cases, columns are v_rms, q, T(z), hist
    ncases = len(cases)
    ncols = 2
    if regime_grid is None:
        regime_grid = [None] * len(cases)
    if includepdf:
        ncols = ncols + 1
    if includeTz:
        ncols = ncols + 1
    if includegraphic:
        ncols = ncols + 1
    fig, axes = plt.subplots(ncases, ncols, figsize=(17, ncases * 2.5))
    if t1 is None:
        t1 = [0] * ncases
    numplotted = 0
    for ii, case in enumerate(cases):
        icol = 0
        sol_df = None  # make sure these are blanked
        ts_df = None
        if os.path.exists(data_path + 'output-' + case) and not (regime_grid[ii] == 'no convection') and not (
                regime_grid[ii] == 'sluggish'):
            print('Plotting summary for', case, 'using t1 =', t1[ii])
            t1_ii = t1[ii]
            if iterable_not_string(load):
                load_ii = load[ii]
            else:
                load_ii = load
            dat = ad.Aspect_Data(directory=data_path + 'output-' + case + '/', read_statistics=True, **kwargs)
            if ii == ncases - 1:  # show x label in bottom row only
                setxlabel = True
            else:
                setxlabel = False
            setylabel = True
            legend = False
            if numplotted == 0:
                legend = True

            if show_sols:  # load df
                if t1[ii] < 1:
                    sol_df = pro.pickleio(case, suffix='_T', t1=t1_ii,
                                          dat_new=dat, load=load_ii, data_path=data_path, fig_path=fig_path, **kwargs)
                else:
                    show_sols = False  # automatically override processing

            # rms velocity plot
            ax = axes[ii, icol]
            fig, ax = plot_evol(case, 'rms_velocity', fig=fig, ax=ax, save=False, mark_used=True, t1=t1_ii, dat=dat,
                                show_sols=show_sols, ylabel='rms velocity', c='k', settitle=False, setxlabel=setxlabel,
                                setylabel=setylabel, legend=False, labelsize=labelsize, labelpad=labelpad,
                                sol_df=sol_df)

            ax.text(0.01, 0.95, labels[ii], horizontalalignment='left', verticalalignment='top',
                    transform=ax.transAxes, fontsize=labelsize, zorder=4)

            # heat flux plot
            icol = icol + 1
            ax = axes[ii, icol]
            fig, ax = plot_evol(case, 'heatflux_top', fig=fig, ax=ax, save=False, mark_used=False, dat=dat,
                                show_sols=False, ylabel='heat flux', c='xkcd:light red', settitle=False,
                                setxlabel=setxlabel, setylabel=setylabel, labelsize=labelsize, labelpad=labelpad,
                                label='top')
            fig, ax = plot_evol(case, 'heatflux_bottom', fig=fig, ax=ax, save=False, mark_used=True, t1=t1_ii, dat=dat,
                                show_sols=show_sols, ylabel='heat flux', yscale=-1, c='xkcd:purple blue',
                                settitle=False, setxlabel=setxlabel, setylabel=setylabel, legend=legend,
                                labelsize=labelsize, labelpad=labelpad, label='bottom', sol_df=sol_df)

            if includeTz:
                icol = icol + 1
                ax = axes[ii, icol]
                if t1_ii < 1:
                    if not show_sols:
                        sol_df = pro.pickleio(case, suffix='_T', t1=t1_ii,
                                              dat_new=dat, load=load_ii, data_path=data_path, fig_path=fig_path,
                                              **kwargs)

                    fig, ax = plot_T_profile(case, T_params=sol_df, n='mean', data_path=data_path, setylabel=False,
                                             setxlabel=setxlabel, save=False, fig_path=fig_path, fig=fig, ax=ax,
                                             legend=legend, labelsize=labelsize)
                else:
                    ax.text(0.01, 0.95, '\n\n\nnot converged?',
                            horizontalalignment='left', verticalalignment='top',
                            transform=ax.transAxes, fontsize=labelsize)
                if setxlabel:
                    ax.set_xlabel('mean temperature', fontsize=labelsize)
                if setylabel:
                    ax.set_ylabel('depth', fontsize=labelsize)

            if includepdf:
                icol = icol + 1
                ax = axes[ii, icol]
                if t1_ii < 1:
                    ts_df = pro.pickleio(case, suffix='_h_all', t1=t1_ii,
                                         dat_new=dat, load=load_ii, data_path=data_path, fig_path=fig_path,
                                         **kwargs)

                    fig, ax = plot_pdf(case, df=ts_df, keys=['h_rms', 'h_peak'], fig=fig, ax=ax, save=False,
                                       settitle=False, setxlabel=setxlabel, legend=legend, labelsize=labelsize,
                                       c_list=[c_rms, c_peak], path=data_path)
                else:
                    ax.text(0.01, 0.95, '\n\n\nnot converged?',
                            horizontalalignment='left', verticalalignment='top',
                            transform=ax.transAxes, fontsize=labelsize)
                ax.set_xlim(dt_xlim[0], dt_xlim[1])  # for fair comparison

            if includegraphic:
                icol = icol + 1
                ax = axes[ii, icol]
                fgraph = fig_path + 'graphical/' + case + '.png'
                try:
                    img = mpimg.imread(fgraph)
                    ax.imshow(img)
                    print('    Grabbing graphical output for', case)
                except FileNotFoundError:
                    try:
                        fig, ax = static_T_field(case, data_path=data_path, labelsize=labelsize, ticksize=10,
                                                 cmap='gist_heat', c='k', cbar=False, title='', fig=fig, ax=ax,
                                                 shading='nearest', return_artists=False, save=False, i_n=-1, avg=False)
                        print('    Plotting graphical output from solution for', case)
                    except Exception as e:
                        print(e)
                        print('    Graphical output not found:', fgraph)
                        ax.text(0.01, 0.95, '\n\n\nno image saved',
                                horizontalalignment='left', verticalalignment='top',
                                transform=ax.transAxes, fontsize=labelsize)
                        # fig.delaxes(ax)

            numplotted += 1
        else:
            print('Case', case, 'not found')
            ax = axes[ii, 0]
            ax.text(0.01, 0.95, labels[ii] + '\n\n' + str(regime_grid[ii]),
                    horizontalalignment='left', verticalalignment='top',
                    transform=ax.transAxes, fontsize=labelsize)
    plt.suptitle(suptitle, fontsize=labelsize * 2, y=1.02)
    if save:
        fig.tight_layout()
        plot_save(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt)
    return fig, axes


def get_h_average(Ra, eta, which_h='rms', end=None, data_path=data_path_bullard, load=True, **kwargs):
    case = 'Ra' + Ra + '-eta' + eta + end
    rms, peak = plot_geth(case=case, averagescheme='timefirst', data_path=data_path, load=load, **kwargs)
    if which_h == 'rms':
        return rms
    elif which_h == 'peak':
        return peak
    else:
        raise Exception('invalid which h')


def plot_fit_parameter_grid(Ra_ls, eta_ls, data_path=data_path_bullard, fig_path=fig_path_bullard, load_grid=None,
                            vmin=None, vmax=None, averagescheme=None, which_x=None, regime_grid=None,
                            include_regimes=['chaotic'],
                            save=True, fname='fit-grid', labelsize=16, fig_fmt='.png', t1_grid=None, end_grid=None,
                            cticklabels=None,
                            cticks=None, title='', lognorm=False, log=False, clabel=r'$h_{rms}$', which_h='rms',
                            nlevels_contour=10, cmap='Greys_r', cmap_contours='spring', postprocess_kwargs={},
                            **kwargs):
    # make grid
    fig, ax = plot_parameter_grid(Ra_ls, eta_ls, function=get_h_average, data_path=data_path, load=True, cmap=cmap,
                                  vmin=vmin, vmax=vmax, save=False, labelsize=16, t1_grid=t1_grid, end=end_grid,
                                  cticklabels=cticklabels,
                                  cticks=cticks, title=title, lognorm=False, log=False, clabel=clabel, overplot_h=False,
                                  which_h=which_h, **kwargs)

    # # get fit things
    const, expon = plot_model_data(Ra_ls, eta_ls, regime_grid=regime_grid, t1_grid=t1_grid, load_grid=load_grid,
                                   end_grid=end_grid, literature_file=None, legend=False, averagescheme=averagescheme,
                                   postprocess_kwargs=postprocess_kwargs,
                                   which_h=which_h, which_x=which_x, data_path=data_path,
                                   save=False, cbar=None, include_regimes=include_regimes,
                                   intercept=False)

    # get fitted h values for contours
    Ra = np.array([float(r) for r in Ra_ls])
    eta = np.array([float(e) for e in eta_ls])
    Rv, ev = np.meshgrid(Ra, eta)
    idx = regime_grid == include_regimes
    Ra_r = np.logspace(np.log10(Rv[idx].min()), np.log10(Rv[idx].max()), endpoint=True)  # in the regime
    eta_r = np.logspace(np.log10(ev[idx].min()), np.log10(ev[idx].max()), endpoint=True)
    X, Y = np.meshgrid(Ra_r, eta_r)
    H = const * X ** expon[0] * Y ** expon[1]

    # normalize to axis limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    logRa2 = minmaxnorm(np.log10(Ra), xlim[0], xlim[1])
    logeta2 = minmaxnorm(np.log10(eta), ylim[0], ylim[1])
    logRv2, logev2 = np.meshgrid(logRa2, logeta2)
    logRa2_r = np.linspace(logRv2[idx].min(), logRv2[idx].max())  # in the regime
    logeta2_r = np.linspace(logev2[idx].min(), logev2[idx].max())
    X2, Y2 = np.meshgrid(logRa2_r, logeta2_r)

    CS = ax.contour(X2, Y2, H, nlevels_contour, cmap=cmap_contours)
    ax.clabel(CS, inline=1, fontsize=8)

    if save:
        plot_save(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt, tight_layout=False)


def plot_parameter_grid(Ra, eta, function, data_path=data_path_bullard, fig_path=fig_path_bullard, load='auto',
                        vmin=None, vmax=None, set_under=None, set_over=None,
                        save=True, fname='grid', labelsize=16, fig_fmt='.png', t1=None, end=None, cticklabels=None,
                        cticks=None, title='', lognorm=False, log=False, clabel=None, discrete=False,
                        cmap='jet', clist=None, **kwargs):
    # plot output of any (scalar-returning) function (with optional topo contours?)

    if t1 is None:
        t1 = np.zeros((len(eta), len(Ra)))
    if not_iterable(load):  #
        load = np.array([[load] * len(Ra)] * len(eta))
    fig, ax = plt.subplots(1, 1)
    plot_grid = np.zeros((len(eta), len(Ra)))

    for jj, eta_str in enumerate(eta):
        cases, _ = pro.get_cases_list(Ra, eta_str, end[jj])
        for ii, Ra_str in enumerate(Ra):
            # calculate value at this parameter-space coordinate
            if os.path.exists(data_path + 'output-' + cases[ii] + '/'):
                plot_grid[jj, ii] = function(Ra=Ra_str, eta=eta_str, ii=ii, jj=jj, load=load[jj][ii],
                                             t1=t1[jj][ii], end=end[jj][ii], data_path=data_path, **kwargs)
            else:
                plot_grid[jj, ii] = np.nan
    print('plot_grid', np.shape(plot_grid))
    if log:
        plot_grid = np.log10(plot_grid)
    m = np.ma.masked_where(np.isnan(plot_grid), plot_grid)
    print('m', m)
    if discrete and vmax is not None and vmin is not None and (vmax - vmin + 1) != len(cticklabels):
        print('Did you mean vmax - vmin =', len(cticklabels))
    if vmax is None:
        vmax = np.max(m)
    elif log:
        vmax = np.log10(vmax)
    if vmin is None:
        vmin = np.min(m)
    elif log:
        vmin = np.log10(vmin)

    if clist is None:
        if discrete:
            cmap = plt.cm.get_cmap(cmap, vmax - vmin)
        else:
            cmap = plt.cm.get_cmap(cmap)
    else:
        cmap = cmap_from_list(clist, cmap_name='regimes')

    cmap_extend = 'neither'
    if set_under is not None:
        cmap.set_under(set_under)
        cmap_extend = 'min'
    if set_over is not None:  # doesn't work for custom LinearSegmented colormap (bug: https://github.com/matplotlib/matplotlib/issues/4117/)
        cmap.set_over(set_over)
        if set_under is None:
            cmap_extend = 'max'
        else:
            cmap_extend = 'both'

    if lognorm:
        im_norm = LogNorm(vmin=vmin, vmax=vmax, clip=False)
    else:
        im_norm = None  # Normalize(vmin=vmin, vmax=vmax, clip=False)

    im = ax.imshow(m, origin='lower', aspect='equal', interpolation='None', cmap=cmap, vmin=vmin, vmax=vmax,
                   norm=im_norm)

    # draw grid lines
    for x, xval in enumerate(Ra):
        ax.axvline(x=x - 0.5, c='k', lw=0.5)
    for y, yval in enumerate(eta):
        ax.axhline(y=y - 0.5, c='k', lw=0.5)

    ax.set_xlabel('Ra', fontsize=labelsize)
    ax.set_ylabel(r'$\Delta \eta$', fontsize=labelsize)
    ax.set_xticks(np.arange(len(Ra)))
    ax.set_yticks(np.arange(len(eta)))
    ax.set_xticklabels(Ra)
    ax.set_yticklabels(eta)

    cbar = plt.colorbar(im, shrink=0.5, extend=cmap_extend)

    if cticks is not None:
        cbar.set_ticks(cticks)
    elif cticks is None and discrete:
        nlabels = len(cticklabels)
        tick_locs = (np.arange(vmin, vmax + 1) + 0.5) * (nlabels - 1) / nlabels
        cbar.set_ticks(tick_locs)
        # cbar.set_ticks(np.arange(vmin, vmax + 1))
    if cticklabels is not None:
        cbar.ax.set_yticklabels(cticklabels)
    if clabel is not None:
        cbar.set_label(clabel, rotation=270, labelpad=17, fontsize=labelsize)

    ax.set_title(title, fontsize=labelsize)
    if save:
        plot_save(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt, tight_layout=False)
    return fig, ax


def plot_velocity_profile(case, dat=None, n=None, xlabel='rms velocity', ylabel='depth', fig=None, ax=None,
                          data_path=data_path_bullard,
                          labelsize=16, fig_path=fig_path_bullard, fname='velocity', save=True, fig_fmt='.png',
                          **kwargs):
    if fig is None:
        fig, ax = plt.subplots(figsize=(4, 4))

    if dat is None:
        dat = ad.Aspect_Data(directory=data_path + 'output-' + case + '/',
                             read_statistics=False, read_parameters=False, **kwargs)
    if n is None:
        n = dat.final_step()

    x, y, _, u, v, _, mag = dat.read_velocity(n, **kwargs)
    mag_av = ad.horizontal_mean(mag, x)
    fig, ax = dat.plot_profile(mag_av, n=n, y=y, label='', ylabel=None, fig=fig, ax=ax, **kwargs)
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_ylabel(ylabel, fontsize=labelsize)
    ax.set_ylim(y.min(), y.max())

    if save:
        plot_save(fig, fig_path=fig_path, fname=fname, fig_fmt=fig_fmt)
    return fig, ax


def plot_T_profile(case, T_params=None, n=-1, dat=None, data_path=data_path_bullard, t1=0,
                   setylabel=True, setxlabel=True, save=True, load='auto',
                   fig_path=fig_path_bullard, fig=None, ax=None, fend='_T-z', fig_fmt='.png',
                   legend=True, labelsize=16, **kwargs):
    if T_params is None:
        T_params = pro.pickleio(case, suffix='_T', t1=t1,
                                dat_new=dat, load=load, data_path=data_path, fig_path=fig_path, **kwargs)
    if fig is None:
        fig, ax = plt.subplots(figsize=(4, 4))

    if n == 'mean':  # avg of all steady state sols
        T_params = T_params.mean(axis=0)  # T params df already only contains steady state values
    else:
        try:
            T_params = T_params.loc[T_params['sol'] == n]
        except IndexError:
            print('No T parameterisation found for solution n =', n)
            return fig, ax

    delta_rh_n = np.array(T_params['delta_rh'])  # ensure not list
    delta_0_n = np.array(T_params['delta_0'])
    D_l_n = np.array(T_params['y_L'])
    T_l_n = np.array(T_params['T_l'])
    T_i_n = np.array(T_params['T_i'])
    T_f = np.array(T_params['T_av'].tolist())
    y_f = np.array(T_params['y'].tolist())

    ax.plot(T_f, y_f, c='k', lw=1)
    ax.axhline(D_l_n, label='$\delta_{L}$', c='xkcd:tangerine', lw=0.5)
    ax.axhline(1 - delta_0_n, label=r'$\delta_0$', c='xkcd:red orange', lw=0.5)
    try:
        ax.text(0, 1 - delta_0_n, r'$\delta_{rh} = $' + '{:04.3f}'.format(delta_rh_n), ha='left', va='top',
                color='xkcd:red orange', fontsize=labelsize - 2)
    except TypeError:
        ax.text(0, 1 - delta_0_n, r'$\delta_{rh} = $' + '{:04.3f}'.format(delta_rh_n.item()), ha='left', va='top',
                color='xkcd:red orange', fontsize=labelsize - 2)
    ax.plot([T_l_n, T_l_n], [0, D_l_n], ls='--', alpha=0.5, lw=0.5, label=r'$T_L$', c='xkcd:tangerine')
    ax.plot([T_i_n, T_i_n], [0, 1 - delta_0_n], ls='--', alpha=0.5, lw=0.5,
            label=r'$T_i$', c='xkcd:red orange')
    if legend:
        ax.legend(frameon=True, fontsize=labelsize - 4, ncol=2)
    if setxlabel:
        ax.set_xlabel('temperature', fontsize=labelsize)
    if setylabel:
        ax.set_ylabel('depth', fontsize=labelsize)
    if save:
        plot_save(fig, case + fend, fig_path=fig_path, fig_fmt=fig_fmt)
    return fig, ax


def plot_pdf(case, df=None, keys=None, fig_path=fig_path_bullard, fig=None, ax=None, save=True, settitle=True,
             setxlabel=True, legend=True, labelsize=16, c_list=None, labels=None, fend='h_hist', fig_fmt='.png',
             xlabel='dynamic topography', **kwargs):
    if c_list is None:
        c_list = ['xkcd:forest green', 'xkcd:periwinkle']
    if labels is None:
        labels = ['rms', 'peak']
    if keys is None:
        keys = ['h_rms', 'h_peak']
    if ax is None:
        fig = plt.figure()
        ax = plt.gca()

    for ii, key in enumerate(keys):
        try:
            x = df[key]
            x = x[~np.isnan(x)]  # get rid of nan
            c = c_list[ii]
            try:
                ax.hist(x, color=c, histtype='step', label=labels[ii])
                extralabel = ''
                if ii == 0:
                    extralabel = 'mean'
                ax.axvline(x=np.mean(x), color='k', ls='-', lw=1, label=extralabel)
                if ii == 0:
                    extralabel = 'median'
                ax.axvline(x=np.median(x), color='k', ls='--', label=extralabel)
            except ValueError as e:
                print('x', x)
                print(e)
            if ii == 0:
                ax.text(0.05, 0.05, 'n={:d}/{:d}'.format(len(x), df.index[-1]), ha='left', va='bottom',
                        transform=ax.transAxes, color='xkcd:tangerine')
        except KeyError:
            print('Key', key, 'not found in', case)

    ax.yaxis.set_ticks([])
    if legend:
        ax.legend(frameon=False, fontsize=labelsize - 4)
    if setxlabel:
        ax.set_xlabel(xlabel, fontsize=labelsize)
    if settitle:
        ax.set_title(case, fontsize=labelsize)
    if save:
        plot_save(fig, case + fend, fig_path=fig_path, fig_fmt=fig_fmt)
    return fig, ax


def subplots_evol_at_sol(Ra_ls, eta_ls, regime_grid=None, save=True, t1_grid=None,
                         load_grid='auto', psuffixes=None,
                         fig_path=fig_path_bullard, fname='evol', fig_fmt='.png', end_grid=None, normtime=True,
                         labelsize=14, xlabel=r'Time', ylabels=None, keys=None, title='', legsize=10,
                         xlabelpad=8, ylabelpad=-2, markers=None, markersize=20, alpha=0.5,
                         fig=None, axes=None, cmap='magma', vmin=None, vmax=None, include_regimes=None,
                         regime_names=None, colour_by='Ra', zoom=True,
                         data_path=data_path_bullard, **kwargs):
    # plot time-evolution of list of keys for all cases in given regime

    if psuffixes is None:
        psuffixes = ['_T']
    Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid) = pro.reshape_inputs(Ra_ls, eta_ls, (
        t1_grid, load_grid, end_grid, regime_grid))
    if include_regimes is None:
        include_regimes = regime_names
    if markers is None:
        markers = ['o', '^', 's', 'D', 'v', 'X']
    if ylabels is None:
        ylabels = keys
    if iterable_not_string(keys):
        nkeys = len(keys)
    elif keys is None:
        raise Exception('No y-axis keys provided!')
    else:
        nkeys = 1
    if fig is None:
        fig, axes = plt.subplots(nkeys, 1, figsize=(7, nkeys * 2), sharex=True)
        if nkeys == 1:
            axes = np.array([axes])
    logeta_fl = [np.log10(float(a)) for a in eta_ls]
    logRa_fl = [np.log10(float(a)) for a in Ra_ls]
    if colour_by == 'eta':
        c_list = colorize(logeta_fl, cmap=cmap, vmin=vmin, vmax=vmax)[0]
    elif colour_by == 'Ra':
        c_list = colorize(logRa_fl, cmap=cmap, vmin=vmin, vmax=vmax)[0]

    mins = np.ones(len(keys))
    maxes = np.zeros(len(keys))
    for jj, eta_str in enumerate(eta_ls):
        cases, Ra_var = pro.get_cases_list(Ra_ls, eta_str, end_grid[jj])
        if colour_by == 'eta':
            c = c_list[jj]

        for ii, case in enumerate(cases):
            t1_ii = t1_grid[jj][ii]
            load_ii = load_grid[jj][ii]
            marker_ii = markers[ii]
            if colour_by == 'Ra':
                c = c_list[ii]

            if (t1_ii != 1) and (os.path.exists(data_path + 'output-' + case)) and (
                    regime_grid[jj][ii] in include_regimes):
                Ra_ii = float(Ra_var[ii])

                # load data
                if load_ii == 'auto':
                    dat = ad.Aspect_Data(directory=data_path + 'output-' + case + '/',
                                         read_statistics=False, read_parameters=False, **kwargs)
                else:
                    dat = None
                dfs = []
                for ip, suffix in enumerate(psuffixes):
                    df1 = pro.pickleio(case, suffix=suffix, t1=t1_ii,
                                       dat_new=dat, data_path=data_path, load=load_ii, **kwargs)
                    dfs.append(df1)
                try:
                    df = pd.concat(dfs, axis=1)
                    df = df.loc[:, ~df.columns.duplicated()]
                except Exception as e:
                    for dfi in dfs:
                        print(dfi)
                    raise e

                # do the plotting on each axis
                x_data = df['time']
                if normtime:
                    # Normalised [0,1]
                    x_data = (x_data - np.min(x_data)) / np.ptp(x_data)
                for k, key in enumerate(keys):
                    ax = axes[k]
                    y_data = df[key]
                    ax.set_ylabel(ylabels[k], fontsize=labelsize, labelpad=ylabelpad)
                    if k == len(keys) - 1:
                        ax.set_xlabel(xlabel, fontsize=labelsize, labelpad=xlabelpad)
                    ax.scatter(x_data, y_data, color=c, s=markersize, marker=marker_ii, alpha=alpha)
                    ax.plot(x_data, y_data, color=c, lw=0.8, alpha=alpha)
                    if np.min(y_data) < mins[k]:
                        mins[k] = np.min(y_data)
                    if np.max(y_data) > maxes[k]:
                        maxes[k] = np.max(y_data)

    # legend proxy artist
    ax = axes[0]
    lines = []
    for jj, leta in enumerate(logeta_fl):
        if colour_by == 'eta':
            c = c_list[jj]
            marker = 'o'
        elif colour_by == 'Ra':
            c = 'k'
            marker = markers[jj]
        p = mlines.Line2D([], [], lw=0, color=c, marker=marker, markersize=markersize / 4, label=leta,
                          alpha=alpha)
        lines.append(p)
    legend1 = ax.legend(lines, [l.get_label() for l in lines], fontsize=legsize, frameon=True, loc="upper left",
                        title=r"log $\Delta \eta$", )
    ax.add_artist(legend1)

    ax = axes[-1]
    lines = []
    for ii, Ra in enumerate(Ra_ls):
        if colour_by == 'eta':
            c = 'k'
            marker = markers[ii]
        elif colour_by == 'Ra':
            c = c_list[ii]
            marker = 'o'
        p = mlines.Line2D([], [], lw=0, color=c, marker=marker, markersize=markersize / 4, label=Ra, alpha=alpha)
        lines.append(p)
    legend2 = ax.legend(lines, [l.get_label() for l in lines], fontsize=legsize, frameon=True, loc="lower right",
                        title="Ra", )
    ax.add_artist(legend2)

    if zoom:
        # zoom to data limits
        for k, ax in enumerate(axes):
            ax.set_ylim((mins[k], maxes[k]))

    plt.suptitle(title, fontsize=labelsize, y=1.02)
    if save:
        plot_save(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt)
    return fig, axes


def subplots_hist(Ra_ls, eta_ls, regime_grid=None, save=True, t1_grid=None, nbins=20,
                  load_grid='auto', psuffixes=None,
                  fig_path=fig_path_bullard, fname='hist-evol', fig_fmt='.png', end_grid=None, normtime=True,
                  labelsize=14, xlabels=None, keys=None, title='', legsize=10,
                  xlabelpad=8, alpha=0.5,
                  fig=None, axes=None, cmap='magma', vmin=None, vmax=None, include_regimes=None,
                  regime_names=None, colour_by='eta', overlap=False,
                  data_path=data_path_bullard, **kwargs):
    # plot histograms of time evolution (implemented for single Ra or single eta - marginilization must be opposite
    # to 'colour_by')

    if psuffixes is None:
        psuffixes = ['_T']
    Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid) = pro.reshape_inputs(Ra_ls, eta_ls, (
        t1_grid, load_grid, end_grid, regime_grid))
    if include_regimes is None:
        include_regimes = regime_names
    if xlabels is None:
        xlabels = keys
    if iterable_not_string(keys):
        nkeys = len(keys)
    elif keys is None:
        raise Exception('No x-axis keys provided!')
    else:
        nkeys = 1
    if overlap:
        nrows = 1
    if colour_by == 'Ra':
        nrows = len(Ra_ls)
    elif colour_by == 'eta':
        nrows == len(eta_ls)
    if fig is None:
        fig, axes = plt.subplots(nrows, nkeys, figsize=(nkeys * 2, nrows * 3), sharex='col')
        if nkeys == 1:
            axes = np.array([axes])

    logeta_fl = [np.log10(float(a)) for a in eta_ls]
    logRa_fl = [np.log10(float(a)) for a in Ra_ls]
    if colour_by == 'eta':
        c_list = colorize(logeta_fl, cmap=cmap, vmin=vmin, vmax=vmax)[0]
    elif colour_by == 'Ra':
        c_list = colorize(logRa_fl, cmap=cmap, vmin=vmin, vmax=vmax)[0]

    for jj, eta_str in enumerate(eta_ls):
        cases, Ra_var = pro.get_cases_list(Ra_ls, eta_str, end_grid[jj])

        for ii, case in enumerate(cases):
            t1_ii = t1_grid[jj][ii]
            load_ii = load_grid[jj][ii]
            if colour_by == 'eta':
                c = c_list[jj]
            elif colour_by == 'Ra':
                c = c_list[ii]
            if overlap:
                axs = axes
            else:
                if colour_by == 'eta':
                    axs = axes[jj]
                elif colour_by == 'Ra':
                    axs = axes[ii]

            if (t1_ii != 1) and (os.path.exists(data_path + 'output-' + case)) and (
                    regime_grid[jj][ii] in include_regimes):
                Ra_ii = float(Ra_var[ii])

                # load data
                if load_ii == 'auto':
                    dat = ad.Aspect_Data(directory=data_path + 'output-' + case + '/',
                                         read_statistics=False, read_parameters=False, **kwargs)
                else:
                    dat = None
                dfs = []
                for ip, suffix in enumerate(psuffixes):
                    df1 = pro.pickleio(case, suffix=suffix, t1=t1_ii,
                                       dat_new=dat, data_path=data_path, load=load_ii, **kwargs)
                    dfs.append(df1)
                try:
                    df = pd.concat(dfs, axis=1)
                    df = df.loc[:, ~df.columns.duplicated()]
                except Exception as e:
                    for dfi in dfs:
                        print(dfi)
                    raise e

                # do the plotting on each axis
                for k, key in enumerate(keys):
                    ax = axs[k]
                    data = df[key].dropna()
                    ax.set_xlabel(xlabels[k], fontsize=labelsize, labelpad=xlabelpad)
                    ax.hist(data, histtype='step', bins=nbins, color=c, density=True)
                    ax.axvline(np.mean(data), c=c, ls='--')
                    ax.axvline(np.median(data), c=c, ls='-')
                    ax.text(0.95, 0.95, 'n = ' + str(len(data)), fontsize=10, c='k', horizontalalignment='right',
                            verticalalignment='top', transform=ax.transAxes)

    # legend proxy artist
    if overlap:
        ax = axes[0]
    else:
        ax = axes[0, 0]
    lines = []
    if colour_by == 'eta':
        for jj, leta in enumerate(logeta_fl):
            p = mlines.Line2D([], [], lw=2, color=c_list[jj], label=leta,
                              alpha=alpha)
            lines.append(p)
        legtitle = r"log $\Delta \eta$"
    elif colour_by == 'Ra':
        for ii, lRa in enumerate(logRa_fl):
            p = mlines.Line2D([], [], lw=2, color=c_list[ii], label=lRa,
                              alpha=alpha)
            lines.append(p)
        legtitle = r"log Ra"

    legend1 = ax.legend(lines, [l.get_label() for l in lines], fontsize=legsize, frameon=True, loc="upper left",
                        title=legtitle)
    ax.add_artist(legend1)

    legend2 = ax.legend(handles=[mlines.Line2D([], [], lw=2, ls='--', color='k', label='mean'),
                                 mlines.Line2D([], [], lw=2, ls='-', color='k', label='median')],
                        bbox_to_anchor=(0.05, 1.05), fontsize=legsize, frameon=True, loc="lower left")
    ax.add_artist(legend2)

    plt.suptitle(title, fontsize=labelsize, y=1.02)
    if save:
        plot_save(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt)
    return fig, axes


def plot_evol(case, col, fig=None, ax=None, save=True, fname='_f', mark_used=True, t1=0, dat=None, show_sols=False,
              ylabel='rms velocity', xlabel='time', yscale=1, c='k', settitle=True, setxlabel=True, fig_fmt='.png',
              setylabel=True, legend=False, labelsize=16, labelpad=5, label=None, sol_df=None,
              fig_path=fig_path_bullard):
    if not setxlabel:
        xlabel = ''
    if not setylabel:
        ylabel = ''
    if ax is None:
        fig = plt.figure()
        ax = plt.gca()

    time, y = pro.read_evol(case, col, dat=dat)
    ax.plot(time, y * yscale, c=c, lw=0.5, label=label)
    ax.set_xlim(0, ax.get_xlim()[1])
    ax.set_xlabel(xlabel, fontsize=labelsize, labelpad=labelpad)
    ax.set_ylabel(ylabel, fontsize=labelsize, labelpad=labelpad)
    if settitle:
        ax.set_title(case, fontsize=labelsize)
    if legend:
        ax.legend(frameon=False, fontsize=labelsize - 4)
    if mark_used:
        # Create a Rectangle patch to overlie "transient" times
        rect = patches.Rectangle((ax.get_xlim()[0], ax.get_ylim()[0]), t1,
                                 ax.get_ylim()[1] - ax.get_ylim()[0],
                                 edgecolor='None', facecolor='xkcd:pale olive green', alpha=0.6, zorder=3)
        ax.add_patch(rect)
    if show_sols and sol_df is not None:
        sol_times = np.array(sol_df['time'])  # steady state sols
        # print('sol times', sol_times)
        for t in sol_times:
            ax.axvline(x=t, c='k', lw=0.5, ls='-', alpha=0.6, zorder=0)
    if save:
        plot_save(fig, case + fname, fig_path=fig_path, fig_fmt=fig_fmt)
    return fig, ax


def plot_topo_profile(case, ts, save=True, fig_path=fig_path_bullard, data_path=data_path_bullard, verbose=True,
                      fig_fmt='.png', ):
    x, h = pro.read_topo_stats(case, ts, data_path=data_path)
    # normalize to 0 mean
    h_norm = pro.trapznorm(h)
    fig = plt.figure()
    plt.plot(x, h_norm)
    plt.xlabel('x')
    plt.ylabel('dynamic topography')
    plt.title(case)
    if verbose:
        print('mean:', pro.trapzmean(h_norm))
        print('max:', np.max(h_norm))
        print('min:', np.min(h_norm))
    if save:
        plot_save(fig, case + '_h_' + '{:05}'.format(ts), fig_path=fig_path, fig_fmt=fig_fmt)


def fit_cases_on_plot(yx_all, ax, yerr=1, xerr=1, legend=True, showallscatter=False, n_fitted=2, c_list=None,
                      c='xkcd:periwinkle', legsize=12, lw=1, legloc='lower left', showchisq=False,
                      **kwargs):
    x = [a[1] for a in yx_all]
    y = [a[0] for a in yx_all]

    if np.array(x[0]).ndim > 0 and np.array(y[0]).ndim > 0:
        flatx = [item for sublist in x for item in sublist]
        flaty = [item for sublist in y for item in sublist]
    else:
        flatx, flaty = x, y

    if n_fitted == 3:  # fit to 3 parameter power law
        flatx0 = [a[0] for a in flatx]
        flatx1 = [a[1] for a in flatx]
        flatx = [flatx0, flatx1]

    const, expon, const_err, expon_err, chisqr, MSE = pro.fit_wrapper(flatx, flaty, yerr=yerr, xerr=xerr,
                                                                      n_fitted=n_fitted, **kwargs)

    # newlabel = 'C = {:.2e} +- {:.2e}'.format(const, const_err)
    newlabel = r'$C = {:.2f} \pm {:.2f}$'.format(const, const_err)
    if expon is not None:
        # newlabel = newlabel + '\np = {:.3f} +- {:.3f}'.format(expon[0], expon_err[0])
        newlabel = newlabel + '\n' + r'$p = {:.2f} \pm {:.2f}$'.format(expon[0], expon_err[0])

    if len(expon) > 1:
        # newlabel = newlabel + '\nq = {:.3f} +- {:.3f}'.format(expon[1], expon_err[1])
        newlabel = newlabel + '\n' + r'$q = {:.2f} \pm {:.2f}$'.format(expon[1], expon_err[1])

    # plot
    xprime = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1])
    if n_fitted == 3:
        z_vec = np.unique(flatx1)
        if c_list is None:
            c_list = colorize(np.log10(z_vec), cmap='winter')[0]
        for ind, z in enumerate(z_vec):
            hprime = const * xprime ** expon[0] * z ** expon[1]
            h2, = ax.plot(xprime, hprime, c=c_list[ind], ls='--', lw=lw, zorder=100, label='dum')

    else:
        hprime = const * xprime ** expon
        h3, = ax.plot(xprime, hprime, c=c, ls='--', lw=lw, zorder=100, label='dum')

        # error
        SE_y = pro.fit_SE(flatx, flaty, [np.log10(const), expon], xn=xprime)

        yn_upper = hprime + SE_y
        yn_lower = hprime - SE_y
        ax.fill_between(xprime, yn_lower, yn_upper, fc='k', alpha=0.3)

    if legend:
        if showchisq:
            newlabel = newlabel + '\n' + r'$\chi^2_\nu = {:.4f}$'.format(MSE)
        # handles, labels = ax.get_legend_handles_labels()
        # try:
        #     labels[-1] = newlabel
        # except IndexError:
        #     labels = newlabel
        ax.text(0.95, 0.95, newlabel, fontsize=legsize,
                horizontalalignment='right',
                verticalalignment='top',
                transform=ax.transAxes)
        # leg = ax.legend(fontsize=legsize, handles=handles, labels=labels, loc=legloc)
        # ax.add_artist(leg)

    if showallscatter:
        ax.scatter(flatx, flaty, c=c, alpha=0.05, s=10)
    return ax


def plot_model_data_errorbars(Ra_ls, eta_ls, regime_grid=None, t1_grid=None, load_grid=None, end_grid=None,
                              literature_file=None, ms=30, fig=None, ax=None, figsize=(7, 5), xlabelpad=10,
                              ylabelpad=10, title=None,
                              legend=True, postprocess_kwargs=None, regime_names=None, which_x='h_components',
                              c_contours='k', fc='w', averagescheme=None, ylim=None, which_h='rms',
                              data_path=data_path_bullard, alpha=1,
                              save=True, fname='model-data', labelsize=16, clist=None, vmin=None, vmax=None,
                              cmap='magma', z_name=None, include_regimes=None, show_cbar=True, clabel=None,
                              cticklabels=None, errortype='time',
                              ylabel='Model', xlabel='Data', errs=None, elw=1, ecapsize=5, crot=0, discrete=True,
                              errorsize=9, sigma=2, **kwargs):
    print('errortype = standard to use errorbars are SE of the mean - TODO check log errorbars')

    if averagescheme is None:
        raise Exception('Averaging scheme not implemeted, must be timefirst or timelast')
    if postprocess_kwargs is None:
        postprocess_kwargs = {}
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid) = pro.reshape_inputs(Ra_ls, eta_ls, (
        t1_grid, load_grid, end_grid, regime_grid))
    if include_regimes is None:
        include_regimes = regime_names
    twocomponent = False
    if iterable_not_string(which_x):
        twocomponent = True

    h_data_all = []
    x_data_all = []
    z_data_all = []
    quants = dict.fromkeys(['h_data', 'x'])
    for jj, eta_str in enumerate(eta_ls):
        cases, Ra_var = pro.get_cases_list(Ra_ls, eta_str, end_grid[jj])
        for ii, case in enumerate(cases):
            t1_ii = t1_grid[jj][ii]
            load_ii = load_grid[jj][ii]
            if (t1_ii != 1) and (os.path.exists(data_path + 'output-' + case)) and (
                    regime_grid[jj][ii] in include_regimes):

                # extract x values for plotting
                if twocomponent:
                    x = [plot_getx(Ra_var[ii], eta_str, case=case, which_x=xx, data_path=data_path,
                                   t1=t1_ii, load=load_ii, postprocess_kwargs=postprocess_kwargs,
                                   averagescheme=averagescheme, **kwargs) for xx in which_x]
                else:
                    x, x_times = plot_getx(Ra_var[ii], eta_str, case=case, which_x=which_x, data_path=data_path,
                                           t1=t1_ii, load=load_ii, postprocess_kwargs=postprocess_kwargs,
                                           averagescheme=averagescheme, return_all=True, **kwargs)

                # get the y values, depending on averaging scheme
                h_rms, h_peak, h_rms_times, h_peak_times = plot_geth(case=case, t1=t1_ii, data_path=data_path,
                                                                     averagescheme=averagescheme,
                                                                     postprocess_kwargs=postprocess_kwargs,
                                                                     return_all=True, **kwargs)
                """^ paste"""

                if which_h == 'rms':
                    h = h_rms
                    h_times = h_rms_times
                elif which_h == 'peak':
                    h = h_peak
                    h_times = h_peak_times
                else:
                    raise Exception('Invalid entry for which_h')

                # append to working
                h_data_all.append(h)
                x_data_all.append(x)
                if z_name == 'eta':
                    z_data_all.append(jj)
                elif z_name == 'regime':
                    z_data_all.append(regime_grid[jj][ii])

                d_times = {'h_data': h_times, 'x': x_times}
                qdict = pro.parameter_percentiles(case, df=d_times,
                                                  keys=quants.keys(), plot=False, sigma=sigma)

                for key in quants.keys():
                    if errortype is 'time':
                        try:
                            quants[key] = np.vstack((quants[key], qdict[key]))  # add to array of errors
                        except ValueError:  # haven't added anything yet
                            quants[key] = qdict[key].reshape((1, 3))  # reshape so it works if you just have one row
                    elif errortype is 'standard':
                        SE_mean = np.std(d_times[key]) / np.sqrt(len(d_times[key]))  # todo: log!
                        avg = np.mean(d_times[key])
                        SE_vec = np.array([avg - SE_mean, avg, avg + SE_mean])
                        try:
                            quants[key] = np.vstack((quants[key], SE_vec))  # add to array of errors
                        except ValueError:  # haven't added anything yet
                            quants[key] = SE_vec.reshape((1, 3))  # reshape so it works if you just have one row

    if twocomponent:
        x0 = [a[0] for a in x_data_all]
        x1 = [a[1] for a in x_data_all]
        h_data = h_data_all
        const, expon, const_err, expon_err, chisqr, MSE = pro.fit_wrapper([x0, x1], h_data, n_fitted=3, **kwargs)
        h_fit = const * np.array(x0) ** expon[0] * np.array(x1) ** expon[1]
        fiterr = None
    else:
        x_data, h_data = [list(tup) for tup in zip(*sorted(zip(x_data_all, h_data_all)))]  # sort according to x
        n_fitted = 2
        if which_x == 'h_components':
            n_fitted = 1
        const, expon, const_err, expon_err, chisqr, MSE = pro.fit_wrapper(x_data, h_data, n_fitted=n_fitted, **kwargs)
        expon = expon[0]
        expon_err = expon_err[0]
        h_fit = const * np.array(x_data) ** expon
        fiterr = pro.fit_SE(x_data, h_data, [np.log10(const), expon], xn=x_data)
    if z_name is None:
        clist = ['k'] * len(h_data)
        z_vec = range(len(h_data))
    else:
        if twocomponent:
            z_vec = z_data_all
        else:
            z_vec = [q for _, q in sorted(zip(x_data_all, z_data_all))]
    ax.set_ylabel(ylabel, fontsize=labelsize, labelpad=ylabelpad)
    ax.set_xlabel(xlabel, fontsize=labelsize, labelpad=xlabelpad)
    if title is None:
        if twocomponent:
            title = 'Fit to h = ({:.2e}'.format(const) + r') Ra' + '^{:.3f}'.format(
                expon[0]) + r' $\Delta \eta$' + '^{:.3f}'.format(expon[1])
        elif which_x == 'h_components':
            title = 'Fit to h = ({:.2f}'.format(const) + r') $\alpha \Delta T_{rh} \delta_{rh}$' + '^{:.3f}'.format(
                expon)
        elif 'Ra' in which_x:
            title = 'Fit to h = ({:.2f}'.format(const) + ') Ra^{:.3f}'.format(expon)
    ax.set_title(title, fontsize=labelsize)

    if clist is None:
        try:
            clist = colorize(z_vec, cmap=cmap, vmin=vmin, vmax=vmax)[0]
        except Exception as e:
            cmap = cmap_from_ascii(cmap, path=cmap_path, end='.txt', ncol=4)
            clist = colorize(z_vec, cmap=cmap, vmin=vmin, vmax=vmax)[0]

    # get errorbars and plot them
    err = dict.fromkeys(quants.keys())
    try:
        for key in quants.keys():
            err[key] = np.asarray([quants[key][:, 1] - quants[key][:, 0], quants[key][:, 2] - quants[key][:, 1]])
        mark = 'o'
        for pp in range(len(h_data)):
            ax.errorbar(h_data[pp], h_fit[pp], yerr=fiterr[pp],
                        xerr=np.asarray([err['h_data'][:, pp]]).T, elinewidth=elw, alpha=alpha,
                        fmt=mark, c=clist[z_vec[pp]], capsize=ecapsize, ms=ms, zorder=10)

    except TypeError as e:
        if quants['h_data'] is None:  # no cases in given regimes as quants is dict of None
            pass
        else:
            raise e

    ax.set_xscale('log')
    ax.set_yscale('log')
    if ylim is not None:
        ax.set_ylim(ylim)
        ax.set_xlim(ylim)
    else:
        ax.axis('equal')
    fig, ax = plot_error_contours(fig, ax, errs=errs, c=c_contours, fc=fc, fontsize=errorsize)
    if save:
        plot_save(fig, fname, **kwargs)
    return fig, ax


def plot_model_data(Ra_ls, eta_ls, regime_grid=None, t1_grid=None, load_grid=None, end_grid=None,
                    literature_file=None, ms=30, fig=None, ax=None, figsize=(7, 5), xlabelpad=10, ylabelpad=10,
                    legend=True, postprocess_kwargs=None, regime_names=None, which_x='h_components',
                    c_contours='k', fc='w', averagescheme=None, ylim=None, which_h='rms', data_path=data_path_bullard,
                    save=True, fname='model-data', labelsize=16, clist=None, vmin=None, vmax=None,
                    cmap='magma', cbar=None, include_regimes=None, show_cbar=True, **kwargs):
    # outdated.........................
    print('WARNING: outdated method!!!!')
    if averagescheme is None:
        raise Exception('Averaging scheme not implemeted, must be timefirst or timelast')
    if postprocess_kwargs is None:
        postprocess_kwargs = {}
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid) = pro.reshape_inputs(Ra_ls, eta_ls, (
        t1_grid, load_grid, end_grid, regime_grid))
    if include_regimes is None:
        include_regimes = regime_names
    twocomponent = False
    if iterable_not_string(which_x):
        twocomponent = True
    if cbar == 'eta':
        clabel = r'$\Delta \eta$'
        cticklabels = None
        if vmin is None:
            vmin = 0.9e5
        if vmax is None:
            vmax = 2e8
        crot = 0
        cnorm = LogNorm(vmin=vmin, vmax=vmax)
        discrete = False
    elif cbar == 'regime':
        clabel = ''  # 'Stationarity'
        cticklabels = ['steady', 'transitional', 'chaotic']
        vmin, vmax = 1, 3
        cnorm = None
        crot = 0  # 70
        cmap = cmap_from_list(clist, cmap_name='regimes')
        discrete = True
    else:
        cnorm = None
        vmin, vmax = None, None

    h_data_all = []
    x_data_all = []
    c_data_all = []
    jj_all = []
    for jj, eta_str in enumerate(eta_ls):
        cases, Ra_var = pro.get_cases_list(Ra_ls, eta_str, end_grid[jj])
        for ii, case in enumerate(cases):
            t1_ii = t1_grid[jj][ii]
            load_ii = load_grid[jj][ii]
            if (t1_ii != 1) and (os.path.exists(data_path + 'output-' + case)) and (
                    regime_grid[jj][ii] in include_regimes):

                # extract x values for plotting
                if twocomponent:
                    x = [plot_getx(Ra_var[ii], eta_str, case=case, which_x=xx, data_path=data_path,
                                   t1=t1_ii, load=load_ii, postprocess_kwargs=postprocess_kwargs,
                                   averagescheme=averagescheme, **kwargs) for xx in which_x]
                else:
                    x = plot_getx(Ra_var[ii], eta_str, case=case, which_x=which_x, data_path=data_path,
                                  t1=t1_ii, load=load_ii, postprocess_kwargs=postprocess_kwargs,
                                  averagescheme=averagescheme, **kwargs)

                # get the y values, depending on averaging scheme
                h_rms, h_peak = plot_geth(case=case, t1=t1_ii, data_path=data_path, averagescheme=averagescheme,
                                          postprocess_kwargs=postprocess_kwargs, **kwargs)
                """^ paste"""

                if which_h == 'rms':
                    h = h_rms
                elif which_h == 'peak':
                    h = h_peak
                else:
                    raise Exception('Invalid entry for which_h')

                h_data_all.append(h)
                x_data_all.append(x)
                jj_all.append(jj)
                if cbar == 'eta':
                    c_data_all.append(float(eta_str))
                elif cbar == 'regime':
                    c_data_all.append(pro.regime_to_digital(ii, jj, regime_grid=regime_grid, regime_names=regime_names))

    if twocomponent:
        x0 = [a[0] for a in x_data_all]
        x1 = [a[1] for a in x_data_all]
        h_data = h_data_all
        expon, const = pro.fit_2log(x1=x0, x2=x1, h=h_data_all)
        h_fit = const * np.array(x0) ** expon[0] * np.array(x1) ** expon[1]

    else:
        x_data, h_data = [list(tup) for tup in zip(*sorted(zip(x_data_all, h_data_all)))]  # sort according to x
        print('x_data: min', np.min(x_data), 'max', np.max(x_data), '(n =', len(x_data), ')')
        expon, const = pro.fit_log(x_data, h_data, weights=None, **kwargs)
        xprime = np.linspace(np.min(x_data), np.max(x_data))
        h_fit = const * np.array(x_data) ** expon
        print('(data, fit)', [(da, fi) for (da, fi) in zip(h_data, h_fit)])
    if not cbar:
        c = 'k'
        vmin, vmax = None, None
    else:
        if twocomponent:
            c = c_data_all
        else:
            c = [q for _, q in sorted(zip(x_data_all, c_data_all))]
    ax.set_ylabel('Model', fontsize=labelsize, labelpad=ylabelpad)
    ax.set_xlabel('Data', fontsize=labelsize, labelpad=xlabelpad)
    if twocomponent:
        title = 'Fit to h = ({:.2e}'.format(const) + r') Ra' + '^{:.3f}'.format(
            expon[0]) + r' $\Delta \eta$' + '^{:.3f}'.format(expon[1])
    elif which_x == 'h_components':
        title = 'Fit to h = ({:.2f}'.format(const) + r') $\alpha \Delta T_{rh} \delta_{rh}$' + '^{:.3f}'.format(expon)
    elif 'Ra' in which_x:
        title = 'Fit to h = ({:.2f}'.format(const) + ') Ra^{:.3f}'.format(expon)
    ax.set_title(title, fontsize=labelsize)
    ax.set_xscale('log')
    ax.set_yscale('log')
    if ylim is not None:
        ax.set_ylim(ylim)
        ax.set_xlim(ylim)
    else:
        ax.axis('equal')

    if (clist is not None) and (not show_cbar):
        for pp in range(len(h_data)):
            scat = ax.scatter(h_data[pp], h_fit[pp], s=ms, zorder=100, c=clist[jj_all[pp]])
    else:
        scat = ax.scatter(h_data, h_fit, s=ms, zorder=100, c=np.squeeze(c), cmap=cmap, norm=cnorm, vmin=vmin, vmax=vmax)

    if show_cbar:
        cbar = colourbar(scat, label=clabel, ticklabels=cticklabels, labelsize=labelsize, discrete=discrete,
                         vmin=vmin, vmax=vmax, rot=crot)

    fig, ax = plot_error_contours(fig, ax, c=c_contours, fc=fc)

    if save:
        plot_save(fig, fname, **kwargs)
    return fig, ax


def plot_error_contours(fig, ax, errs=None, c='k', fc='w', fontsize=9, labels=True, log=True):
    if errs is None:
        errs = [0.5, 0.2, 0.1]
    x0 = np.array(ax.get_xlim())
    y0 = np.array(ax.get_ylim())
    # set 1:1 line
    ax.plot(x0, y0, c=c, lw=2)
    for err in errs:
        y_plus = y0 + err * y0
        y_minus = y0 - err * y0
        for i, y in enumerate([y_plus, y_minus]):
            l, = ax.plot(x0, y, c=c, lw=1, ls='--')
            if labels:
                fY = interp1d(x0, y)
                fX = interp1d(y, x0)
                if log:
                    if i == 0:  # top error
                        pos = [10 ** ((np.log10(x0[0]) + np.log10(fX(y0[-1]))) / 2.),
                               10 ** ((np.log10(fY(x0[0])) + np.log10(y0[-1])) / 2.)]
                    elif i == 1:
                        pos = [10 ** ((np.log10(fX(y0[0])) + np.log10(x0[-1])) / 2.),
                               10 ** ((np.log10(y0[0]) + np.log10(fY(x0[-1]))) / 2.)]
                else:
                    pos = [(x0[-2] + x0[-1]) / 3., (y[-2] + y[-1]) / 3.]

                # transform data points to screen space
                xscreen = ax.transData.transform(np.array((x0[-2::], y[-2::])))
                rot = np.rad2deg(np.arctan2(*np.abs(np.gradient(xscreen)[0][0][::-1])))
                # if (x0[0] < pos[0] < x0[1]) and (y0[0] < pos[1] < y0[1]):
                ltex = ax.text(pos[0], pos[1], '{0:.0f}%'.format(err * 100), size=fontsize, rotation=rot,
                               color=l.get_color(),
                               ha="center", va="center", bbox=dict(boxstyle='square,pad=-0.0', ec=fc, fc=fc))
    return fig, ax


def plot_norm_spectra(Ra_ls, eta_ls, cmap='rainbow', end_grid=None, regime_grid=None, include_regimes=None,
                      data_path=data_path_bullard, pend='_sph', fend='.pkl', figname='h_spectra_stacked',
                      marker='.', lw=0.5, xlabel='Wavenumber', ylabel='Normalised power spectral density',
                      x2label='spherical harmonic degree', save=True, alpha=1, labelsize=16, ticksize=12,
                      fig=None, ax=None, figsize=(5, 5), z_name='Ra', vmin=None, vmax=None, clabel=None,
                      norm='min_l', dim=False, d=1, alpha_m=1, dT=1, R_p=None, cbar=False, show_degrees=True,
                      add_files=None, add_label=None, legsize=12, xlim=None, ylim=None, show_beta_guide=False,**kwargs):
    import pickle as pkl
    import sh_things as sh
    global R
    if R_p is None:
        R_p = d

    if fig is None and ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    if clabel is None:
        xlabel = z_name

    # get colouring
    ji_use = []
    z_vec = []
    for jj, eta_str in enumerate(eta_ls):
        cases, Ra_var = pro.get_cases_list(Ra_ls, eta_str, end_grid[jj])
        for ii, case in enumerate(cases):
            if regime_grid[jj][ii] in include_regimes:
                fname = data_path + 'output-' + case + '/pickle/' + case + pend + fend
                if os.path.isfile(fname):
                    ji_use.append((jj, ii))
                    if z_name == 'Ra':
                        z_vec.append(np.log10(float(Ra_ls[ii])))
                    elif z_name == 'eta':
                        z_vec.append(np.log10(float(eta_str)))
                    elif z_name == 'Ra_i_eff':
                        z_vec.append(np.log10(
                            plot_getx(Ra_ls[ii], eta_str, case=case, which_x='Ra_i_eff', averagescheme='timefirst',
                                      data_path=data_path, load=True, postprocess_kwargs=None, return_all=False, t1=0,
                                      alpha_m=alpha_m, **kwargs)))
                else:
                    print(fname, 'not found')
    clist = colorize(z_vec, cmap=cmap, vmin=vmin, vmax=vmax)[0]

    # load spectra
    zz = 0
    for (jj, ii) in ji_use:
        case = 'Ra' + Ra_ls[ii] + '-eta' + eta_ls[jj] + end_grid[jj, ii]
        print('case', case)
        fname = data_path + 'output-' + case + '/pickle/' + case + pend + fend

        S, k = pkl.load(open(fname, "rb"))
        if k[0] == 0:  # only wavenumbers greater than 0
            k = k[1:]
            S = S[1:]
        # l = sh.k_to_l(k, R=d)  # should be l=1.9674 at the top

        if dim:
            k = k * d**-1
            S = S * d**3 * dT**2 * alpha_m**2

        # wavenumber range where spectrum makes sense
        wl_min, wl_max = sh.nat_scales(case, ax=None, alpha=alpha_m, d=d, dim=dim, data_path=data_path, **kwargs)
        # wl_min, wl_max = sh.nat_scales(case, dim=False, alpha=alpha_m, data_path=data_path, ax=None, **kwargs)
        k_min, k_max = 2*np.pi / wl_max, 2*np.pi / wl_min
        print('k max (dimensional)', k_max)
        if k_min is not None and (k_min > np.min(k)):
            i_min = np.argmax(k >= k_min)
        if k_max is not None and (k_max < np.max(k)):
            i_max = np.argmax(k >= k_max)
        try:
            kv = k[i_min:i_max + 1]
            Sv = S[i_min:i_max + 1]
        except UnboundLocalError:
            raise Exception('kmin, kmax out of bounds')
        kv, Sv_norm = sh.norm_spectrum(kv, Sv, k_min=k_min, norm=norm)

        ax.plot(kv, Sv_norm, c=clist[zz], alpha=alpha, lw=lw, marker=marker)
        # ax.plot(k, S_norm, c=clist[zz], alpha=0.1, lw=lw, marker=marker)  # not in range

        zz = zz + 1

    ax.loglog()
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_ylabel(ylabel, fontsize=labelsize)
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    if show_beta_guide:
        ax = sh.show_beta_guide(ax, x0=kv[1], y0=Sv_norm[1] * 1e-3, x1=kv[3], m=-2, c='xkcd:slate', lw=3,
                                legsize=legsize)

    if cbar:
        if cmap is None:
            print('cbar not implemented without cmap')
        else:
            cax = colourbar(vector=z_vec, ax=ax, vmin=vmin, vmax=vmax, label=clabel, labelsize=labelsize,
                            ticksize=ticksize, ticks=None, ticklabels=None, labelpad=17,
                            rot=None, discrete=False, cmap=cmap, tickformatter=None, pad=0.05, log=False)
    if show_degrees:
        if dim:
            R = R_p
        else:
            R = 1

        def to_deg(k):
            return k * 2 * np.pi * R - 0.5

        def to_wn(l):
            return (l + 0.5) / (2 * np.pi * R)

        secax = ax.secondary_xaxis('top', functions=(to_deg, to_wn))
        secax.set_xlabel(x2label, fontsize=labelsize)
        secax.tick_params(axis='both', which='major', labelsize=ticksize)

    if add_files is not None:
        for ii, f in enumerate(add_files):
            ax = plot_from_txt(f, ax, label=add_label[ii], header=0,
                               additional_mod_fn=sh.mod_loaded_spectrum,
                               plot_kwargs={'c': 'k'}, is_2D=True, is_wl=True, normalise=True, norm=norm)
        ax.legend(frameon=False, fontsize=legsize)

    if save:
        plot_save(fig, fname=figname, **kwargs)
    return fig, ax


def plot_from_txt(filepath, ax, label=None, header=0, additional_mod_fn=None, plot_kwargs=None, **kwargs):
    if plot_kwargs is None:
        plot_kwargs = {}
    df = pd.read_csv(filepath, header=header, names=['x', 'y'], index_col=False, comment='#')
    x, y = df['x'].to_numpy(), df['y'].to_numpy()
    if additional_mod_fn is not None:
        x, y = additional_mod_fn(x, y, **kwargs)
    ax.plot(x, y, label=label, **plot_kwargs)
    return ax
