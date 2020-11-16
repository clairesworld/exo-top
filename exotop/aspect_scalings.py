import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.image as mpimg
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
import pickle as pkl
from collections import Iterable
from six import string_types
from scipy.optimize import curve_fit
from scipy import stats

sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop import aspect_postprocessing2 as post
from exotop.mpl_tools import colorize

data_path_bullard = '/raid1/cmg76/aspect/model-output/'
fig_path_bullard = '/raid1/cmg76/aspect/figs/'

def savefig(fig, fname, fig_path=fig_path_bullard, fig_fmt='.png', bbox_inches='tight', **kwargs):
    path = fig_path + fname + fig_fmt
    directory = os.path.dirname(path)
    os.makedirs(directory, exist_ok=True)
    fig.savefig(path, bbox_inches=bbox_inches, **kwargs)
    print('Saved to ', path, '!')

def read_topo_stats(case, ts, data_path=data_path_bullard):
    df = pd.read_csv(data_path + 'output-' + case + '/dynamic_topography_surface.' + '{:05}'.format(ts), header=None,
                     names=['x', 'y', 'h'], skiprows=1, index_col=False, delimiter=r"\s+", engine='python')
    return df['x'], df['h']


def pickleio(case, suffix, postprocess_functions, t1=0, load='auto', dat_new=None, at_sol=True,
             data_path=data_path_bullard, fend='.pkl', **kwargs):
    # do pickling strategy
    case_path = data_path + 'output-' + case + '/'
    fname = case + suffix + fend
    df = pd.DataFrame()
    dump_flag = False
    reprocess_flag = False
    t1_new = 0

    if os.path.exists(case_path):  # do nothing if case doesn't exist
        os.makedirs(case_path + 'pickle/', exist_ok=True)
        if load == 'auto' or load is True:
            if os.path.exists(case_path + 'pickle/' + fname):
                # open pickled file
                try:
                    df = pkl.load(open(case_path + 'pickle/' + fname, "rb"))
                except ValueError:  # python 2?
                    df = pkl.load(open(case_path + 'pickle/' + fname, "rb"), protocol=2)
                print('    Found', fname)

                if load == 'auto':  # check for additional timesteps
                    print('      Checking for new timesteps...')
                    time_f_old = df.time.iat[-1]
                    if dat_new is None:
                        dat_new = post.Aspect_Data(directory=case_path, verbose=False,
                                                   read_statistics=True, read_parameters=False)
                    time_new = dat_new.stats_time
                    t1_new = time_new[np.argmin(time_new > time_f_old)]
                    if t1_new > 0:  # new timesteps
                        reprocess_flag = True
                        print('      Updating', fname)

            elif load == 'auto':  # pkl file not found
                reprocess_flag = True
                dat_new = post.Aspect_Data(directory=case_path, verbose=False,
                                           read_statistics=True, read_parameters=False)
                print('    File', fname, 'not found, processing...')

        else:  # load is False so automatically calculate shit
            reprocess_flag = True
            dat_new = post.Aspect_Data(directory=case_path, verbose=False,
                                       read_statistics=True, read_parameters=False)

        if reprocess_flag:
            if at_sol:
                sol_files_new = dat_new.read_stats_sol_files()
                df = process_at_solutions(case, postprocess_functions=postprocess_functions, dat=dat_new,
                                          t1=np.maximum(t1, t1_new), data_path=data_path, sol_files=sol_files_new,
                                          df_to_extend=df, **kwargs)
            else:
                df = process_steadystate(case, postprocess_functions=postprocess_functions, dat=dat_new,
                                         t1=np.maximum(t1, t1_new), data_path=data_path,
                                         df_to_extend=df, **kwargs)
            dump_flag = True  # always save if you did something

        if dump_flag:
            pkl.dump(df, open(case_path + 'pickle/' + fname, "wb"))

    return df


def concat_pickles(case, keys=None, suffixes=None, new_suffix=None, fend='.pkl', data_path=data_path_bullard):
    case_path = data_path + 'output-' + case + '/'
    dfs = []
    for suffix in suffixes:
        fname = case + suffix + fend
        if os.path.exists(case_path + 'pickle/' + fname):
            df_new = pd.DataFrame()
            df_loaded = pkl.load(open(case_path + 'pickle/' + fname, "rb"))  # open pickled file
            for key in keys:
                if key in df_loaded.columns:
                    df_new[key] = df_loaded[key]  # add this column to new df
                    print('Copied column', key, 'from', fname)
            dfs.append(df_new)
        else:
            print('File', fname, 'does not exist')
    pkl.dump(pd.concat(dfs), open(case_path + 'pickle/' + case + new_suffix + fend, "wb"))




def get_cases_list(Ra, eta):
    # either Ra or eta is iterable
    if isinstance(Ra, Iterable) and not isinstance(Ra, string_types):
        x_var = Ra
        cases = ['Ra' + r + '-eta' + eta + '-wide' for r in Ra]
    elif isinstance(eta, Iterable) and not isinstance(eta, string_types):
        x_var = eta
        cases = ['Ra' + Ra + '-eta' + e + '-wide' for e in eta]
    else:
        raise Exception('Ra or eta must be iterable')
    return cases, x_var


def trapznorm(A):
    mean = np.trapz(A) / (len(A) - 1)
    #     print('original mean:',mean)
    return A - mean


def trapzmean(A):
    return np.trapz(A) / (len(A) - 1)


def peak_and_rms(h):
    return np.max(h), np.sqrt(trapzmean(h ** 2))


def read_evol(case, col, dat=None, data_path=data_path_bullard):
    # return time, column i, (from statistics file)
    if dat is None:
        dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', read_statistics=True, verbose=False)
    return dat.stats_time, eval('dat.stats_' + col)


# def read_evol_old(case, i, data_path=data_path_bullard, skiprows=None):
#     # e.g. return time, column i, nsteps (from statistics file)
#     df = pd.read_csv(data_path + 'output-' + case + '/statistics', header=None,
#                      skiprows=skiprows,
#                      index_col=False, delimiter=r"\s+", engine='python', comment='#')
#     return np.array(df.iloc[:, 1]), np.array(df.iloc[:, i - 1]), len(df.index)


def process_at_solutions(case, postprocess_functions, dat=None, t1=0, data_path=data_path_bullard,
                         df_to_extend=None, sol_files=None, **kwargs):
    if df_to_extend is None:
        df_to_extend = pd.DataFrame()
    if dat is None:
        dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False, read_statistics=True,
                               read_parameters=False)
    time = dat.stats_time
    i_time = np.argmax(time > t1)  # index of first timestep to process
    if i_time > 0:
        if sol_files is None:
            try:
                sol_files = dat.sol_files
            except AttributeError:
                sol_files = dat.read_stats_sol_files()
        sols_in_time = sol_files[i_time:]
        n_quasi, n_indices = np.unique(sols_in_time, return_index=True)  # find graphical snapshots within time range
        n_ts = n_indices + i_time  # TODO: not off by 1 ?
        for ii, n in enumerate(n_quasi):
            ts = n_ts[ii]  # timestep at this solution
            for fn in postprocess_functions:
                new_params_dict = fn(case, n=n, ts=ts, dat=dat, **kwargs)
                new_params_dict['sol'] = n
                new_params_dict['time'] = time[ts]
                new_params = pd.DataFrame(new_params_dict, index=[ts])
                df_to_extend = pd.concat([df_to_extend, new_params])
                print('        Processed', fn, 'for solution', n, '/', int(n_quasi[-1]))

    else:
        new_params = pd.DataFrame({'sol':[None], 'time':[0]})
        df_to_extend = pd.concat([df_to_extend, new_params])
        print('    No timesteps after t =', time[i_time], '(tf =', time[-1], ')')
        # print('    No solutions after t =', time[i_time], '(tf =', time[-1], ')')
    return df_to_extend


def process_steadystate(case, postprocess_functions, dat=None, t1=0, data_path=data_path_bullard,
                        df_to_extend=None, **kwargs):
    if df_to_extend is None:
        df_to_extend = pd.DataFrame()
    if dat is None:
        dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False, read_statistics=True,
                               read_parameters=False)
    time = dat.stats_time
    i_time = np.argmax(time > t1)  # index of first timestep to process

    if (t1 != 1) and (i_time > 0):
        print('        Processing', postprocess_functions, 'for ', case, ', from timestep', i_time, 'to', len(time))
        for ii in range(i_time, len(time)):
            ts = ii  # timestep at this solution
            for fn in postprocess_functions:
                new_params_dict = fn(case, n=None, ts=ts, dat=dat, **kwargs)
                new_params_dict['time'] = time[ts]
                new_params = pd.DataFrame(new_params_dict, index=[ts])
                df_to_extend = pd.concat([df_to_extend, new_params])

    else:
        new_params = pd.DataFrame({'time':[0]})
        df_to_extend = pd.concat([df_to_extend, new_params])
        print('    No timesteps after t =', time[i_time], '(tf =', time[-1], ')')
    return df_to_extend


def h_at_ts(case, ts=None, hscale=1, **kwargs):
    h_params_n = {}
    try:
        x, h = read_topo_stats(case, ts)
        h_norm = trapznorm(h)
        peak, rms = peak_and_rms(h_norm)
        h_params_n['h_peak'] = peak * hscale
        h_params_n['h_rms'] = rms * hscale

    except FileNotFoundError as e:
        print('    No dynamic topography found at ts =', ts)
        h_params_n['h_peak'] = np.nan
        h_params_n['h_rms'] = np.nan

    for key in h_params_n.keys():
        h_params_n[key] = [h_params_n[key]]
    return h_params_n


def Nu_at_ts(case, ts=None, dat=None, data_path=data_path_bullard, **kwargs):
    if dat is None:
        dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False, read_statistics=True)
    Nu = dat.Nu(k=1)
    return {'Nu': Nu[ts]}


def T_parameters_at_sol(case, n, dat=None, data_path=data_path_bullard, **kwargs):
    if dat is None:
        dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False,
                               read_statistics=True, read_parameters=True)
    x, y, z, u, v, _ = dat.read_velocity(n, verbose=False)
    x, y, z, T = dat.read_temperature(n, verbose=False)
    T_params_n = dat.T_components(n, T=T, u=u, v=v, cut=True)
    for key in T_params_n.keys():
        T_params_n[key] = [T_params_n[key]]
    return T_params_n


# def get_T_params_old(case, t1=0, data_path=data_path_bullard, pickleto=None, picklefrom=None, plotTz=False,
#                  setylabel=True, setxlabel=True, savefig=True, verbose=False,
#                  fig_path=fig_path_bullard, fig=None, ax=None,
#                  legend=True, labelsize=16, **kwargs):
#     if (os.data_path.exists(data_path + 'output-' + case)):
#         flag = False
#         if (picklefrom is not None) and (os.data_path.exists(fig_path + 'data/' + picklefrom)):
#             try:
#                 T_params = pkl.load(open(fig_path + 'data/' + picklefrom, "rb"))
#                 dummy = T_params['delta_u']
#                 print('loaded T profiles from', case)
#             except ValueError:
#                 T_params = pkl.load(open(fig_path + 'data/' + picklefrom, "rb"), protocol=2)
#                 print('loaded T profiles from', case)
#             except KeyError:
#                 flag = True
#                 pickleto = picklefrom
#                 print(picklefrom, 'incomplete, re-calculating...')
#         #             if not T_params['dT_rh']: # if empty
#         #                 flag=True
#         else:
#             flag = True
#             if picklefrom is not None:  # save if tried to load but failed
#                 pickleto = picklefrom
#                 print(picklefrom, 'not found, re-calculating...')
#
#         if flag and (t1 != 1):
#             print('calculating T profiles from', case)
#             T_params = {}
#             dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False, read_statistics=True)
#             time = dat.stats_time
#             snaps = dat.read_stats_sol_files()
#             i_time = np.argmax(time > t1)  # index of first timestep in quasi-ss
#             if i_time > 0:
#                 n_quasi = snaps[i_time:]  # find graphical snapshots within time range
#                 for n in np.unique(n_quasi):
#                     n = int(n)
#                     x, y, z, u, v, _ = dat.read_velocity(n, verbose=False)
#                     x, y, z, T = dat.read_temperature(n, verbose=False)
#                     T_params_n = dat.T_components(n, T=T, u=u, v=v, cut=True)
#                     print('    read', n, '/', n_quasi[-1])
#                     for key in T_params_n.keys():
#                         try:
#                             T_params[key].append(T_params_n[key])
#                         except:  # key does not exist yet
#                             T_params[key] = []
#                             T_params[key].append(T_params_n[key])
#             else:
#                 times_at_sols = dat.find_time_at_sol(sol_files=snaps, return_indices=False)
#                 print('no quasi-steady state solutions')
#                 print('  sol files:', snaps)
#                 print('  times at sols:', times_at_sols)
#
#         if pickleto is not None:
#             pkl.dump(T_params, open(fig_path + 'data/' + pickleto, "wb"))
#
#         if plotTz and ((not flag) or (len(n_quasi) > 2)):
#             dT_rh_f = T_params['dT_rh'][-1]
#             delta_rh_f = T_params['delta_rh'][-1]
#             D_l_f = T_params['delta_L'][-1]
#             T_l_f = T_params['T_l'][-1]
#             T_f = T_params['T_av'][-1]
#             fig, ax = dat.plot_profile(T_f, fig=fig, ax=ax, xlabel='', ylabel='', c='k', lw=1)
#             ax.axhline(D_l_f, label='$z_{lid}$', c='xkcd:tangerine', lw=0.5)
#             ax.axhline(D_l_f - delta_rh_f, label=r'$z_\delta$', c='xkcd:red orange', lw=0.5)
#             ax.text(0, D_l_f - delta_rh_f, r'$\delta = $' + '{:04.2f}'.format(delta_rh_f), ha='left', va='top',
#                     color='xkcd:red orange', fontsize=labelsize - 2)
#             ax.plot([T_l_f, T_l_f], [0, D_l_f], ls='--', alpha=0.5, lw=0.5, c='xkcd:tangerine')
#             ax.plot([T_l_f + dT_rh_f, T_l_f + dT_rh_f], [0, D_l_f - delta_rh_f], ls='--', alpha=0.5, lw=0.5,
#                     c='xkcd:red orange')
#             if legend:
#                 ax.legend(frameon=False, fontsize=labelsize - 2)
#             if setxlabel:
#                 ax.set_xlabel('temperature', fontsize=labelsize)
#             if setylabel:
#                 ax.set_ylabel('depth', fontsize=labelsize)
#             if savefig:
#                 fig.savefig(fig_path + case + '-T_z.png', bbox_inches='tight')
#
#         return T_params, fig, ax
#
#     else:
#         print('case', case, 'not found')
#         return {}, fig, ax


# def get_h(case, t1=0, data_path=data_path_bullard, dat=None, load=True,
#              fig_path=fig_path_bullard, hscale=1, fend='_pdtop.pkl', **kwargs):
#     flag = False
#     pfile = fig_path + 'data/' + case + fend
#     if load and (os.data_path.exists(pfile)):
#         try:
#             peak_list, rms_list = pkl.load(open(pfile, "rb"))
#             print('loaded h for case', case)
#         except ValueError:
#             peak_list, rms_list = pkl.load(open(pfile, "rb"), protocol=2)
#             print('loaded h for case', case)
#         if (not peak_list) or (not rms_list):  # if stored stuff is empty
#             flag = True
#     else:
#         flag = True
#
#     if flag:  # load
#         if dat is None:
#             dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False, read_statistics=True,
#                                    read_parameters=False)
#         time = dat.stats_time
#         i_time = np.argmax(time > t1)
#         rms_list = []
#         peak_list = []
#         t_used = []
#         if t1 != 1:
#             print('building distribution of h for case', case, ', n =', len(range(i_time, len(time))))
#             for ii in range(i_time, len(time)):
#                 try:
#                     x, h = read_topo_stats(case, ts=ii)
#                     h_norm = trapznorm(h)
#                     peak, rms = peak_and_rms(h_norm)
#                     rms_list.append(rms)
#                     peak_list.append(peak)
#                     t_used.append(ii)
#                 except FileNotFoundError as e:
#                     print('file not found:', e)
#
#         pkl.dump((peak_list, rms_list), open(pfile, "wb"))
#
#     peak_list_scaled = [a * hscale for a in peak_list]
#     rms_list_scaled = [a * hscale for a in rms_list]
#     return peak_list_scaled, rms_list_scaled


# def pd_top(case, t1=0, data_path=data_path_bullard, fig_path=fig_path_bullard, sigma=2,
#            plot=True, load='auto', peak_list=None, rms_list=None):
#     if sigma == 2:
#         qs = [2.5, 50, 97.5]
#     elif sigma == 1:
#         qs = [16, 50, 84]
#
#     if (peak_list is None) or (rms_list is None):
#         peak_list, rms_list = get_h(case, t1, data_path, pickleto, picklefrom)
#
#     if (not peak_list) or (not rms_list):  # empty
#         print(case, '- h list is empty')
#         return np.array([np.nan, np.nan, np.nan]), np.array([np.nan, np.nan, np.nan]), fig, ax
#
#     if plot:
#         fig, ax = plot_h_pdf(case, rms_list, peak_list, fig_path=fig_path,
#           fig=None, ax=None, savefig=True, settitle=True, setxlabel=True, legend=True,
#            labelsize=16, fend='.png', c_peak='xkcd:forest green', c_rms='xkcd:periwinkle', **kwargs)
#         return np.percentile(peak_list, qs), np.percentile(rms_list, qs), fig, ax
#
#     return np.percentile(peak_list, qs), np.percentile(rms_list, qs)

def parameter_percentiles(case, df=None, keys=None, plot=False, sigma=2, **kwargs):
    # probability distribution for a number of parameters with df containing time evolutions
    if sigma == 2:
        qs = [2.5, 50, 97.5]
    elif sigma == 1:
        qs = [16, 50, 84]
    else:
        print('Unrecognized sigma value')
    qdict = {}
    for key in keys:
        try:
            qdict[key] = np.percentile(df[key], qs)
        except KeyError:
            raise (key, 'not processed yet for', case)
        except Exception as e:
            qdict[key] = np.array([np.nan] * len(qs))
            raise e

    if plot:
        fig, ax = plot_pdf(case, df=df, keys=keys, **kwargs)
        return qdict, fig, ax
    return qdict


def fit_log(x, h):
    x1 = np.log10(np.array(x))  # this should work for time-series of all x corresponding to h
    h1 = np.log10(np.array(h))
    try:
        slope, intercept, r_value, p_value, std_err = stats.linregress(x1, h1)
    except ValueError:
        raise ('error | x', np.shape(x1), 'h', np.shape(h1))
    return slope, 10 ** intercept


def fit_h_sigma(x, h, h_err=None, fn='line'):
    def line(x, a, b):
        return a * x + b

    #     def expon(x, C, n):
    #         return C * x**n

    idx = np.nonzero(np.isnan(x) is False)[0]
    x_fit = np.log10(x[idx])
    h_fit = np.log10(h[idx])
    if h_err is not None:
        h_err = np.log10(h_err[idx])
    print('fitting x =', x_fit, 'h =', h_fit)
    if fn == 'line':
        popt, pcov = curve_fit(line, x_fit, h_fit, sigma=h_err)
        print('slope:', popt[0], 'intercept:', popt[1])

    return 10 ** (popt[1] + popt[0] * x)  # h evaluated at x


def plot_h_vs_Ra(Ra=None, eta=None, t1=None, data_path=data_path_bullard, fig_path=fig_path_bullard,
                 load='auto', showallscatter=False,
                 save=True, fname='h_vs_Ra', sigma=2, fig_fmt='.png',
                 labelsize=16, xlabel='', ylabel='dynamic topography', title='',
                 c_peak='xkcd:forest green', c_rms='xkcd:periwinkle',
                 fit=False, cases=None, x_var=None, logx=True, logy=True,
                 fig=None, ax=None, dt_ylim=(3e-3, 7e-2), xlim=None, hscale=1, **kwargs):
    # Ra or eta is list of strings, t1 is a list of numbers the same length
    if cases is None:
        cases, x_var = get_cases_list(Ra, eta)
    if t1 is None:
        t1 = [0] * len(x_var)

    quants_h_peak = np.zeros((len(x_var), 3))
    quants_h_rms = np.zeros((len(x_var), 3))
    x = np.zeros((len(x_var), 3))
    peak_all = []
    rms_all = []

    for ii, case in enumerate(cases):
        df = pickleio(case, suffix='_h_all', postprocess_functions=['h_at_ts'], t1=t1[ii], load=load, dat_new=None,
                      data_path=data_path, hscale=hscale, at_sol=False, **kwargs)
        qdict = parameter_percentiles(case, df=df, sigma=sigma, keys=['h_peak', 'h_rms'], plot=False)
        quants_h_peak[ii, :] = qdict['h_peak']
        quants_h_rms[ii, :] = qdict['h_rms']

        x[ii, :] = float(x_var[ii])
        h_peak = df['h_peak']
        h_rms = df['h_rms']
        peak_all.append((h_peak, x[ii, :]))
        rms_all.append((h_rms, x[ii, :]))

    yerr_peak = [quants_h_peak[:, 1] - quants_h_peak[:, 0], quants_h_peak[:, 2] - quants_h_peak[:, 1]]
    yerr_rms = [quants_h_rms[:, 1] - quants_h_rms[:, 0], quants_h_rms[:, 2] - quants_h_rms[:, 1]]
    xerr = None

    if ax is None:
        fig = plt.figure()
        ax = plt.gca()

    if fit:
        if len(x_var) > 1:
            fitx = [[a[1]] * len(a[0]) for a in rms_all]
            fith = [a[0] for a in rms_all]
            flatfitx = [item for sublist in fitx for item in sublist]
            flatfith = [item for sublist in fith for item in sublist]
            expon, const = fit_log(flatfitx, flatfith)
            xprime = [a[1] for a in rms_all]
            hprime = const * np.array(xprime) ** expon
            h3, = ax.plot(xprime, hprime, c=c_rms, ls='--', lw=1, zorder=100,
                          label='{:.2e} Ra^{:.3f}'.format(const, expon))
            ax.legend(
                # handles=[h3], labels=[],
                loc='lower left')  # always show what fit is
        else:
            print('not enough points to fit -- Ra', Ra, 'eta', eta)

    ax.errorbar(x[:, 1], quants_h_peak[:, 1], yerr=yerr_peak, xerr=xerr,
                fmt='^', c=c_peak, alpha=0.9, capsize=5)
    ax.errorbar(x[:, 1], quants_h_rms[:, 1], yerr=yerr_rms, xerr=xerr,
                fmt='o', c=c_rms, alpha=0.9, capsize=5)

    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')
    if dt_ylim is not None:
        ax.set_ylim(dt_ylim[0], dt_ylim[1])  # for fair comparison
    if xlim is not None:
        ax.set_xlim(xlim)
    ax.set_ylabel(ylabel, fontsize=labelsize)
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_title(title, fontsize=labelsize)

    if save:
        savefig(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt)
    return fig, ax


def plot_h_vs_Td(Ra=None, eta=None, t1=None, data_path=data_path_bullard, fig_path=fig_path_bullard,
                 load='auto', fig_fmt='.png',
                 save=True, fname='h_T', sigma=2, showallscatter=False,
                 labelsize=16, xlabel=r'$\delta_rh \Delta T_{rh}$', ylabel='dynamic topography', title='',
                 c_peak='xkcd:forest green', c_rms='xkcd:periwinkle', legend=True,
                 fit=False, cases=None, x_var=None, logx=True, logy=True,
                 fig=None, ax=None, dt_ylim=(3e-3, 7e-2), xlim=None, hscale=1, **kwargs):
    # Ra or eta is list of strings, t1 is a list of numbers the same length
    # instead of plotting vs Ra or eta, plot vs theoretical components of scaling relationship

    if cases is None:
        cases, x_var = get_cases_list(Ra, eta)
    if t1 is None:
        t1 = [0] * len(cases)

    quants_h_peak = np.zeros((len(x_var), 3))
    quants_h_rms = np.zeros((len(x_var), 3))
    quants_h_components = np.zeros((len(x_var), 3))
    peak_all = []
    rms_all = []

    for ii, case in enumerate(cases):
        dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False, read_statistics=True)

        # load h and T
        df1 = pickleio(case, suffix='_h', postprocess_functions=['h_at_sol'], t1=t1[ii], load=load, dat_new=dat,
                       data_path=data_path, hscale=hscale, **kwargs)
        df2 = pickleio(case, suffix='_T', postprocess_functions=['T_parameters_at_sol'], t1=t1[ii], load=load,
                       dat_new=dat,
                       data_path=data_path, hscale=hscale, **kwargs)
        df = pd.concat([df1, df2])

        alpha = dat.parameters['Material model']['Simple model']['Thermal expansion coefficient']
        h_components = alpha * (np.array(df['dT_rh']) / np.array(df['dT_m'])) * (
                np.array(df['delta_rh']) / np.array(df['d_m']))
        df['h_components'] = h_components

        qdict = parameter_percentiles(case, df=df, sigma=sigma, keys=['h_peak', 'h_rms', 'h_components'], plot=False)
        quants_h_peak[ii, :] = qdict['h_peak']
        quants_h_rms[ii, :] = qdict['h_rms']
        quants_h_components[ii, :] = qdict['h_components']

        h_peak = df['h_peak']
        h_rms = df['h_rms']

        # old way of accounting for loading all h instead of just at sols
        # t1_idx = np.argmax(dat.stats_time > t1[ii])  # timestep corresponding to t1
        # sol_idx = dat.find_time_at_sol()
        # sol_idx = np.array(sol_idx[sol_idx >= t1_idx])  # cut any values of idx below t1_idx
        # # account for the stored peak_list and rms_list starting at t1 for indexing
        # if np.shape(sol_idx) != np.shape(h_components):
        #     print('inconsistent times: sol_idx', np.shape(sol_idx), 'delta', np.shape(df['delta_rh']))
        # try:
        #     peak_list = [peak_list[j - t1_idx] for j in sol_idx]
        #     rms_list = [rms_list[j - t1_idx] for j in sol_idx]
        # except IndexError:
        #     print('sol_idx - t1_idx', [j - t1_idx for j in sol_idx])

        peak_all.append((h_peak, h_components))
        rms_all.append((h_rms, h_components))

    yerr_peak = [quants_h_peak[:, 1] - quants_h_peak[:, 0], quants_h_peak[:, 2] - quants_h_peak[:, 1]]
    yerr_rms = [quants_h_rms[:, 1] - quants_h_rms[:, 0], quants_h_rms[:, 2] - quants_h_rms[:, 1]]
    xerr = [quants_h_components[:, 1] - quants_h_components[:, 0],
            quants_h_components[:, 2] - quants_h_components[:, 1]]

    if ax is None:
        fig = plt.figure()
        ax = plt.gca()

    fitx = [a[1] for a in rms_all]
    fith_rms = [a[0] for a in rms_all]
    flatfitx = [item for sublist in fitx for item in sublist]
    flatfith_rms = [item for sublist in fith_rms for item in sublist]
    fith_peak = [a[0] for a in peak_all]
    flatfith_peak = [item for sublist in fith_peak for item in sublist]

    if fit:
        #         fitx = x_var # todo: only fit subset?
        #         fitidx = np.where(np.intersect1d(x[:,1], fitx))[0]
        if len(x_var) > 1:  # can only fit if at least 2 data
            expon, const = fit_log(flatfitx, flatfith_rms)
            xprime = np.linspace(np.min(flatfitx), np.max(flatfitx))
            hprime = const * xprime ** expon
            h3, = ax.plot(xprime, hprime, c=c_rms, ls='--', lw=1, zorder=100,
                          label='{:.2e} x^{:.3f}'.format(const, expon))
            ax.legend(
                # handles=[h3], labels=[],
                loc='lower left')  # always show what fit is
        else:
            print('not enough points to fit -- Ra', Ra, 'eta', eta)

    ax.errorbar(quants_h_components[:, 1], quants_h_peak[:, 1], yerr=yerr_peak, xerr=xerr,
                fmt='^', c=c_peak, alpha=0.9, capsize=5)
    ax.errorbar(quants_h_components[:, 1], quants_h_rms[:, 1], yerr=yerr_rms, xerr=xerr,
                fmt='o', c=c_rms, alpha=0.9, capsize=5)

    if showallscatter:
        ax.scatter(flatfitx, flatfith_rms, c=c_rms, alpha=0.1, s=20)
        ax.scatter(flatfitx, flatfith_peak, c=c_peak, alpha=0.1, s=20)

    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')
    if dt_ylim is not None:
        ax.set_ylim(dt_ylim[0], dt_ylim[1])  # for fair comparison
    if xlim is not None:
        ax.set_xlim(xlim)
    ax.set_ylabel(ylabel, fontsize=labelsize)
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_title(title, fontsize=labelsize)

    if save:
        savefig(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt)
    return fig, ax


def subplots_h_vs(Ra_ls, eta_ls, regime_grid, c_regimes, load='auto', save=True,
                  sigma=2, t1=None, fit=False, nrows=2, ncols=2, T_components=False,
                  data_path=data_path_bullard, fig_path=fig_path_bullard, fname='h_Ra_all', fig_fmt='.png',
                  ylim=(6e-3, 7e-2), labelsize=14, xlim=None, xlabel='Ra', ylabel='dynamic topography',
                  logx=True, logy=True, showallscatter=False, xlabelpad=12, ylabelpad=2, hscale=1, **kwargs):
    # subplots for different eta
    #     fig, axes = plt.subplots(2,2, figsize=(7,7))
    #     flaxes = axes.flatten()
    if T_components:
        plot_fn = plot_h_vs_Td
    else:
        plot_fn = plot_h_vs_Ra

    if t1 is None:
        t1 = np.zeros(len(eta_ls), len(Ra_ls))
    fig = plt.figure(figsize=(7, 7))
    bigax = fig.add_subplot(111)  # The big subplot
    bigax.spines['top'].set_color('none')
    bigax.spines['bottom'].set_color('none')
    bigax.spines['left'].set_color('none')
    bigax.spines['right'].set_color('none')
    bigax.tick_params(labelcolor='w', which='both', bottom=False, left=False, right=False, top=False)
    bigax.set_xlabel(xlabel, fontsize=labelsize, labelpad=xlabelpad)
    bigax.set_ylabel(ylabel, fontsize=labelsize, labelpad=ylabelpad)
    legend = False
    for ii, eta in enumerate(eta_ls):
        z = int(str(nrows) + str(ncols) + str(ii + 1))
        ax = fig.add_subplot(z)
        #         ax = flaxes[ii]
        #         if (ii % 2) == 0:
        #             ylabel='dynamic topography'
        #         else:
        #             ylabel=''
        #         if ii==0:
        #             legend=True
        #         else:
        #             legend=False
        #         if ii>1:
        #             xlabel='Ra'
        #         else:
        #             xlabel=''
        Ra_steady = [Ra_ls[j] for j in np.nonzero(regime_grid[ii] == 'steady')[0]]
        Ra_trans = [Ra_ls[j] for j in np.nonzero(regime_grid[ii] == 'transitional')[0]]
        Ra_chaos = [Ra_ls[j] for j in np.nonzero(regime_grid[ii] == 'chaotic')[0]]
        Ra_steady_idx = [j for j in np.nonzero(regime_grid[ii] == 'steady')[0]]
        Ra_trans_idx = [j for j in np.nonzero(regime_grid[ii] == 'transitional')[0]]
        Ra_chaos_idx = [j for j in np.nonzero(regime_grid[ii] == 'chaotic')[0]]

        # steady
        if not (not Ra_steady):
            fig, ax = plot_fn(Ra=Ra_steady, eta=eta, t1=t1[ii, Ra_steady_idx], sigma=sigma, fig=fig, ax=ax,
                              data_path=data_path,
                              c_rms=c_regimes[0], c_peak=c_regimes[0], dt_ylim=ylim, legend=legend,
                              load=load, plotpd=False, xlim=xlim,
                              save=False, fit=True, ylabel='', xlabel='', labelsize=labelsize,
                              fig_path=fig_path, logx=logx, logy=logy,
                              showallscatter=showallscatter, hscale=hscale, **kwargs
                              )
        # trans
        if not (not Ra_trans):
            fig, ax = plot_fn(Ra=Ra_trans, eta=eta, t1=t1[ii, Ra_trans_idx], sigma=sigma, fig=fig, ax=ax,
                              data_path=data_path,
                              c_rms=c_regimes[1], c_peak=c_regimes[1], dt_ylim=ylim, legend=legend,
                              load=load, plotpd=False, xlim=xlim,
                              save=False, fit=True, ylabel='', xlabel='', labelsize=labelsize,
                              fig_path=fig_path, logx=logx, logy=logy,
                              showallscatter=showallscatter, hscale=hscale, **kwargs
                              )
            # chaotic
        if not (not Ra_chaos):
            fig, ax = plot_fn(Ra=Ra_chaos, eta=eta, t1=t1[ii, Ra_chaos_idx], sigma=sigma, fig=fig, ax=ax,
                              data_path=data_path,
                              c_rms=c_regimes[2], c_peak=c_regimes[2], dt_ylim=ylim, legend=legend,
                              load=load, plotpd=False, xlim=xlim,
                              save=False, fit=True, xlabel='', ylabel='', labelsize=labelsize,
                              fig_path=fig_path, logx=logx, logy=logy,
                              showallscatter=showallscatter, hscale=hscale, **kwargs
                              )
        ax.text(0.5, 0.95, '$\Delta \eta$=' + eta, fontsize=labelsize, ha='center', va='top',
                transform=ax.transAxes)
        if ii % ncols != 0:
            ax.yaxis.tick_right()

        # show regime boundaries
    #         try:
    #             ax.axvline(float(Ra_steady[-1])*2, c='k', lw=0.5, alpha=0.6, ls='--')
    #             ax.text(ax.get_xlim()[0], ylim[0], 'steady', fontsize=8, va='bottom', ha='left')
    #         except:
    #             print('no steady regime for eta', eta)
    #         try:
    #             ax.axvline(float(Ra_trans[-1])*2, c='k', lw=0.5, alpha=0.6, ls='--')
    #             ax.text(float(Ra_steady[-1])*2, ylim[0], 'trans.', fontsize=8, va='bottom', ha='left')
    #         except:
    #             print('no transitional regime for eta', eta)
    #         try:
    #             ax.text(float(Ra_chaos[0]), ylim[0], 'chaotic', fontsize=8, va='bottom', ha='left')
    #         except:
    #             print('no chaotic regime for eta', eta)

    axes = fig.axes
    ax = bigax  # axes[0]
    h1 = ax.scatter([], [], label='peak', marker='^', c='k', alpha=0.9)
    h2 = ax.scatter([], [], label='rms', marker='o', c='k', alpha=0.9)
    h3 = ax.scatter([], [], label='steady', marker='o', c=c_regimes[0], alpha=0.9)
    h4 = ax.scatter([], [], label='trans.', marker='o', c=c_regimes[1], alpha=0.9)
    h5 = ax.scatter([], [], label='chaotic', marker='o', c=c_regimes[2], alpha=0.9)
    outer_legend = ax.legend(handles=[h1, h2, h3, h4, h5], labels=['peak', 'rms', 'steady', 'trans.', 'chaotic'],
                             borderaxespad=0., ncol=5,
                             bbox_to_anchor=(0., 1.02, 1., .102), loc=3, frameon=False,
                             # mode="expand",
                             )
    ax.add_artist(outer_legend)
    #     fig.tight_layout()
    fig.subplots_adjust(wspace=0.05, hspace=0.15)
    if save:
        savefig(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt, bbox_inches=None, bbox_extra_artists=(outer_legend,))
    return fig, axes


def scales_with_Ra(Ra_data=None, y_data=None, fig_path=fig_path_bullard,
                   save=True, fname='claire', showallscatter=False,
                   labelsize=16, ylabel='', xlabel='Ra', title='', fig_fmt='.png',
                   c_scatter='xkcd:forest green', legend=True,
                   fit=False, logx=True, logy=True,
                   fig=None, ax=None, ylim=None, xlim=None, **kwargs):
    if fig is None:
        fig = plt.figure()
        ax = plt.gca()

    ax.plot(Ra_data, y_data, '-o', c=c_scatter)

    if fit:
        if (len(Ra_data) > 1):  # can only fit if at least 2 data
            # bl
            expon, const = fit_log(Ra_data, y_data)
            xprime = np.logspace(np.log10(np.min(Ra_data)), np.log10(np.max(Ra_data)))
            yprime = const * xprime ** expon
            h3, = ax.plot(xprime, yprime, c=c_scatter, ls='--', lw=1, zorder=100,
                          label='{:.2e} x^{:.3f}'.format(const, expon))
            ax.legend(
                # handles=[h3], labels=[],
                # loc='lower left'
            )  # always show what fit is
        else:
            print('not enough points to fit -- Ra', Ra_data)

    #     if showallscatter:
    #         ax.scatter(flatfitx, flatfith_rms, c=c_rms, alpha=0.1, s=20)
    #         ax.scatter(flatfitx, flatfith_peak, c=c_peak, alpha=0.1, s=20)

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
        savefig(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt)
    return fig, ax


def plot_multi_Ra_scaling(Ra=None, eta=None, t1=None, keys=None, data_path=data_path_bullard, fig_path=fig_path_bullard,
                          load='auto', save=True, fname='bl-Nu', labelsize=16, ylabel=None, title='',
                          legend=True, cmap='magma', compare_pub=None, compare_label=None, vmin=4, vmax=9,
                          fig=None, axes=None,fig_fmt='.png', **kwargs):
    # Ra or eta is list of strings, t1 is a list of numbers the same length
    # instead of plotting vs Ra or eta, plot vs theoretical components of scaling relationship
    if ylabel is None:
        ylabel = [r'$\delta$', 'Nu']
    if fig is None:
        fig, axes = plt.subplots(2, 1)
    logeta_fl = [np.log10(float(a)) for a in eta]
    c_list = colorize(logeta_fl, cmap=cmap)[0]

    for jj, eta_str in enumerate(eta):
        cases, Ra_var = get_cases_list(Ra, eta_str)
        if t1 is None:
            t1_eta = [0] * len(cases)
        else:
            t1_eta = t1[jj]
        c_scatter = c_list[jj]

        plot_data = {'Ra': []}
        for key in keys:
            plot_data[key] = []

        for ii, case in enumerate(cases):
            t1_ii = t1_eta[ii]

            if (t1_ii != 1) and (os.path.exists(data_path + 'output-' + case)):
                plot_data[Ra].append(float(Ra_var[ii]))

                # load T components  
                dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False,
                                       read_statistics=True)
                df1 = pickleio(case, suffix='_T', postprocess_functions=['T_parameters_at_sol'], t1=t1[ii],
                               load=load, dat_new=dat,
                               data_path=data_path, **kwargs)

                # extract Nu
                df2 = pickleio(case, suffix='_Nu', postprocess_functions=['Nu_at_ts'], t1=t1[ii],
                               load=load, dat_new=dat,
                               data_path=data_path, **kwargs)
                df = pd.concat([df1, df2])

                for key in keys:
                    plot_data[key].append(np.median(df[key]))

                if compare_pub is not None:
                    cmplabel = compare_label
                    if (jj > 0) and (ii > 0):
                        cmplabel = None
                    d_compare = compare_pub(case, dat=dat, Ra=plot_data['Ra'][ii], d_eta=float(eta_str), df=df,
                                            load=load)
                    for k, key in enumerate(keys):
                        try:
                            axes[k].plot(d_compare['Ra_i'], d_compare[key], '^', alpha=0.7, c=c_scatter, label=cmplabel)
                        except KeyError:
                            print('Key', key, 'not returned by', compare_pub)

        for k, key in enumerate(keys):
            fig, axes[0] = scales_with_Ra(Ra_data=plot_data['Ra'], y_data=plot_data[key],
                                          fig_path=fig_path, save=False, labelsize=labelsize, ylabel=ylabel[k],
                                          c_scatter=c_scatter, fig=fig, ax=axes[k], **kwargs)

    scat = axes[-1].scatter(logeta_fl, logeta_fl, visible=False, c=np.array(logeta_fl), cmap=cmap,
                            vmin=vmin, vmax=vmax)  # dummy
    cbar = fig.colorbar(scat, ax=[axes[0], axes[1]])
    cbar.set_label(r'log($\Delta \eta$)', fontsize=labelsize, rotation=270, labelpad=18)

    if save:
        savefig(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt)
    return fig, axes


def solomatov95(Ra=None, d_eta=None, df=None, case=None, dat=None,
                data_path=data_path_bullard, load='auto', **kwargs):
    if df is None:
        df = pickleio(case, postprocess_functions=[T_parameters_at_sol], suffix='_T', dat_new=dat, data_path=data_path,
                      load=load)
    T0 = 1
    dT = 1
    T_i = np.median(df['T_i'])
    gamma = np.log(d_eta)  # gamma for this delta eta
    p = gamma * dT
    eta_0 = np.exp(-gamma * T0)
    eta_i = np.exp(-gamma * T_i)
    Ra_i = np.array(Ra) * eta_0 / eta_i
    delta_0 = 1.85 * p ** 1.3185 * Ra_i ** -0.3185
    delta_1 = p ** -1 * delta_0
    Nu = (delta_0 + delta_1) ** -1
    return {'Ra_i': Ra_i, 'delta_0': delta_0, 'Nu': Nu}


def moresi95(Ra=None, d_eta=None, df=None, dat=None, case=None,
             data_path=data_path_bullard, load='auto', **kwargs):
    if df is None:
        df = pickleio(case, postprocess_functions=[T_parameters_at_sol], suffix='_T', dat_new=dat, data_path=data_path,
                      load=load)
    T0 = 1
    dT = 1
    T_i = np.median(df['T_i'])
    gamma = np.log(d_eta)  # gamma for this delta eta
    p = gamma * dT
    eta_0 = np.exp(-gamma * T0)
    eta_i = np.exp(-gamma * T_i)
    Ra_i = np.array(Ra) * eta_0 / eta_i
    delta_1 = 0.58 * p ** 0.29 * Ra_i ** -0.24
    Nu = 1.89 * p ** -1.02 * Ra_i ** 0.2
    T_i_scaling = 1 - (1.1 * p ** -0.73 * Ra_i ** -0.04)
    delta_0 = T_i_scaling / Nu
    return {'Ra_i': Ra_i, 'delta_0': delta_0, 'Nu': Nu, 'T_i': T_i_scaling, 'delta_1': delta_1}


def subplots_cases(cases, labels=None, labelsize=16, labelpad=5, t1=None, save=True, dt_xlim=(0.0, 0.065),
                   fname='cases', data_path=data_path_bullard, fig_path=fig_path_bullard, fig_fmt='.png',
                   load='auto', includegraphic=False, c_rms='xkcd:forest green', c_peak='xkcd:periwinkle',
                   suptitle='', includepdf=True, includeTz=True, regime_grid=None, show_sols=False, **kwargs):
    # rows are cases, columns are v_rms, q, T(z), hist
    ncases = len(cases)
    ncols = 2
    if regime_grid is None:
        regime_grid = [None]*len(cases)
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
        if os.path.exists(data_path + 'output-' + case):
            print('Plotting summary for', case)
            dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False, read_statistics=True)
            if ii == ncases - 1:  # show x label in bottom row only
                setxlabel = True
            else:
                setxlabel = False
            setylabel = True
            legend = False
            if numplotted == 0:
                legend = True

            if show_sols:  # load df
                sol_df = pickleio(case, suffix='_T', postprocess_functions=[T_parameters_at_sol], t1=t1[ii],
                                  dat_new=dat, load=load, data_path=data_path, fig_path=fig_path, **kwargs)

            ax = axes[ii, icol]
            fig, ax = plot_evol(case, 'rms_velocity', fig=fig, ax=ax, save=False, mark_used=True, t1=t1[ii], dat=dat,
                                show_sols=show_sols, ylabel='rms velocity', c='k', settitle=False, setxlabel=setxlabel,
                                setylabel=setylabel, legend=False, labelsize=labelsize, labelpad=labelpad,
                                sol_df=sol_df)

            ax.text(0.01, 0.95, labels[ii], horizontalalignment='left', verticalalignment='top',
                    transform=ax.transAxes, fontsize=labelsize)

            icol = icol + 1
            ax = axes[ii, icol]
            fig, ax = plot_evol(case, 'heatflux_top', fig=fig, ax=ax, save=False, mark_used=False, dat=dat,
                                show_sols=False, ylabel='heat flux', c='xkcd:light red', settitle=False,
                                setxlabel=setxlabel, setylabel=setylabel, labelsize=labelsize, labelpad=labelpad,
                                label='top')
            fig, ax = plot_evol(case, 'heatflux_bottom', fig=fig, ax=ax, save=False, mark_used=True, t1=t1[ii], dat=dat,
                                show_sols=show_sols, ylabel='heat flux', yscale=-1, c='xkcd:purple blue',
                                settitle=False, setxlabel=setxlabel, setylabel=setylabel, legend=legend,
                                labelsize=labelsize, labelpad=labelpad, label='bottom', sol_df=sol_df)

            if includeTz and t1[ii] < 1:  # final timestep only
                icol = icol + 1
                ax = axes[ii, icol]

                if not show_sols:
                    sol_df = pickleio(case, suffix='_T', postprocess_functions=[T_parameters_at_sol], t1=t1[ii],
                                  dat_new=dat, load=load, data_path=data_path, fig_path=fig_path, **kwargs)

                fig, ax = plot_T_params(case, T_params=sol_df, data_path=data_path, n=-1, save=False,
                                        setxlabel=setxlabel, setylabel=False, legend=False, fig_path=fig_path, fig=fig,
                                        ax=ax)
                if legend:
                    ax.legend(frameon=False)
                if setxlabel:
                    ax.set_xlabel('temperature', fontsize=labelsize)
                if setylabel:
                    ax.set_ylabel('depth', fontsize=labelsize)

            if includepdf and t1[ii] < 1:
                icol = icol + 1
                ax = axes[ii, icol]
                sol_df = pickleio(case, suffix='_h_all', postprocess_functions=[h_at_ts], t1=t1[ii],
                                  dat_new=dat, load=load, data_path=data_path, fig_path=fig_path,
                                  at_sol=False, **kwargs)

                fig, ax = plot_pdf(case, keys=['h_rms', 'h_peak'], df=sol_df, path=data_path,
                                   fig=fig, ax=ax, save=False, settitle=False, setxlabel=setxlabel,
                                   legend=legend, labelsize=labelsize, c_list=[c_rms, c_peak])
                ax.set_xlim(dt_xlim[0], dt_xlim[1])  # for fair comparison

            if includegraphic:
                icol = icol + 1
                ax = axes[ii, icol]
                fgraph = fig_path + 'graphical/' + case + '.png'
                try:
                    img = mpimg.imread(fgraph)
                    ax.imshow(img)
                    print('Plotting graphical output for', case)
                except FileNotFoundError:
                    print('Graphical output not found:', fgraph)
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
        savefig(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt)
    return fig, axes


def plot_convection_regimes(Ra, eta, regime_grid, data_path=data_path_bullard, fig_path=fig_path_bullard, load='auto',
                            save=True, fname='regimes', labelsize=16, fig_fmt='.png', t1=None,
                            overploth=False, nlevels=10, clist=None, cmap_contours='spring', **kwargs):
    # Ra and eta are lists of strings
    if t1 is None:
        t1 = np.zeros((len(eta), len(Ra)))
    if clist is None:
        cmap = plt.cm.get_cmap('jet', 3)
    else:
        cmap = cmap_from_list(clist, cmap_name='regimes')
    fig, ax = plt.subplots(1, 1)

    plot_grid = np.zeros((len(eta), len(Ra)))
    for x, xval in enumerate(Ra):
        for y, yval in enumerate(eta):
            desc = regime_grid[y, x]
            if desc == 'chaotic':
                plot_grid[y, x] = 3
            elif desc == 'transitional':
                plot_grid[y, x] = 2
            elif desc == 'steady':
                plot_grid[y, x] = 1
            else:  # i.e. no convection
                plot_grid[y, x] = np.nan  # label grid with numbers
    m = np.ma.masked_where(np.isnan(plot_grid), plot_grid)
    im = ax.imshow(m, origin='bottom', aspect='equal', interpolation='None',
                   vmin=1, vmax=3, cmap=cmap)

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

    cbar = plt.colorbar(im, ticks=[1.5, 2, 2.5], shrink=0.5)
    cbar.ax.set_yticklabels(['steady', 'transitional', 'chaotic'])

    if overploth:  # do h rms contours, only if already stored rms
        #         for iir in range(1,4): # separate contours for each regime (don't link)
        #             print('plotting contours for regime', iir)
        h_grid = np.zeros((len(eta), len(Ra)))
        h_grid[:] = np.nan
        for x, Raval in enumerate(Ra):
            for y, etaval in enumerate(eta):
                if not np.isnan(plot_grid[y, x]):
                    #                 if plot_grid[y,x] == iir:
                    case = 'Ra' + Raval + '-eta' + etaval + '-wide'
                    sol_df = pickleio(case, suffix='_h', postprocess_functions=[h_at_ts], t1=t1[y, x],
                                      load=load, data_path=data_path, fig_path=fig_path, **kwargs)
                    rms = np.median(sol_df['h_rms'])
                    h_grid[y, x] = rms
        CS = ax.contour(h_grid, nlevels, cmap=cmap_contours)
        ax.clabel(CS, inline=1, fontsize=10)

    if save:
        fig.tight_layout()
        savefig(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt)


def plot_T_params(case, T_params, n=-1,
                  setylabel=True, setxlabel=True, save=True,
                  fig_path=fig_path_bullard, fig=None, ax=None, fname='_T-z', fig_fmt='.png',
                  legend=True, labelsize=16,  **kwargs):
    # take nth row
    if fig is None:
        fig, ax = plt.subplots(figsize=(4, 4))
    try:
        T_params = T_params.iloc[n]
    except IndexError:
        print('No T parameterisation found for solution n =', n)
        return fig, ax

    dT_rh_f = T_params['dT_rh']
    delta_rh_f = T_params['delta_rh']
    D_l_f = T_params['delta_L']
    T_l_f = T_params['T_l']
    T_f = np.array(T_params['T_av'].tolist())
    y_f = np.array(T_params['y'].tolist())

    ax.plot(T_f, y_f, c='k', lw=1)
    ax.axhline(D_l_f, label='$\delta_{L}$', c='xkcd:tangerine', lw=0.5)
    ax.axhline(D_l_f - delta_rh_f, label=r'$\delta_0$', c='xkcd:red orange', lw=0.5)
    ax.text(0, D_l_f - delta_rh_f, r'$\delta_{rh} = $' + '{:04.2f}'.format(delta_rh_f), ha='left', va='top',
            color='xkcd:red orange', fontsize=labelsize - 2)
    ax.plot([T_l_f, T_l_f], [0, D_l_f], ls='--', alpha=0.5, lw=0.5, c='xkcd:tangerine')
    ax.plot([T_l_f + dT_rh_f, T_l_f + dT_rh_f], [0, D_l_f - delta_rh_f], ls='--', alpha=0.5, lw=0.5,
            c='xkcd:red orange')
    if legend:
        ax.legend(frameon=False, fontsize=labelsize - 2)
    if setxlabel:
        ax.set_xlabel('temperature', fontsize=labelsize)
    if setylabel:
        ax.set_ylabel('depth', fontsize=labelsize)
    if save:
        savefig(fig, case+fname, fig_path=fig_path, fig_fmt=fig_fmt)
    return fig, ax


def plot_pdf(case, df=None, keys=None, fig_path=fig_path_bullard, fig=None, ax=None, save=True, settitle=True,
             setxlabel=True, legend=True, labelsize=16, c_list=None, labels=None, fname='h_hist', fig_fmt='.png',
             **kwargs):
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
            ax.text(0.05, 0.05, 'n={:d}/{:d}'.format(len(x), df.index[-1]), ha='left', va='bottom',
                    transform=ax.transAxes, color='b')
        except KeyError:
            print('Key', key, 'not found in', case)

    ax.yaxis.set_ticks([])
    if legend:
        ax.legend(frameon=False, fontsize=labelsize - 2)
    if setxlabel:
        ax.set_xlabel('dynamic topography', fontsize=labelsize)
    if settitle:
        ax.set_title(case, fontsize=labelsize)
    if save:
        savefig(fig, case + fname, fig_path=fig_path, fig_fmt=fig_fmt)
    return fig, ax


# def pd_h_components(case, t1=0, data_path=data_path_bullard, fig_path=fig_path_bullard, sigma=2,
#                     pickleto=None, picklefrom=None, plotTz=False, savefig=False,
#                     settitle=True, setxlabel=True, c='xkcd:pale purple', params_list=None,
#                     legend=True, plotpd=False, labelsize=16, fig=None, ax=None, alpha=None):
#     # probability distribution of h' = f(x) for single case
#     if sigma == 2:
#         qs = [2.5, 50, 97.5]
#     elif sigma == 1:
#         qs = [16, 50, 84]
#
#     if params_list is None:
#         T_params, fig, ax = get_T_params(case, t1=t1, data_path=data_path,
#                                          setxlabel=setxlabel,
#                                          pickleto=pickleto, picklefrom=picklefrom,
#                                          savefig=savefig, fig=fig, ax=ax,
#                                          plotTz=plotTz, fig_path=fig_path)
#     else:
#         T_params = params_list
#
#     if alpha is None:
#         dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False)
#         alpha = dat.parameters['Material model']['Simple model']['Thermal expansion coefficient']
#
#     x_list = alpha * (np.array(T_params['dT_rh']) / np.array(T_params['dT_m'])) * (
#             np.array(T_params['delta_rh']) / np.array(T_params['d_m']))
#     if plotpd and (not plotTz):
#         if ax is None:
#             fig = plt.figure()
#             ax = plt.gca()
#         ax.hist(x_list, color=c, histtype='step')
#         ax.axvline(x=np.median(x_list), color='k', ls='--', label='median')
#         ax.axvline(x=np.mean(x_list), color='k', ls='-', lw=1, label='mean')
#         ax.yaxis.set_ticks([])
#         ax.text(0.95, 0.95, 'n = {:d}'.format(len(x_list)), ha='right', va='top', transform=ax.transAxes)
#         if legend:
#             ax.legend(frameon=False, fontsize=labelsize - 2)
#         if setxlabel:
#             ax.set_xlabel(r'$\Delta T_{rh}/\Delta T_m \; \delta/d_m$', fontsize=labelsize)
#         if settitle:
#             ax.set_title(case, fontsize=labelsize)
#         if savefig:
#             if not os.data_path.exists(fig_path):
#                 os.makedirs(fig_path)
#             fig.savefig(fig_path + case + '_x_hist.png')
#
#     if (not T_params['dT_rh']):  # empty
#         print(case, '- T list is empty')
#         return np.array([np.nan, np.nan, np.nan]), fig, ax
#     return np.percentile(x_list, qs), fig, ax


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
    time, y = read_evol(case, col, dat=dat)
    ax.plot(time, y * yscale, c=c, lw=0.5, label=label)
    ax.set_xlim(0, ax.get_xlim()[1])
    ax.set_xlabel(xlabel, fontsize=labelsize, labelpad=labelpad)
    ax.set_ylabel(ylabel, fontsize=labelsize, labelpad=labelpad)
    if settitle:
        ax.set_title(case, fontsize=labelsize)
    if legend:
        ax.legend(frameon=False, fontsize=labelsize - 2)
    if mark_used:
        # Create a Rectangle patch to mark "transient" times
        rect = patches.Rectangle((ax.get_xlim()[0], ax.get_ylim()[0]), t1,
                                 ax.get_ylim()[1] - ax.get_ylim()[0],
                                 edgecolor='None', facecolor='k', alpha=0.2, zorder=0)
        ax.add_patch(rect)
    if show_sols and sol_df is not None:
        # find steady state sols
        sol_times = np.array(sol_df['time'])
        for t in sol_times:
            ax.axvline(x=t, lw=0.5, ls='--', alpha=0.5, zorder=0)
    if save:
        savefig(fig, case + fname, fig_path=fig_path, fig_fmt=fig_fmt)
    return fig, ax


def plot_top_profile(case, ts, save=True, fig_path=fig_path_bullard, data_path=data_path_bullard, verbose=True,
                     fig_fmt='.png', ):
    x, h = read_topo_stats(case, ts, data_path=data_path)
    # normalize to 0 mean
    h_norm = trapznorm(h)
    fig = plt.figure()
    plt.plot(x, h_norm)
    plt.xlabel('x')
    plt.ylabel('dynamic topography')
    plt.title(case)
    if verbose:
        print('mean:', trapzmean(h_norm))
        print('max:', np.max(h_norm))
        print('min:', np.min(h_norm))
    if save:
        savefig(fig, case + '_h_' + '{:05}'.format(ts), fig_path=fig_path, fig_fmt=fig_fmt)


def cmap_from_list(clist, n_bin=None, cmap_name=''):
    if n_bin is None:
        n_bin = len(clist)
    cm = LinearSegmentedColormap.from_list(cmap_name, clist, N=n_bin)
    return cm
