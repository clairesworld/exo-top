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


def read_topo_stats(case, ts, path='model-output/'):
    df = pd.read_csv(path + 'output-' + case + '/dynamic_topography_surface.' + '{:05}'.format(ts), header=None,
                     names=['x', 'y', 'h'],
                     skiprows=1,
                     index_col=False, delimiter=r"\s+", engine='python')
    return df['x'], df['h']


def pickleio(case, suffix, postprocess_functions, t1=0, load='auto', dat_new=None, datsuffix='_dat',
             data_path=data_path_bullard, fend='.pkl', **kwargs):
    # do pickling strategy
    case_path = data_path + 'output-' + case + '/'
    fname = case + suffix + fend
    run_dict = {}
    dump_flag = False
    reprocess_flag = False
    t1_new = t1

    if os.path.exists(case_path):  # do nothing if case doesn't exist
        os.makedirs(case_path + 'pickle/', exist_ok=True)
        if (load == 'auto' or load is True):
            if os.path.exists(case_path + 'pickle/' + fname):
                # open pickled file
                try:
                    run_dict = pkl.load(open(case_path + 'pickle/' + fname, "rb"))
                except ValueError:  # python 2?
                    run_dict = pkl.load(open(case_path + 'pickle/' + fname, "rb"), protocol=2)
                print('Found', fname)

                if load == 'auto':  # check for additional timesteps
                    time_old = run_dict['time']  # should all be after t1
                    if dat_new is None:
                        try:
                            dat_new = pkl.load(open(case_path + 'pickle/' + datsuffix + fend, "rb"))
                        except FileNotFoundError:
                            dat_new = post.Aspect_Data(directory=case_path, verbose=False,
                                                       read_statistics=True, read_parameters=False)
                            pkl.dump(dat_new, open(case_path + 'pickle/' + datsuffix + fend, "wb"))
                    try:
                        time_new = dat_new.stats_time
                    except AttributeError:
                        dat_new.read_statistics()
                    t1_new = np.argmax(time_new > time_old[-1])
                    if t1_new > 0:  # new timesteps
                        reprocess_flag = True
                        print('Updating', fname)

            elif load == 'auto':  # pkl file not found
                reprocess_flag = True
                print(fname, 'not found, processing...')

        else:  # load is False so automatically calculate shit
            reprocess_flag = True
            dat_new = post.Aspect_Data(directory=case_path, verbose=False,
                                       read_statistics=True, read_parameters=False)

        if reprocess_flag and (t1 != 1):
            try:
                sol_files = dat_new.sol_files
            except AttributeError:
                sol_files = dat_new.read_stats_sol_files()
            run_dict = process_at_solutions(case, postprocess_functions=postprocess_functions, t1=t1_new,
                                            data_path=data_path, dat=dat_new, dict_to_extend=run_dict,
                                            sol_files=sol_files, **kwargs)
            dump_flag = True  # always save if you did something

        if dump_flag:
            pkl.dump(run_dict, open(case_path + 'pickle/' + fname, "wb"))
            pkl.dump(dat_new, open(case_path + 'pickle/' + datsuffix + fend, "wb"))

    return run_dict  # will returning None break?


def get_cases_list(Ra, eta):
    # either Ra or eta is iterable
    if isinstance(Ra, Iterable) and not isinstance(Ra, string_types):
        x_var = Ra
        cases = ['Ra' + r + '-eta' + eta + '-wide' for r in Ra]
    elif (isinstance(eta, Iterable) and not isinstance(eta, string_types)):
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


def read_evol(case, i, path=data_path_bullard, skiprows=None):
    # e.g. return time, column i, nsteps (from statistics file)
    df = pd.read_csv(path + 'output-' + case + '/statistics', header=None,
                     skiprows=skiprows,
                     index_col=False, delimiter=r"\s+", engine='python', comment='#')
    return np.array(df.iloc[:, 1]), np.array(df.iloc[:, i - 1]), len(df.index)


def process_at_solutions(case, postprocess_functions, dat=None, t1=0, data_path=data_path_bullard,
                         dict_to_extend=None, sol_files=None, **kwargs):
    if dict_to_extend is None:
        dict_to_extend = {}
    try:
        dict_to_extend['sol'].extend([])
    except KeyError:  # haven't processed time stamps yet
        dict_to_extend['sol'] = []
        dict_to_extend['time'] = []
        dict_to_extend['timestep'] = []
        dict_to_extend['nsols'] = 0
    if dat is None:
        dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False, read_statistics=True)
    time = dat.stats_time
    if sol_files is None:
        try:
            sol_files = dat.sol_files
        except AttributeError:
            sol_files = dat.read_stats_sol_files()
    i_time = np.argmax(time > t1)  # index of first timestep to process
    if i_time > 0:  # probably don't want to process every timestep
        sols_in_time = sol_files[i_time:]
        n_quasi, n_indices = np.unique(sols_in_time, return_index=True)  # find graphical snapshots within time range
        for ii, n in enumerate(n_quasi):
            n = int(n)
            ts = n_indices[ii]  # timestep at this solution
            for fn in postprocess_functions:
                dict_to_extend = fn(case, n=n, dict_to_append=dict_to_extend, ts=ts, dat=dat, **kwargs)
                print('    calculated', fn, 'for solution', n, '/', int(n_quasi[-1]))

        # always store time tags and AspectData object
        dict_to_extend['sol'].extend(n_quasi)
        dict_to_extend['time'].extend(time[n_indices])
        dict_to_extend['timestep'].extend(n_indices)
        dict_to_extend['nsols'] = dict_to_extend['nsols'] + len(n_quasi)

    else:
        times_at_sols = dat.find_time_at_sol(sol_files=sol_files, return_indices=False)
        print('no quasi-steady state solutions up to', times_at_sols[-1])
    return dict_to_extend


def get_h(case=None, dict_to_extend={}, ts=None, hscale=1, **kwargs):
    p_dict = dict_to_extend
    try:
        h_params_n = {}
        x, h = read_topo_stats(case, ts)
        h_norm = trapznorm(h)
        peak, rms = peak_and_rms(h_norm)
        h_params_n['h_peak'] = peak * hscale
        h_params_n['h_rms'] = rms * hscale

    except FileNotFoundError as e:
        print('    file not found:', e)
        print('    ts =', ts)
        h_params_n['h_peak'] = None
        h_params_n['h_rms'] = None

    for key in h_params_n.keys():
        try:
            p_dict[key].append(h_params_n[key])
        except KeyError:  # key does not exist yet
            p_dict[key] = []
            p_dict[key].append(h_params_n[key])
    return p_dict


def get_T_params(case=None, n=None, dict_to_append={}, dat=None, data_path=data_path_bullard, **kwargs):
    T_params = dict_to_append
    if dat is None:
        try:
            dat = T_params['dat']
        except KeyError:
            dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False, read_statistics=True)
            T_params['dat'] = dat
    n = int(n)
    x, y, z, u, v, _ = dat.read_velocity(n, verbose=False)
    x, y, z, T = dat.read_temperature(n, verbose=False)
    T_params_n = dat.T_components(n, T=T, u=u, v=v, cut=True)
    for key in T_params_n.keys():
        try:
            T_params[key].append(T_params_n[key])
        except KeyError:  # key does not exist yet
            T_params[key] = []
            T_params[key].append(T_params_n[key])
    return T_params, dat


def plot_T_params(case, T_params, n=-1, dat=None,
                  setylabel=True, setxlabel=True, savefig=True,
                  fig_path=fig_path_bullard, fig=None, ax=None,
                  legend=True, labelsize=16, data_path=data_path_bullard, **kwargs):
    if dat is None:
        try:
            dat = T_params['dat']
        except KeyError:
            dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False, read_statistics=True)
            T_params['dat'] = dat
    dT_rh_f = T_params['dT_rh'][n]
    delta_rh_f = T_params['delta_rh'][n]
    D_l_f = T_params['delta_L'][n]
    T_l_f = T_params['T_l'][n]
    T_f = T_params['T_av'][n]
    fig, ax = dat.plot_profile(T_f, fig=fig, ax=ax, xlabel='', ylabel='', c='k', lw=1)
    ax.axhline(D_l_f, label='$z_{lid}$', c='xkcd:tangerine', lw=0.5)
    ax.axhline(D_l_f - delta_rh_f, label=r'$z_\delta$', c='xkcd:red orange', lw=0.5)
    ax.text(0, D_l_f - delta_rh_f, r'$\delta = $' + '{:04.2f}'.format(delta_rh_f), ha='left', va='top',
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
    if savefig:
        fig.savefig(fig_path + case + '-T_z.png', bbox_inches='tight')
    return fig, ax


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





# def get_h_old(case, t1=0, data_path=data_path_bullard, dict_to_extend={}, dat=None,
#           fig_path=fig_path_bullard, hscale=1, **kwargs):
#     flag = False
#     if (picklefrom is not None) and (os.data_path.exists(fig_path + 'data/' + picklefrom)):
#         try:
#             peak_list, rms_list = pkl.load(open(fig_path + 'data/' + picklefrom, "rb"))
#             print('loaded h for case', case)
#         except ValueError:
#             peak_list, rms_list = pkl.load(open(fig_path + 'data/' + picklefrom, "rb"), protocol=2)
#             print('loaded h for case', case)
#         if (not peak_list) or (not rms_list):  # if stored stuff is empty
#             flag = True
#             pickleto = picklefrom
#     else:
#         flag = True
#         if picklefrom is not None:  # save if tried to load but failed
#             pickleto = picklefrom
#             print(picklefrom, 'not found or empty, re-calculating')
#
#     if flag:  # load
#         time, v_rms, nsteps = read_evol(case, i=11, data_path=data_path)
#         # what is the probability distribution of i from t1 to end?
#         i_time = np.argmax(time > t1)
#         rms_list = []
#         peak_list = []
#         t_used = []
#         if t1 != 1:
#             print('building distribution of h for case', case, ', n =', len(range(i_time, len(time))))
#             for ii in range(i_time, len(time)):
#                 try:
#                     x, h = read_topo_stats(case, ii)
#                     h_norm = trapznorm(h)
#                     peak, rms = peak_and_rms(h_norm)
#                     rms_list.append(rms)
#                     peak_list.append(peak)
#                     t_used.append(ii)
#                 except FileNotFoundError as e:
#                     print('file not found:', e)
#
#     try:
#         if pickleto is not None:
#             pkl.dump((peak_list, rms_list), open(fig_path + 'data/' + pickleto, "wb"))
#
#         peak_list_scaled = [a * hscale for a in peak_list]
#         rms_list_scaled = [a * hscale for a in rms_list]
#         return peak_list_scaled, rms_list_scaled
#     except:
#         return [np.nan], [np.nan]


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

def parameter_percentiles(case, p_dict=None, keys=None, t1=0, plot=False, sigma=2, **kwargs):
    # probability distribution for a number of parameters with p_dict containing time evolutions
    if sigma == 2:
        qs = [2.5, 50, 97.5]
    elif sigma == 1:
        qs = [16, 50, 84]
    pdf_dict = {}
    for key in keys:
        try:
            pdf_dict[key] = np.percentile(p_dict[key], qs)
        except KeyError as e:
            raise (key, 'not processed yet for', case)
        except Exception as e:
            raise (e)
            pdf_dict[key] = np.array([np.nan] * len(qs))

    if plot:
        fig, ax = plot_pdf(case, p_dict=p_dict, keys=keys, **kwargs)
        return pdf_dict, fig, ax
    return pdf_dict


def plot_pdf(case, p_dict=None, keys=None, fig_path=fig_path_bullard, fig=None, ax=None, savefig=True, settitle=True,
             setxlabel=True, legend=True, labelsize=16, fend='.png', c_list=None, labels=None, fname='h_hist',
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
        x = p_dict[key]
        c = c_list[ii]
        ax.hist(x, color=c, histtype='step', label=labels[ii])
        extralabel = ''
        if ii == 0:
            extralabel = 'mean'
        ax.axvline(x=np.mean(x), color='k', ls='-', lw=1, label=extralabel)
        if ii == 0:
            extralabel = 'median'
        ax.axvline(x=np.median(x), color='k', ls='--', label=extralabel)

    ax.yaxis.set_ticks([])
    ax.text(0.05, 0.05, 'n={:d}'.format(len(x)), ha='left', va='bottom', transform=ax.transAxes)
    if legend:
        ax.legend(frameon=False, fontsize=labelsize - 2)
    if setxlabel:
        ax.set_xlabel('dynamic topography', fontsize=labelsize)
    if settitle:
        ax.set_title(case, fontsize=labelsize)
    if savefig:
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path + case + fname + fend, bbox_inches='tight')
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


def plot_evol(case, i, fig=None, ax=None, savefig=True, fend='_f.png', mark_used=True, t1=0,
              ylabel='rms velocity', xlabel='time', yscale=1, c='k', settitle=True, setxlabel=True,
              setylabel=True, legend=False, labelsize=16, labelpad=5, label=None, fig_path=fig_path_bullard):
    if not setxlabel:
        xlabel = ''
    if not setylabel:
        ylabel = ''
    if ax is None:
        fig = plt.figure()
        ax = plt.gca()
    time, y, nsteps = read_evol(case, i=i)
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
    if savefig:
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path + case + fend)
    return fig, ax


def case_subplots(cases, labels=None, labelsize=16, labelpad=5, t1=None, save=True, dt_xlim=(0.0, 0.065),
                  fname='cases.png', data_path=data_path_bullard, fig_path=fig_path_bullard,
                  loadh='auto', loadT='auto', includegraphic=False, c_rms='xkcd:forest green', c_peak='xkcd:periwinkle',
                  suptitle='', includepdf=True, includeTz=True, **kwargs):
    # rows are cases, columns are v_rms, q, T(z), hist
    ncases = len(cases)
    ncols = 2
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
    delrow = []
    for ii, case in enumerate(cases):
        icol = 0
        if os.path.exists(data_path + 'output-' + case):
            print('plotting summary for', case)
            if ii == ncases - 1:  # show x label in bottom row only
                setxlabel = True
            else:
                setxlabel = False
            setylabel = True
            legend = False
            if numplotted == 0:
                legend = True

            ax = axes[ii, icol]
            fig, ax = plot_evol(case, i=11, fig=fig, ax=ax, savefig=False, ylabel='rms velocity',
                                c='k', settitle=False, setxlabel=setxlabel, setylabel=setylabel,
                                labelsize=labelsize, labelpad=labelpad, legend=False, mark_used=True, t1=t1[ii])

            ax.text(0.01, 0.95, labels[ii], horizontalalignment='left', verticalalignment='top',
                        transform=ax.transAxes, fontsize=labelsize)

            icol = icol + 1
            ax = axes[ii, icol]
            fig, ax = plot_evol(case, i=20, fig=fig, ax=ax, savefig=False, ylabel='heat flux',
                                c='xkcd:light red', settitle=False, setxlabel=setxlabel, setylabel=setylabel,
                                labelsize=labelsize, labelpad=labelpad, label='top', mark_used=True, t1=t1[ii])
            fig, ax = plot_evol(case, i=19, fig=fig, ax=ax, savefig=False, ylabel='heat flux', yscale=-1,
                                c='xkcd:purple blue', settitle=False, setxlabel=setxlabel, setylabel=setylabel,
                                labelsize=labelsize, labelpad=labelpad, label='bottom', legend=legend, mark_used=False)

            if includeTz:  # final timestep only
                icol = icol + 1
                ax = axes[ii, icol]

                run_dict = pickleio(case, suffix='_parameters', postprocess_functions=[get_T_params], t1=t1[ii],
                                    load=loadT, data_path=data_path, fig_path=fig_path, **kwargs)

                fig, ax = plot_T_params(case, T_params=run_dict, data_path=data_path, savefig=False,
                                        setxlabel=setxlabel,
                                        setylabel=False, legend=False, fig_path=fig_path, fig=fig, ax=ax)
                if legend:
                    ax.legend(frameon=False)
                if setxlabel:
                    ax.set_xlabel('temperature', fontsize=labelsize)
                if setylabel:
                    ax.set_ylabel('depth', fontsize=labelsize)

            if includepdf:
                icol = icol + 1
                ax = axes[ii, icol]
                try:
                    run_dict['h_rms']
                except:
                    run_dict = pickleio(case, suffix='_parameters', postprocess_functions=[get_h], t1=t1[ii],
                                        load=loadh, data_path=data_path, fig_path=fig_path, **kwargs)
                fig, ax = plot_pdf(case, keys=['h_rms', 'h_peak'], p_dict=run_dict, path=data_path,
                                   fig=fig, ax=ax, savefig=False, settitle=False, setxlabel=setxlabel,
                                   legend=legend, labelsize=labelsize, c_list=[c_rms, c_peak])
                ax.set_xlim(dt_xlim[0], dt_xlim[1])  # for fair comparison

            if includegraphic:
                icol = icol + 1
                ax = axes[ii, icol]
                try:
                    img = mpimg.imread(fig_path + 'graphical/' + case + '.png')
                    ax.imshow(img)
                    print('plotting graphical output for', case)
                except FileNotFoundError:
                    print('file graphical/', case, '.png not found')
                    fig.delaxes(ax)

            numplotted += 1
        else:
            print(case, 'not found')
            ax = axes[ii, 0]
            ax.text(0.01, 0.95, labels[ii] + '\n\nno stagnant lid convection',
                    horizontalalignment='left', verticalalignment='top',
                    transform=ax.transAxes, fontsize=labelsize)
            delrow.append(ii)
    #     for mm in delrow: # delete unused subplot rows
    #         for nn in range(icol+1):
    #             try:
    #                 fig.delaxes(axes[mm][nn])
    #             except:
    #                 print('could not delete axis', mm, nn)
    plt.suptitle(suptitle, fontsize=labelsize * 2, y=1.02)
    fig.tight_layout()
    if save:
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path + fname, bbox_inches='tight')
    return fig, axes


def plot_h_vs_Ra(Ra=None, eta=None, t1=None, data_path=data_path_bullard, fig_path=fig_path_bullard,
                 load='auto', showallscatter=False,
                 save=True, fname='h.png', sigma=2,
                 labelsize=16, xlabel='', ylabel='dynamic topography', title='',
                 c_peak='xkcd:forest green', c_rms='xkcd:periwinkle',
                 fit=False, fitRa=None, fitfn='line', cases=None, x_var=None, logx=True, logy=True,
                 fig=None, ax=None, dt_ylim=(3e-3, 7e-2), xlim=None, hscale=1, **kwargs):
    # Ra or eta is list of strings, t1 is a list of numbers the same length
    if cases is None:
        cases, x_var = get_cases_list(Ra, eta)
    if t1 is None:
        t1 = [0] * len(x_var)
    if fitRa is None:
        fitRa = Ra

    quants_h_peak = np.zeros((len(x_var), 3))
    quants_h_rms = np.zeros((len(x_var), 3))
    x = np.zeros((len(x_var), 3))
    peak_all = []
    rms_all = []

    for ii, case in enumerate(cases):
        run_dict = pickleio(case, suffix='_parameters', postprocess_functions=['get_h'], t1=0, load=load, dat_new=None,
                            datsuffix='_dat', data_path=data_path_bullard, fend='.pkl', **kwargs)
        p_dict = get_h(case, t1=t1[ii], data_path=data_path, hscale=hscale, dict_to_extend=run_dict,)
        h_peak = p_dict['h_peak']
        h_rms = p_dict['h_rms']

        quants_h_peak[ii, :], quants_h_rms[ii, :], _, _ = parameter_percentiles(case, p_dict=run_dict, sigma=sigma,
                                                                                keys=['h_peak', 'h_rms'], plot=False)

        x[ii, :] = float(x_var[ii])
        peak_all.append((h_peak, x[ii, :] ))
        rms_all.append((h_rms, x[ii, :] ))

    yerr_peak = [quants_h_peak[:, 1] - quants_h_peak[:, 0], quants_h_peak[:, 2] - quants_h_peak[:, 1]]
    yerr_rms = [quants_h_rms[:, 1] - quants_h_rms[:, 0], quants_h_rms[:, 2] - quants_h_rms[:, 1]]
    xerr = None

    if ax is None:
        fig = plt.figure()
        ax = plt.gca()

    if fit:
        if (len(x_var) > 1):
            fitx = [[a[1]] * len(a[0]) for a in rms_all]
            fith = [a[0] for a in rms_all]
            flatfitx = [item for sublist in fitx for item in sublist]
            flatfith = [item for sublist in fith for item in sublist]
            expon, const = fit_log(flatfitx, flatfith)
            xprime = [a[1] for a in rms_all]
            hprime = const * xprime ** expon
            h3, = ax.plot(xprime, hprime, c=c_rms, ls='--', lw=1, zorder=100,
                          label='{:.2e} Ra^{:.3f}'.format(const, expon))
            ax.legend(
                # handles=[h3], labels=[],
                loc='lower left')  # always show what fit is
        else:
            print('not enough points to fit -- Ra', Ra, 'eta', eta)

    ax.errorbar(x[:, 1], h_peak[:, 1], yerr=yerr_peak, xerr=xerr,
                fmt='^', c=c_peak, alpha=0.9, capsize=5)
    ax.errorbar(x[:, 1], h_rms[:, 1], yerr=yerr_rms, xerr=xerr,
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
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path + fname, bbox_inches='tight')
    return fig, ax


def plot_h_vs_Td(Ra=None, eta=None, t1=None, data_path=data_path_bullard, fig_path=fig_path_bullard,
                 loadpickle=False, dumppickle=False, loadpicklex=False,
                 save=True, fname='h_T.png', plotpd=False, sigma=2, showallscatter=False,
                 labelsize=16, xlabel=r'$\delta_rh \Delta T_{rh}$', ylabel='dynamic topography', title='',
                 c_peak='xkcd:forest green', c_rms='xkcd:periwinkle', legend=True,
                 fit=False, fitRa=None, fitfn='line', cases=None, x_var=None, logx=True, logy=True,
                 fig=None, ax=None, dt_ylim=(3e-3, 7e-2), xlim=None, hscale=1, **kwargs):
    # Ra or eta is list of strings, t1 is a list of numbers the same length
    # instead of plotting vs Ra or eta, plot vs theoretical components of scaling relationship

    if cases is None:
        cases, x_var = get_cases_list(Ra, eta)
    if t1 is None:
        t1 = [0] * len(cases)
    if fitRa is None:
        fitRa = Ra

    h_peak = np.zeros((len(x_var), 3))
    h_rms = np.zeros((len(x_var), 3))
    x = np.zeros((len(x_var), 3))
    peak_all = []
    rms_all = []

    for ii, case in enumerate(cases):
        t1_ii = t1[ii]
        # load h
        picklefile = case + '_pdtop.pkl'
        if loadpickle:
            picklefrom = picklefile
        else:
            picklefrom = None
        if dumppickle:
            pickleto = picklefile
        else:
            pickleto = None
        # assume time-dependent convection (if steady state then shouldn't matter)
        peak_list, rms_list = get_h(case, t1=t1_ii, path=data_path, pickleto=pickleto, picklefrom=picklefrom,
                                    fig_path=fig_path, hscale=hscale)

        try:
            h_peak[ii, :], h_rms[ii, :], _, _ = pdf_h(case, plot=False, t1=t1_ii, path=data_path, fig_path=fig_path,
                                                      peak_list=peak_list, rms_list=rms_list)
        except Exception as e:
            print('aspect_scalings.py:', e, '\n setting h all nan for case', case)
            h_peak[ii, :], h_rms[ii, :] = ([np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan])

            # load T components
        dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False)
        picklefile = case + '_pdx.pkl'
        if loadpicklex:
            picklefrom = picklefile
        else:
            picklefrom = None
        if dumppickle:
            pickleto = picklefile
        else:
            pickleto = None
        # assume time-dependent convection (if steady state then shouldn't matter)
        T_params, _, _ = get_T_params(case, t1=t1_ii, data_path=data_path,
                                      pickleto=pickleto, picklefrom=picklefrom, plotTz=False,
                                      fig_path=fig_path)
        alpha = dat.parameters['Material model']['Simple model']['Thermal expansion coefficient']
        x_list = alpha * (np.array(T_params['dT_rh']) / np.array(T_params['dT_m'])) * (
                np.array(T_params['delta_rh']) / np.array(T_params['d_m']))

        try:
            x[ii, :], _, _ = pdf_h_components(case, t1=t1_ii, data_path=data_path, sigma=sigma,
                                              pickleto=pickleto,
                                              picklefrom=picklefrom, fig_path=fig_path, plotTz=False,
                                              params_list=T_params, alpha=alpha)
        except Exception as e:
            x[ii, :] = (np.nan, np.nan, np.nan)

        #         if fit:
        # extract timesteps in h where you have T snapshot - excluding initial transient

        dat.read_statistics(verbose=False)
        t1_idx = np.argmax(dat.stats_time > t1_ii)  # timestep corresponding to t1
        sol_idx = dat.find_time_at_sol()
        sol_idx = np.array(sol_idx[sol_idx >= t1_idx])  # cut any values of idx below t1_idx
        # account for the stored peak_list and rms_list starting at t1 for indexing 
        if np.shape(sol_idx) != np.shape(x_list):
            print('inconsistent times: sol_idx', np.shape(sol_idx), 'delta', np.shape(T_params['delta_rh']))
        try:
            peak_list = [peak_list[j - t1_idx] for j in sol_idx]
            rms_list = [rms_list[j - t1_idx] for j in sol_idx]
        except IndexError:
            print('sol_idx - t1_idx', [j - t1_idx for j in sol_idx])

        peak_all.append((peak_list, x_list))
        rms_all.append((rms_list, x_list))

    yerr_peak = [h_peak[:, 1] - h_peak[:, 0], h_peak[:, 2] - h_peak[:, 1]]
    yerr_rms = [h_rms[:, 1] - h_rms[:, 0], h_rms[:, 2] - h_rms[:, 1]]
    xerr = [x[:, 1] - x[:, 0], x[:, 2] - x[:, 1]]

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
        if (len(x_var) > 1):  # can only fit if at least 2 data
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

    ax.errorbar(x[:, 1], h_peak[:, 1], yerr=yerr_peak, xerr=xerr,
                fmt='^', c=c_peak, alpha=0.9, capsize=5)
    ax.errorbar(x[:, 1], h_rms[:, 1], yerr=yerr_rms, xerr=xerr,
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
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path + fname, bbox_inches='tight')
    return fig, ax


def fit_log(x, h, plot=True):
    x1 = np.log10(np.array(x))  # this should work for time-series of all x corresponding to h
    h1 = np.log10(np.array(h))
    try:
        slope, intercept, r_value, p_value, std_err = stats.linregress(x1, h1)
    except ValueError:
        print('error | x', np.shape(x1), 'h', np.shape(h1))
    return slope, 10 ** intercept


def fit_h_sigma(x, h, h_err=None, fn='line'):
    def line(x, a, b):
        return a * x + b

    #     def expon(x, C, n):
    #         return C * x**n

    idx = np.nonzero(np.isnan(x) == False)[0]
    x_fit = np.log10(x[idx])
    h_fit = np.log10(h[idx])
    if h_err is not None:
        h_err = np.log10(h_err[idx])
    print('fitting x =', x_fit, 'h =', h_fit)
    if fn == 'line':
        popt, pcov = curve_fit(line, x_fit, h_fit, sigma=h_err)
        print('slope:', popt[0], 'intercept:', popt[1])

    return 10 ** (popt[1] + popt[0] * x)  # h evaluated at x


def subplots_h_vs(Ra_ls, eta_ls, regime_grid, c_regimes, loadpickle=True, dumppickle=False, save=True,
                  sigma=2, t1=None, fit=False, loadpicklex=False, nrows=2, ncols=2, x_components=False,
                  data_path=data_path_bullard, fig_path=fig_path_bullard, fname='h_Ra_all.png',
                  ylim=(6e-3, 7e-2), labelsize=14, xlim=None, xlabel='Ra', ylabel='dynamic topography',
                  logx=True, logy=True, showallscatter=False, xlabelpad=12, ylabelpad=2, hscale=1):
    # subplots for different eta
    #     fig, axes = plt.subplots(2,2, figsize=(7,7))
    #     flaxes = axes.flatten()
    if x_components:
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
                              loadpickle=loadpickle, dumppickle=dumppickle, plotpd=False, xlim=xlim,
                              save=False, fit=True, ylabel='', xlabel='', labelsize=labelsize,
                              fig_path=fig_path, loadpicklex=loadpicklex, logx=logx, logy=logy,
                              showallscatter=showallscatter, hscale=hscale,
                              )
        # trans
        if not (not Ra_trans):
            fig, ax = plot_fn(Ra=Ra_trans, eta=eta, t1=t1[ii, Ra_trans_idx], sigma=sigma, fig=fig, ax=ax,
                              data_path=data_path,
                              c_rms=c_regimes[1], c_peak=c_regimes[1], dt_ylim=ylim, legend=legend,
                              loadpickle=loadpickle, dumppickle=dumppickle, plotpd=False, xlim=xlim,
                              save=False, fit=True, ylabel='', xlabel='', labelsize=labelsize,
                              fig_path=fig_path, loadpicklex=loadpicklex, logx=logx, logy=logy,
                              showallscatter=showallscatter, hscale=hscale,
                              )
            # chaotic
        if not (not Ra_chaos):
            fig, ax = plot_fn(Ra=Ra_chaos, eta=eta, t1=t1[ii, Ra_chaos_idx], sigma=sigma, fig=fig, ax=ax,
                              data_path=data_path,
                              c_rms=c_regimes[2], c_peak=c_regimes[2], dt_ylim=ylim, legend=legend,
                              loadpickle=loadpickle, dumppickle=dumppickle, plotpd=False, xlim=xlim,
                              save=False, fit=True, xlabel='', ylabel='', labelsize=labelsize,
                              fig_path=fig_path, loadpicklex=loadpicklex, logx=logx, logy=logy,
                              showallscatter=showallscatter, hscale=hscale,
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
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path + fname,  # bbox_inches='tight',
                    bbox_extra_artists=(outer_legend,))
    return fig, axes


def Ra_scaling(Ra_data=None, y_data=None, t1=None, path=data_path_bullard, fig_path=fig_path_bullard,
               save=True, fname='claire.png', sigma=2, showallscatter=False,
               labelsize=16, ylabel='', xlabel='Ra', title='',
               c_scatter='xkcd:forest green', legend=True,
               fit=False, cases=None, x_var=None, logx=True, logy=True,
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
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path + fname, bbox_inches='tight')
    return fig, ax


def plot_bl_Nu_scaling(Ra=None, eta=None, t1=None, data_path=data_path_bullard, fig_path=fig_path_bullard,
                       loadpicklex=False, save=True, fname='bl-Nu.png', sigma=2, showallscatter=False,
                       labelsize=16, ylabel=[r'$\delta$', 'Nu'], xlabel='Ra', title='',
                       c_scatter='xkcd:forest green', legend=True, cmap='magma', compare_pub=None, compare_label=None,
                       fitdelta=False, fitNu=False, x_var=None, logx=True, logy=True, vmin=4, vmax=9,
                       fig=None, axes=None, ylim=None, xlim=None, **kwargs):
    # Ra or eta is list of strings, t1 is a list of numbers the same length
    # instead of plotting vs Ra or eta, plot vs theoretical components of scaling relationship
    if sigma == 2:
        qs = [2.5, 50, 97.5]
    elif sigma == 1:
        qs = [16, 50, 84]
    if fig is None:
        fig, axes = plt.subplots(2, 1)
    logeta_fl = [np.log10(float(a)) for a in eta]
    c_list = colorize(logeta_fl, cmap=cmap)[0]

    for jj, eta_str in enumerate(eta):
        cases, Ra_var = get_cases_list(Ra, eta_str)
        if t1 is None:
            t1 = [0] * len(cases)
        else:
            t1_eta = t1[jj]
        c_scatter = c_list[jj]

        Ra_plot = []
        Nu_plot = []
        delta_0_plot = []

        for ii, case in enumerate(cases):
            t1_ii = t1_eta[ii]

            if (t1_ii != 1) and (os.path.exists(data_path + 'output-' + case)):
                Ra_plot.append(float(Ra_var[ii]))

                # load T components  
                dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False)
                picklefile = case + '_pdx.pkl'
                if loadpicklex:
                    picklefrom = picklefile
                else:
                    picklefrom = None

                # assume time-dependent convection (if steady state then shouldn't matter)
                T_params, _, _ = get_T_params(case, t1=t1_ii, data_path=data_path,
                                              pickleto=None, picklefrom=picklefrom, plotTz=False,
                                              fig_path=fig_path)

                # extract Nu
                dat.read_statistics(verbose=False)
                Nu = dat.Nu(k=1)
                t1_idx = np.argmax(dat.stats_time > t1_ii)  # timestep corresponding to t1
                Nu = Nu[t1_idx:]
                Nu_plot.append(np.median(Nu))
                delta_0_plot.append(np.median(T_params['delta_0']))  # compare delta_0 and MS95

            if compare_pub is not None:
                cmplabel = compare_label
                if (jj > 0) and (ii > 0):
                    cmplabel = None
                Ra_i, delta_cmp, Nu_cmp = compare_pub(Ra=Ra_plot[ii], d_eta=float(eta_str), T_params=T_params)
                axes[0].plot(Ra_i, delta_cmp, '^', alpha=0.7, c=c_scatter, label=cmplabel)
                axes[1].plot(Ra_i, Nu_cmp, '^', alpha=0.7, c=c_scatter, label=cmplabel)

        fig, axes[0] = Ra_scaling(Ra_data=Ra_plot, y_data=delta_0_plot, t1=t1,
                                  path=data_path, fig_path=fig_path,
                                  save=False, sigma=sigma, showallscatter=False,
                                  labelsize=labelsize, ylabel=ylabel[0], xlabel='Ra',
                                  c_scatter=c_scatter, fit=fitdelta, logx=logx, logy=logy,
                                  fig=fig, ax=axes[0], ylim=ylim, xlim=xlim, vmin=vmin, vmax=vmax)
        fig, axes[1] = Ra_scaling(Ra_data=Ra_plot, y_data=Nu_plot, t1=t1,
                                  path=data_path, fig_path=fig_path,
                                  save=False, sigma=sigma, showallscatter=False,
                                  labelsize=labelsize, ylabel=ylabel[1], xlabel='Ra',
                                  c_scatter=c_scatter, fit=fitNu, logx=logx, logy=logy,
                                  fig=fig, ax=axes[1], ylim=ylim, xlim=xlim, vmin=vmin, vmax=vmax)

    scat = axes[1].scatter(logeta_fl, logeta_fl, visible=False, c=np.array(logeta_fl), cmap=cmap,
                           vmin=vmin, vmax=vmax)  # dummy
    cbar = fig.colorbar(scat, ax=[axes[0], axes[1]])
    cbar.set_label(r'log($\Delta \eta$)', fontsize=labelsize, rotation=270, labelpad=18)

    if save:
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path + fname, bbox_inches='tight')
    return fig, axes


def solomatov95(Ra=None, d_eta=None, T_params=None, case=None,
                path=data_path_bullard, fig_path=fig_path_bullard, picklefrom=None):
    if T_params is None:
        T_params, _, _ = get_T_params(case=case, data_path=path, fig_path=fig_path, picklefrom=picklefrom)
    T0 = 1
    dT = 1
    T_i = np.median(T_params['T_i'])
    gamma = np.log(d_eta)  # gamma for this delta eta
    p = gamma * dT
    eta_0 = np.exp(-gamma * T0)
    eta_i = np.exp(-gamma * T_i)
    Ra_i = np.array(Ra) * eta_0 / eta_i
    delta_0 = 1.85 * p ** 1.3185 * Ra_i ** -0.3185
    delta_1 = p ** -1 * delta_0
    Nu = (delta_0 + delta_1) ** -1
    return Ra_i, delta_0, Nu


def moresi95(Ra=None, d_eta=None, T_params=None, case=None,
             data_path=data_path_bullard, fig_path=fig_path_bullard, picklefrom=None):
    if T_params is None:
        T_params, _, _ = get_T_params(case=case, data_path=data_path, fig_path=fig_path, picklefrom=picklefrom)
    T0 = 1
    dT = 1
    T_i = np.median(T_params['T_i'])
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


def plot_convection_regimes(Ra, eta, regime_grid, path=data_path_bullard, fig_path=fig_path_bullard, loadpickle=False,
                            dumppickle=False, save=True, fname='regimes.png', labelsize=16, sigma=2,
                            overploth=False, nlevels=10, clist=None, cmap_contours='spring', **kwargs):
    # Ra and eta are lists of strings
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
                    peak, rms, _, _ = pdf_h(case, path=path, fig_path=fig_path,
                                            sigma=sigma, plot=False, picklefrom=case + '_pdtop.pkl')
                    h_grid[y, x] = rms[1]
        CS = ax.contour(h_grid, nlevels, cmap=cmap_contours)
        ax.clabel(CS, inline=1, fontsize=10)

    fig.tight_layout()
    if save:
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path + fname, bbox_inches='tight')


def plot_top_profile(case, savefig=True, fig_path=fig_path_bullard, path=data_path_bullard, verbose=True):
    time, y, nsteps = read_evol(case, i=2, path=path)
    snap = nsteps - 2  # honestly not sure why it's -2 but seems to work haha
    x, h = read_topo_stats(case, snap)
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
    if savefig:
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path + case + '_h_' + '{:05}'.format(snap) + '.png')


def plot_pd_steadystate(case, i, t1, path=data_path_bullard):
    time, y, nsteps = read_evol(case, i, path=path)
    x, h = read_topo_stats(case, nsteps - 2)
    h_norm = trapznorm(h)
    peak, rms = peak_and_rms(h_norm)
    # what is the probability distribution of i from t1 to end?
    i_time = np.nonzero(time > t1)
    print('transience ends at timestep', i_time[0][0])
    fig = plt.figure()
    plt.gca().hist(y[i_time])


def cmap_from_list(clist, n_bin=None, cmap_name=''):
    if n_bin is None:
        n_bin = len(clist)
    cm = LinearSegmentedColormap.from_list(cmap_name, clist, N=n_bin)
    return cm
