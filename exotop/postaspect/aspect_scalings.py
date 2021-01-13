""" ASPECT runs: functions for bulk postprocessing, storing/loading data, extracting scaling relationships, and plotting """

import numpy as np
import pandas as pd
import pickle as pkl
from scipy.optimize import curve_fit
from scipy import stats
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.image as mpimg
from matplotlib.colors import LogNorm, Normalize
import matplotlib.lines as mlines  # noqa: E402
from sklearn import linear_model
# import statsmodels.api as sm
import sys

sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.postaspect import aspect_postprocessing2 as post  # noqa: E402
from exotop.useful_and_bespoke import colorize, iterable_not_string, cmap_from_list, printe, not_iterable, \
    colourbar, not_string  # noqa: E402

data_path_bullard = '/raid1/cmg76/aspect/model-output/'
fig_path_bullard = '/raid1/cmg76/aspect/figs/'
highlight_colour = 'xkcd:coral'


def plot_save(fig, fname, fig_path=fig_path_bullard, fig_fmt='.png', bbox_inches='tight', tight_layout=True, **kwargs):
    path = fig_path + fname + fig_fmt
    directory = os.path.dirname(path)
    os.makedirs(directory, exist_ok=True)
    if tight_layout:
        fig.tight_layout()
    fig.savefig(path, bbox_inches=bbox_inches, **kwargs)
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  saved to ', path, '!\n')


def read_topo_stats(case, ts, data_path=data_path_bullard):
    df = pd.read_csv(data_path + 'output-' + case + '/dynamic_topography_surface.' + '{:05}'.format(ts), header=None,
                     names=['x', 'y', 'h'], skiprows=1, index_col=False, delimiter=r"\s+", engine='python')
    return df['x'], df['h']


def reshape_one_input(A, proper, default):
    if not_iterable(A) and not_string(A):
        B = np.empty_like(proper, dtype=object)
        if A is None:
            B[:] = default
        else:
            B[:] = A
    elif not np.shape(A) == np.shape(proper):
        B = np.reshape(A, np.shape(proper))
    else:
        B = A
    return B


def reshape_inputs(Ra_ls, eta_ls, grids, defaults=None):
    if defaults is None:
        defaults = (0, 'auto', '', '')

    n_Ra = np.size(Ra_ls)
    n_eta = np.size(eta_ls)
    proper = np.zeros((n_eta, n_Ra))  # model for arrays with the right shapeRa

    grids_new = []
    for grid, default in zip(grids, defaults):
        grids_new.append(reshape_one_input(grid, proper, default))

    # make sure Ra and eta are iterable the right way
    if not not_string(Ra_ls):  # if a string
        Ra_ls = [Ra_ls]
    if not not_string(eta_ls):
        eta_ls = [eta_ls]

    return Ra_ls, eta_ls, grids_new


def pickleio(case, suffix, postprocess_functions, t1=0, load='auto', dat_new=None, at_sol=True,
             data_path=data_path_bullard, fend='.pkl', time_average=False, **kwargs):
    # do pickling strategy. only saving processed data for runs in quasi-steady state (past their t1 value)
    case_path = data_path + 'output-' + case + '/'
    fname = case + suffix + fend
    df = pd.DataFrame()
    dump_flag = False
    dump_flag2 = False
    reprocess_flag = False
    t1_new = t1

    if t1 < 1:
        if os.path.exists(case_path):  # do nothing if case doesn't exist
            os.makedirs(case_path + 'pickle/', exist_ok=True)
            if (load == 'auto' or load) and os.path.exists(case_path + 'pickle/' + fname):
                # open pickled file
                try:
                    df = pkl.load(open(case_path + 'pickle/' + fname, "rb"))
                except ValueError:  # python 2?
                    df = pkl.load(open(case_path + 'pickle/' + fname, "rb"), protocol=2)
                print('    Found', fname)

                if load == 'auto':  # check for additional timesteps
                    if dat_new is None:
                        dat_new = post.Aspect_Data(directory=case_path, verbose=False,
                                                   read_statistics=True, read_parameters=False)
                    try:
                        if at_sol:
                            print('      Checking for new solutions...')
                            sol_f_old = df.sol.iat[-1]
                            sol_new = dat_new.read_stats_sol_files()
                            sol1_new = sol_new[np.argmax(sol_new > sol_f_old)]  # first solution after latest saved
                            t1_new = dat_new.find_time_at_sol(n=sol1_new, sol_files=sol_new, return_indices=False)
                        else:
                            print('      Checking for new timesteps...')
                            time_f_old = df.time.iat[-1]
                            time_new = dat_new.stats_time
                            t1_new = time_new[
                                np.argmax(time_new > time_f_old)]  # first time after latest saved time
                    except AttributeError as e:  # i.e. sol not found in df (because it's empty?)
                        reprocess_flag = True
                    if t1_new > 0:  # new timesteps
                        reprocess_flag = True
                        print('      Updating', fname, 'from t = {:4f}'.format(t1_new))

            elif load:
                print('    File', fname, 'not found')
            else:
                if not load:  # load is False so automatically calculate shit
                    print('    Declined to load', fname, ', processing afresh...')
                elif load == 'auto':
                    print('    File', fname, 'not found, processing...')
                else:
                    raise Exception('load value not understood:', load, type(load))
                reprocess_flag = True
                dat_new = post.Aspect_Data(directory=case_path, verbose=False,
                                           read_statistics=True, read_parameters=False)
                if time_average:
                    # some necessary things u don't have yet
                    i_time = np.argmax(dat_new.stats_time >= t1)
                    sol_new = dat_new.read_stats_sol_files()
                    sol1_new = sol_new[i_time]  # n0 in process_at_average

            if reprocess_flag:
                if not hasattr(dat_new, 'stats_time'):
                    dat_new.read_statistics()
                if at_sol:
                    if not hasattr(dat_new, 'sol_files'):
                        dat_new.read_stats_sol_files()
                    sol_new = dat_new.sol_files
                    if time_average is not 'only':
                        df = process_at_solutions(case, postprocess_functions=postprocess_functions, dat=dat_new,
                                                  t1=np.maximum(t1, t1_new),  # whichever comes later in time
                                                  data_path=data_path, sol_files=sol_new, df_to_extend=df, **kwargs)
                        dump_flag = True  # always save if you did something
                    # if time_average or time_average == 'only':
                    #     try:
                    #         df2 = pkl.load(open(case_path + 'pickle/' + fname + '_average', "rb"))
                    #     except FileNotFoundError:
                    #         df2 = None
                    #     df2 = process_at_average(case, postprocess_functions=postprocess_functions,
                    #                              n0=sol1_new, nf=sol_new[-1], dat=dat_new,
                    #                              data_path=data_path, df_old=df2, **kwargs)
                    #     dump_flag2 = True
                else:
                    df = process_steadystate(case, postprocess_functions=postprocess_functions, dat=dat_new,
                                             t1=np.maximum(t1, t1_new),
                                             data_path=data_path, df_to_extend=df, **kwargs)
                    dump_flag = True  # always save if you did something


            if dump_flag:
                pkl.dump(df, open(case_path + 'pickle/' + fname, "wb"))
            # if dump_flag2:
            #     pkl.dump(df2, open(case_path + 'pickle/' + fname + '_average', "wb"))
    else:
        print('Skipping case', case, 'for t1 <= 1')
    return df


def pickle_and_postprocess(cases, suffix, postprocess_functions, t1=0, **kwargs):
    for ii, case in enumerate(cases):
        pickleio(case, suffix, postprocess_functions, t1=t1, **kwargs)


def pickle_drop_duplicate_row(case, suffix, which='sol', fend='.pkl', data_path=data_path_bullard):
    # remove duplicate rows or nan (e.g. for solution or timestep) - for when you fucked up storing
    case_path = data_path + 'output-' + case + '/'
    fname = case + suffix + fend

    if os.path.exists(case_path + 'pickle/' + fname):
        df = pkl.load(open(case_path + 'pickle/' + fname, 'rb'))
        try:
            series = df[which]
            unique = ~series.duplicated()  # boolean array of duplicates
            df_new = df[unique]
            if df_new.equals(df):
                print('pickle_remove_duplicate(', case, suffix, which, '): No duplicate rows found')
            pkl.dump(df_new, open(case_path + 'pickle/' + fname, 'wb'))
        except KeyError:
            print('pickle_remove_duplicate(', case, suffix, which, '):', fname, 'does not contain column', which)
    else:
        print('pickle_remove_duplicate(', case, suffix, which, '): File', fname, 'not found')


def pickle_drop(case, suffix, keys=None, index=None, fend='.pkl', errors='ignore', data_path=data_path_bullard,
                **kwargs):
    case_path = data_path + 'output-' + case + '/'
    fname = case + suffix + fend
    df = pkl.load(open(case_path + 'pickle/' + fname, "rb"))  # open pickled file
    if keys is not None:  # drop columns
        bad = []
        for key in keys:
            if key not in df.columns:
                bad.append(key)
        if not not bad:
            print('pickle_drop(', case, '): Keys', bad, ' not found to drop')
        df2 = df.drop(labels=keys, axis=1, errors=errors)
    elif index is not None:  # drop rows
        df2 = df.drop(labels=index, axis=0, errors=errors)
    else:
        raise Exception('pickle_drop(): Must provide keys or index to drop')
    if not df2.equals(df):
        pkl.dump(df2, open(case_path + 'pickle/' + fname, "wb"))


def pickle_concat(case, keys=None, suffixes=None, new_suffix=None, fend='.pkl', data_path=data_path_bullard):
    if new_suffix is None:
        new_suffix = '_'
    if keys is None:
        copy_all_keys = True
    else:
        copy_all_keys = False
    case_path = data_path + 'output-' + case + '/'
    dfs = []
    for suffix in suffixes:
        fname = case + suffix + fend
        if os.path.exists(case_path + 'pickle/' + fname):
            df_new = pd.DataFrame()
            df_loaded = pkl.load(open(case_path + 'pickle/' + fname, "rb"))  # open pickled file
            bad = []
            if copy_all_keys:
                keys = df_loaded.columns.values
            for key in keys:
                if key in df_loaded.columns:
                    df_new[key] = df_loaded[key]  # add this column to new df
                    # print('Copied column', key, 'from', fname)
                else:
                    bad.append(key)
            dfs.append(df_new)
            if not not bad:
                print('pickle_concat(): File', fname, 'does not contain', bad)
        else:
            print('pickle_concat(): File', fname, 'does not exist')
    try:
        df_new = pd.concat(dfs, axis=1)  # concatenate along col axis
    except ValueError as e:
        print([d['sol'] for d in dfs])
        raise e
    is_dup_col = df_new.columns.duplicated()
    if is_dup_col.any():  # check for duplicate cols
        for i, col in enumerate(df_new.columns[df_new.columns.duplicated(keep=False)]):
            # raise error if duplicated columns do not match
            df_dups = df_new.loc[:, col]
            if ~df_dups.eq(df_dups.iloc[:, 0], axis=0).all(1).any():  # if any rows do not all match first col
                raise Exception(
                    'Attempting to delete duplicate columns which do not match in ' + case_path + ', ' + suffixes)
        df_new = df_new.loc[:, ~is_dup_col]  # remove duplicate cols if any
    pkl.dump(df_new, open(case_path + 'pickle/' + case + new_suffix + fend, "wb"))


def print_solution_data(case, suffix='_T', keys=None, data_path=data_path_bullard, fend='.pkl'):
    case_path = data_path + 'output-' + case + '/'
    fname = case + suffix + fend
    df_print = pd.DataFrame()
    if os.path.exists(case_path + 'pickle/' + fname):
        badkeys = []
        df = pkl.load(open(case_path + 'pickle/' + fname, "rb"))  # open pickled file
        if keys is None:
            keys = df.columns.values
        for key in keys:
            if key in df.columns:
                df_print[key] = df[key]
            else:
                badkeys.append([key])
        print(df_print)
        if not not badkeys:
            print('Keys not found:', badkeys)
        print('File name:', fname, '| length:', len(df_print.index))
        print('Columns:', df_print.columns.values)
    else:
        print('File', fname, 'does not exist')
    return df_print


def get_cases_list(Ra, eta, end=None):
    # either Ra or eta is iterable
    if iterable_not_string(Ra) and iterable_not_string(eta):
        # end is 2d grid
        x_var = None
        cases = []
        for e, ennn in zip(eta, end):
            cases.extend(['Ra' + r + '-eta' + e + en for r, en in zip(Ra, ennn)])
    elif iterable_not_string(Ra):
        x_var = Ra
        if end is None:
            end = [''] * len(x_var)
        cases = ['Ra' + r + '-eta' + eta + en for r, en in zip(Ra, end)]
    elif iterable_not_string(eta):
        x_var = eta
        if end is None:
            end = [''] * len(x_var)
        cases = ['Ra' + Ra + '-eta' + eta + en for e, en in zip(eta, end)]
    else:
        raise Exception('At least one of Ra or eta must be iterable')
    return cases, x_var


def trapznorm(A):
    mean = np.trapz(A) / (len(A) - 1)
    #     print('original mean:',mean)
    return A - mean


def trapzmean(A):
    return np.trapz(A) / (len(A) - 1)


def peak_and_rms(h):
    return np.max(h), np.sqrt(trapzmean(h ** 2))


def time_averaged_profile(case, n0, nf, which='temperature', dat=None, data_path=data_path_bullard, verbose=False, **kwargs):
    " get time-averaged profile for T, u, v etc. "
    profs_time = []
    if dat is None:
        dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', read_statistics=False, verbose=False)
    for n in range(n0, nf + 1):
        print('Loading', which, 'profile at n =', n)
        if which == 'temperature':
            x, y, _, T = dat.read_temperature(n, verbose=verbose)
            prof = post.horizontal_mean(T, x)
        elif which == 'velocity':
            x, y, _, _, _, _, uv_mag = dat.read_velocity(n, verbose=verbose)
            prof = post.horizontal_mean(uv_mag, x)
        profs_time.append(prof)
    profs_time = np.array(profs_time)
    profs_mean = np.mean(profs_time, axis=0)
    return profs_mean, y


def time_averaged_profile_from_df(df, col):
    nsols = len(df)
    srs = df[col].to_numpy()
    y = df.y.to_numpy()[0]
    profs = np.zeros((nsols, len(srs[0])))
    for ii in range(nsols):
        profs[ii, :] = srs[ii]
    av = np.mean(profs, axis=0)
    return av, y


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
                         df_to_extend=None, sol_files=None, postprocess_kwargs={}, **kwargs):
    if df_to_extend is None:
        df_to_extend = pd.DataFrame()
    if dat is None:
        dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False, read_statistics=True,
                               read_parameters=False)

    time = dat.stats_time
    i_time = np.argmax(time >= t1)  # index of first timestep to process

    if i_time > 0:
        if sol_files is None:
            try:
                sol_files = dat.sol_files
            except AttributeError:
                sol_files = dat.read_stats_sol_files()
        sols_in_time = sol_files[i_time:]
        n_quasi, n_indices = np.unique(sols_in_time, return_index=True)  # find graphical snapshots within time range
        n_ts = n_indices + i_time
        if not isinstance(postprocess_functions, list):
            postprocess_functions = [postprocess_functions]
        # print('df_to_extend before reset', df_to_extend)
        df_to_extend = df_to_extend.reset_index()
        # print('df_to_extend after reset', df_to_extend)
        for ii, n in enumerate(n_quasi):

            ts = n_ts[ii]  # timestep at this solution
            # print(ts, 'timestep at solution', n)
            for fn in postprocess_functions:
                new_params_dict = fn(case, n=n, ts=ts, dat=dat, **postprocess_kwargs, **kwargs)
                new_params_dict['sol'] = int(n)
                new_params_dict['time'] = time[ts]
                new_params_dict['ts'] = int(ts)
                # try:  # need to do this bullshit because adding array to single row breaks df init
                #     new_params = pd.DataFrame(new_params_dict, index=[ts])
                # except ValueError as e:
                #     print('adding to df as list :((((')
                #     print(e)
                #     new_params = pd.DataFrame({k: [v] for k, v in new_params_dict.items()}, index=[ts])
                df_to_extend = df_to_extend.append(new_params_dict, ignore_index=True)
                # print('appending\n', new_params_dict)
                # df_to_extend = pd.concat([df_to_extend, new_params])  # concat row axis (may cause duplicate index)
                print('        Processed', fn, 'for solution', n, '/', int(n_quasi[-1]), '@ ts', ts)
        try:
            df_to_extend.ts = df_to_extend.ts.astype(int)
            df_to_extend.sol = df_to_extend.sol.astype(int)
        except ValueError as e:
            print(df_to_extend.ts)
            print(df_to_extend.sol)
            raise e
        df_to_extend = df_to_extend.set_index('ts')
    else:
        # new_params = pd.DataFrame({'sol':[None], 'time':[None]}, index=[0])
        # df_to_extend = df_to_extend.combine_first(new_params)
        if t1 < 1:
            print('    No timesteps after t = {:.2f} (tf = {:.2f}, t1 = {:.2f})'.format(time[i_time], time[-1], t1))
        else:
            print('    Skipping case with t1 > 1')
    return df_to_extend

#
# def process_at_average(case, postprocess_functions, n0, nf, data_path=data_path_bullard, dat=None,
#                        df_old=None, postprocess_kwargs={}, **kwargs):
#     print('nf', nf)
#     print('n0', n0)
#     n_sols_new = nf - n0
#     if dat is None:
#         dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False, read_statistics=False,
#                                read_parameters=False)
#
#     # get time averages for new solutions
#     T_av, y = time_averaged_profile(case, n0, nf, which='temperature', dat=dat, data_path=data_path, verbose=False, **kwargs)
#     uv_mag_av, y = time_averaged_profile(case, n0, nf, which='velocity', dat=dat, data_path=data_path, verbose=False, **kwargs)
#
#     if df_old is None:
#         n_sols_old = 0
#         T_av_old = T_av
#         uv_mag_av_old = uv_mag_av
#
#     else:
#         n_sols_old = df_old['n_sols']
#         T_av_old = df_old['T_av']
#         uv_mag_av_old = df_old['uv_mag_av']
#
#     # find overall average including older timesteps
#     T_av = np.average(np.vstack((T_av_old, T_av)), axis=0, weights=(n_sols_old, n_sols_new))
#     uv_mag_av = np.average(np.vstack((uv_mag_av_old, uv_mag_av)), axis=0, weights=(n_sols_old, n_sols_new))
#
#     dfs = []
#     for fn in postprocess_functions:
#         try:
#             new_params_dict = fn(case, n=None, T_av=T_av, uv_mag_av=uv_mag_av, dat=dat, **postprocess_kwargs, **kwargs)
#         except Exception:
#             print('         Could not process time-average for', fn, 'getting stats average instead (TODO)')
#         new_params_dict['n_sols'] = n_sols_old + nf - n0
#         dfs.append(pd.from_dict(new_params_dict))
#         print('        Processed', fn, 'for time-average, solutions', n0, 'to', nf)
#     df = pd.concat(dfs, axis=1)
#     df = df.loc[:, ~df.columns.duplicated()]
#
#     print(df)
#     return df


def process_steadystate(case, postprocess_functions, dat=None, t1=0, data_path=data_path_bullard,
                        df_to_extend=None, postprocess_kwargs={}, **kwargs):
    if df_to_extend is None:
        df_to_extend = pd.DataFrame()

    if dat is None:
        dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False, read_statistics=True,
                               read_parameters=False)
    time = dat.stats_time
    i_time = np.argmax(time > t1)  # index of first timestep to process (because argmax returns 1st occurrence of 1)

    if i_time > 0:
        print('        Processing', postprocess_functions, 'from timestep', i_time, 'to', len(time))
        if not isinstance(postprocess_functions, list):
            postprocess_functions = [postprocess_functions]
        for ii in range(i_time, len(time)):
            ts = ii  # timestep at this solution
            for fn in postprocess_functions:
                new_params_dict = fn(case, n=None, ts=ts, dat=dat, **postprocess_kwargs, **kwargs)
                new_params_dict['time'] = time[ts]
                new_params = pd.DataFrame(new_params_dict, index=[ts])
                df_to_extend = pd.concat([df_to_extend, new_params])

    else:
        # new_params = pd.DataFrame({'sol':[None], 'time':[None]}, index=[0])
        # df_to_extend = df_to_extend.combine_first(new_params)
        if t1 < 1:
            print('    No timesteps after t = {:.2f} (tf = {:.2f})'.format(time[i_time], time[-1]))
        else:
            print('    Skipping case with t1 > 1')
    return df_to_extend


def h_at_ts(case, ts=None, **kwargs):
    h_params_n = {}
    try:
        x, h = read_topo_stats(case, ts)
        h_norm = trapznorm(h)
        peak, rms = peak_and_rms(h_norm)
        h_params_n['h_peak'] = peak
        h_params_n['h_rms'] = rms

    except FileNotFoundError:
        print('    No dynamic topography found at ts =', ts)
        # h_params_n['h_peak'] = np.nan
        # h_params_n['h_rms'] = np.nan

    # for key in h_params_n.keys():
    #     h_params_n[key] = [h_params_n[key]]
    return h_params_n


def h_timeaverage(case, ts0, tsf=1e50):
    h_params = {}
    h_all = []
    flag = True
    print('  Reading', tsf - ts0, 'files...')
    while flag and ts0 <= tsf:
        try:
            x, h = read_topo_stats(case, ts0)
            h_norm = trapznorm(h)
            h_all.append(h_norm)
            ts0 = ts0 + 1
        except FileNotFoundError:
            print('    No dynamic topography found at ts =', ts0)
            flag = False
    if not h_all:
        raise Exception(case+'ts0 does not catch topography')
    h_all = np.vstack(h_all)
    peak, rms = peak_and_rms(np.mean(h_all, axis=0))
    h_params['h_peak'] = [peak]
    h_params['h_rms'] = [rms]
    h_params['n'] = [np.shape(h_all)[0]]
    return pd.DataFrame.from_dict(h_params)


def Nu_at_ts(case, ts=None, dat=None, data_path=data_path_bullard, **kwargs):
    if dat is None:
        dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False, read_statistics=True,
                               read_parameters=True)
    N = dat.nusselt(k=1, **kwargs)
    return {'Nu': N[ts]}


def T_components_of_h(case, df=None, dat=None, psuffix='_T', data_path=data_path_bullard, update=False,
                      fend='.pkl', alpha_m=None, postprocess_kwargs={}, **kwargs):
    # calculate T components in h heuristic for all processed solutions
    if 'alpha_m' in postprocess_kwargs:
        if (alpha_m is not None) and (alpha_m != postprocess_kwargs['alpha_m']):
            raise Exception('Error: competing alpha_m values passed to T_componenents_of_h. What the heck is going on!')
        alpha_m = postprocess_kwargs['alpha_m']
    if df is None:
        df = pickleio(case, suffix=psuffix, postprocess_functions=[T_parameters_at_sol],
                      dat_new=dat, data_path=data_path, at_sol=True, fend=fend, **kwargs)

    try:
        h_components = alpha_m * (np.array(df['dT_rh']) / np.array(df['dT_m'])) * (
                np.array(df['delta_rh']) / np.array(df['d_m']))
    except KeyError:
        print(df)

    if update:
        df['h_components'] = h_components
        pkl.dump(df, open(data_path + 'output-' + case + '/pickle/' + case + psuffix + fend, 'wb'))

    return h_components


def T_parameters_at_sol(case, n, dat=None, T_av=None, uv_mag_av=None, data_path=data_path_bullard, **kwargs):
    if dat is None:
        dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False,
                               read_statistics=True, read_parameters=False)
    if uv_mag_av is None:
        x, y, z, u, v, _, uv_mag = dat.read_velocity(n, verbose=False)
        uv_mag_av = post.horizontal_mean(uv_mag, x)
    if T_av is None:
        x, y, z, T = dat.read_temperature(n, verbose=False)
        T_av = post.horizontal_mean(T, x)
    d_n = dat.T_components(n, T_av=T_av, uv_mag_av=uv_mag_av, data_path=data_path, **kwargs)  # dict of components just at solution n
    d_n['h_components'] = T_components_of_h(case, df=d_n, dat=dat, data_path=data_path, **kwargs)

    return d_n


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

def parameter_percentiles(case=None, df=None, keys=None, plot=False, sigma=2, **kwargs):
    # probability distribution for a number of parameters with df containing time evolutions
    # df can also be a dict...
    if sigma == 2:
        qs = [2.5, 50, 97.5]
    elif sigma == 1:
        qs = [16, 50, 84]
    else:
        print('Unrecognized sigma value')
    qdict = {}
    for key in keys:
        vals = df[key]  #.values
        try:
            qdict[key] = np.percentile(vals, qs)
        except TypeError:
            vals = [np.array(a).item() for a in vals]
            qdict[key] = np.percentile(vals, qs)
        except KeyError as e:
            print(key, 'not processed yet for', case)
            raise e
        except Exception as e:
            print(case, 'df[', key, ']')
            print(vals)
            raise e

    if plot:
        fig, ax = plot_pdf(case, df=df, keys=keys, **kwargs)
        return qdict, fig, ax
    return qdict


def fit_log(x, h, intercept=False, weights=None, **kwargs):
    def coefficient(x, b):
        return x + b

    try:
        x1 = np.log10(np.array(x))  # this should work for time-series of all x corresponding to h
        h1 = np.log10(np.array(h))
    except TypeError as e:
        print('h', h, type(h))
        print('x', x, type(x))
        raise e
    except ValueError as e:
        print('h', h, type(h))
        print('x', x, type(x))
        raise e
    if intercept:
        try:
            popt, pcov = curve_fit(coefficient, x1, h1)
            slope = 1.0
            intercept = popt[0]
        except ValueError as e:
            print('x', x1, np.shape(x1))
            print('h', h1, np.shape(h1))
            raise e
    elif weights is None:
        try:
            slope, intercept, r_value, p_value, std_err = stats.linregress(x1, h1)
        except ValueError as e:
            print('x', np.shape(x1), 'h', np.shape(h1))
            raise e
    else:
        # use np polyfit for weighted least squares
        intercept, slope = np.polynomial.polynomial.polyfit(x1, h1, deg=1, w=weights)
    return slope, 10 ** intercept


def fit_2log(x, y, h, **kwargs):

    try:
        df = pd.DataFrame({'x':np.log10(np.array(x)), 'y': np.log10(np.array(y)), 'h': np.log10(np.array(h))})
    except (ValueError, TypeError) as e:
        print('h', h, type(h))
        print('x', x, type(x))
        print('y', y, type(y))
        raise e

    X = df[['x','y']]
    Y = df['h']

    # with sklearn
    regr = linear_model.LinearRegression()
    regr.fit(X, Y)

    print('Intercept: \n', regr.intercept_)
    print('Coefficients: \n', regr.coef_)

    return regr.coef_, 10 ** regr.intercept_



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


def plot_h_vs(Ra=None, eta=None, t1_grid=None, end_grid=None, load_grid='auto', data_path=data_path_bullard,
              fig_path=fig_path_bullard, averagescheme=None, p_dimensionals=None,
              fig_fmt='.png', which_x=None, include_regimes=None, regime_grid=None,
              save=True, fname='h',
              labelsize=16, xlabel='', ylabel='dynamic topography', y2label='', title='',
              c_peak='xkcd:forest green', c_rms='xkcd:periwinkle',
              fit=False, logx=True, logy=True, hscale=1, show_isoviscous=False,
              fig=None, ax=None, ylim=None, xlim=None, postprocess_kwargs={}, regime_names=None, **kwargs):
    # either Ra or eta is 1D list of strings (other is singular), t1, end must match shape
    Ra, eta, (t1_grid, load_grid, end_grid, regime_grid) = reshape_inputs(Ra, eta, (t1_grid, load_grid, end_grid, regime_grid))
    if include_regimes is None:
        include_regimes = regime_names
    multifit = False
    if iterable_not_string(which_x) and ('Ra_i' in which_x) or ('Ra_i_eff' in which_x):
        multifit = True
        psuffixes = ['_T', '_h']
        at_sol = True
        postprocess_functions = [T_parameters_at_sol, h_at_ts]
    elif which_x in ['Ra_i', 'Ra_i_eff', 'h_components']:
        psuffixes = ['_T', '_h']
        at_sol = True
        postprocess_functions = [T_parameters_at_sol, h_at_ts]
    elif which_x == 'Ra':
        psuffixes = ['_h_all']
        at_sol = False
        postprocess_functions = [h_at_ts]
    else:
        raise Exception('plot_h_vs(): Invalid variable for x-axis / not implemented')
    if ax is None:
        fig = plt.figure()
        ax = plt.gca()

    if multifit:
        quants = dict.fromkeys(['h_rms', 'h_peak', *which_x])
    else:
        quants = dict.fromkeys(['h_rms', 'h_peak', which_x])

    yx_peak_all, yx_rms_all, n_sols_all = [], [], []

    # loop over cases
    for jj, etastr in enumerate(eta):
        cases, cases_var = get_cases_list(Ra, etastr, end_grid[jj])
        for ii, case in enumerate(cases):
            if regime_grid[jj][ii] in include_regimes:
                t1_ii = t1_grid[jj][ii]
                load_ii = load_grid[jj][ii]
                # dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False, read_statistics=True)

                # load outputs
                dfs = []
                for ip, ps in enumerate(psuffixes):
                    df1 = pickleio(case, suffix=ps, postprocess_functions=postprocess_functions[ip], t1=t1_ii, load=load_ii,
                                   data_path=data_path, at_sol=at_sol, postprocess_kwargs=postprocess_kwargs, **kwargs)
                    dfs.append(df1)
                df = pd.concat(dfs, axis=1)
                df = df.loc[:, ~df.columns.duplicated()]
                df.dropna(axis=0, inplace=True, subset=['h_rms', 'h_peak'])  # double check

                # extract x values for plotting
                if 'h_components' in which_x:
                    if averagescheme == 'timelast':
                        print('    plot_h_vs(): Averaging T components calcualted at each timestep')
                        x = T_components_of_h(case, df=df.mean(axis=0), data_path=data_path, t1=t1_ii, load=load_ii,
                                                         update=False, **postprocess_kwargs, **kwargs)
                    elif averagescheme == 'timefirst':
                        print('    plot_h_vs(): Calculating T components using time-averaged profiles')
                        T_av, y = time_averaged_profile_from_df(df, 'T_av')
                        uv_mag_av, y = time_averaged_profile_from_df(df, 'uv_mag_av')
                        df_av = T_parameters_at_sol(case, n=None, T_av=T_av, uv_mag_av=uv_mag_av, **postprocess_kwargs, **kwargs)
                        x = T_components_of_h(case, df=df_av, data_path=data_path, t1=t1_ii, load=load_ii,
                                                         update=False, **postprocess_kwargs, **kwargs)
                    elif (which_x not in df.columns) or ((which_x in df.columns) and df[which_x].isnull().values.any()):
                        print('    plot_h_vs(): Calculating T components')
                        x = T_components_of_h(case, df=df, data_path=data_path, t1=t1_ii, load=load_ii,
                                                         update=False, **postprocess_kwargs, **kwargs)
                    else:
                        x = df['h_components']
                elif 'Ra_i_eff' in which_x:  # calculate effective Ra using time-mean of T field params
                    if averagescheme == 'timelast':
                        x = Ra_i_eff(Ra_1=float(cases_var[ii]), d_eta=float(eta), T_i=df['T_i'].mean(),
                                     T_l=df['T_l'].mean(), delta_L=df['delta_L'].mean())
                    elif averagescheme == 'timefirst':
                        x = Ra_i_eff(Ra_1=float(cases_var[ii]), d_eta=float(eta), T_i=df_av['T_i'],
                                     T_l=df_av['T_l'], delta_L=df_av['delta_L'])
                    else:
                        raise Exception(
                            'Ra_i_eff not implemented yet if using h output over all timesteps without averaging')
                        # x = Ra_i_eff(Ra_1=float(cases_var[ii]), d_eta=float(eta), T_i=df['T_i'],
                        #              T_l=df['T_l'], delta_L=df['delta_L'])
                elif 'Ra_i' in which_x:
                    if averagescheme == 'timelast':
                        x = Ra_interior(Ra_1=float(cases_var[ii]), d_eta=float(eta), T_i=df['T_i'].mean())
                    elif averagescheme == 'timefirst':
                        x = Ra_interior(Ra_1=float(cases_var[ii]), d_eta=float(eta), T_i=df_av['T_i'])
                    else:
                        raise Exception(
                            'Ra_i not implemented yet if using h output over all timesteps without averaging')
                elif 'Ra' in which_x:
                    if averagescheme == 'timelast' or averagescheme == 'timefirst':
                        x = float(Ra[ii])
                    else:
                        x = float(Ra[ii]) * np.ones(
                            len(df.index))  # normally this is equal to Ra (constant along index)

                if multifit:
                    if 'eta' in which_x:
                        x1 = float(etastr)
                    else:
                        raise Exception(
                            'multifit not implemented yet for second parameter other than delta eta')
                    try:
                        df_plot = pd.DataFrame({which_x[0]: x, which_x[1]: x1})
                    except ValueError:
                        df_plot = pd.DataFrame({which_x[0]: [x], which_x[1]: [x1]})
                else:
                    try:
                        df_plot = pd.DataFrame({which_x: x})
                    except ValueError:
                        df_plot = pd.DataFrame({which_x: [x]})

                # figure out the rest of the plotting stuff, mostly y values, depending on averaging scheme
                if averagescheme == 'timelast':
                    n_sols_all.append(len(df.index))
                    # df_plot = pd.concat([df_plot, df[['h_rms', 'h_peak']]], axis=1)
                    df_plot['h_rms'] = df.h_rms.mean()
                    df_plot['h_peak'] = df.h_peak.mean()
                elif averagescheme == 'timefirst':
                    df_h = pickleio_average(case, suffix='_h_mean', postprocess_fn=h_timeaverage, t1=t1_ii, load=True,
                                            data_path=data_path, postprocess_kwargs=postprocess_kwargs, **kwargs)
                    df_plot = pd.concat([df_plot, df_h], axis=1)
                    n_sols_all.append(1)
                else:
                    # use each xy point (y=h) for fitting to
                    df_plot = pd.concat([df_plot, df], axis=1)
                    # df_plot.dropna(axis=0, how='any', subset=quants.keys(), inplace=True)  # remove any rows with nans
                    n_sols_all.extend([len(df.index)] * len(df.index))

                # append to working
                if multifit:
                    yx_peak_all.append((np.array(df_plot['h_peak'].values) * hscale, np.array(df_plot[[*which_x]].values)))
                    yx_rms_all.append((np.array(df_plot['h_rms'].values) * hscale, np.array(df_plot[[*which_x]].values)))
                else:
                    yx_peak_all.append((np.array(df_plot['h_peak'].values) * hscale, np.array(df_plot[which_x].values)))
                    yx_rms_all.append((np.array(df_plot['h_rms'].values) * hscale, np.array(df_plot[which_x].values)))
                qdict = parameter_percentiles(case, df=df_plot, keys=quants.keys(), plot=False)
                for key in quants.keys():
                    try:
                        quants[key] = np.vstack((quants[key], qdict[key]))
                    except ValueError:  # haven't added anything yet
                        quants[key] = qdict[key].reshape((1,3))  # reshape so it works if you just have one row

    # get errorbars and plot them
    err = dict.fromkeys(quants.keys())
    try:
        for key in quants.keys():
            err[key] = [quants[key][:, 1] - quants[key][:, 0], quants[key][:, 2] - quants[key][:, 1]]
        ax.errorbar(quants[which_x][:, 1], quants['h_peak'][:, 1], yerr=err['h_peak'], xerr=err[which_x], elinewidth=0.5,
                    fmt='d', c=c_peak, alpha=0.8, capsize=5, markeredgecolor=highlight_colour)
        ax.errorbar(quants[which_x][:, 1], quants['h_rms'][:, 1], yerr=err['h_rms'], xerr=err[which_x], elinewidth=0.5,
                    fmt='o', c=c_rms, capsize=5)
    except TypeError:  # no cases in given regimes
        pass

    if fit:
        if which_x == 'h_components':
            intercept=True
        else:
            intercept=False
        ax = fit_cases_on_plot(yx_rms_all, ax, weights=1 / np.array(n_sols_all), c=c_rms, labelsize=labelsize,
                               intercept=intercept, multifit=multifit, **kwargs)

    if show_isoviscous:
        df_JFR = read_JFR('2Dcart_fixed_T_stats_updated.csv', path='/raid1/cmg76/aspect/benchmarks/JFR/')
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
    ax.set_ylabel(ylabel, fontsize=labelsize)
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_title(title, fontsize=labelsize)

    if p_dimensionals is not None:
        ax2 = ax.twinx()
        # ax2.set_ylabel(y2label, fontsize=labelsize)
        ymin, ymax = ax.get_ylim()
        # apply function and set transformed values to right axis limits
        ax2.set_ylim((dimensionalise_h(ymin, p_dimensionals), dimensionalise_h(ymax, p_dimensionals)))
        # set an invisible artist to twin axes
        # to prevent falling back to initial values on rescale events
        ax2.plot([], [])
    if save:
        plot_save(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt)
    return fig, ax


def dimensionalise_h(hprime, p):
    try:
        return hprime * (p['alpha_m'] * p['dT_m'] * p['d_m'])
    except KeyError:
        raise Exception('Need alpha_m, dT_m, and d_m in p_dimensionals to dimensionalise')


def nondimensionalise_h(h, p):
    try:
        return h / (p['alpha_m'] * p['dT_m'] * p['d_m'])
    except KeyError:
        raise Exception('Need alpha_m, dT_m, and d_m in p_dimensionals to nondimensionalise')


def fit_cases_on_plot(yx_all, ax, legend=True, showallscatter=False, weights=None, multifit=False,
                      c='xkcd:periwinkle', legsize=8, legloc='lower left', **kwargs):
    x = [a[1] for a in yx_all]
    y = [a[0] for a in yx_all]
    if np.array(x[0]).ndim > 0 and np.array(y[0]).ndim > 0:
        flatx = [item for sublist in x for item in sublist]
        flaty = [item for sublist in y for item in sublist]
    else:
        flatx, flaty = x, y
    if len(x) > 1:  # can only fit if at least 2 data

        if multifit:
            flatx0 = [a[0] for a in flatx]
            flatx1 = [a[1] for a in flatx]
            x0prime = np.linspace(np.min(flatx0), np.max(flatx0))
            x1prime = np.linspace(np.min(flatx1), np.max(flatx1))
            x0v, x1v = np.meshgrid(x0prime, x1prime)
            const, expon = fit_2log(x=flatx0, y=flatx1, h=flaty)
            hprime = const * x0v ** expon[0] * x1v ** expon[1]
            print('hprime', np.shape(hprime))
        else:
            xprime = np.linspace(np.min(flatx), np.max(flatx))
            expon, const = fit_log(flatx, flaty, weights=weights, **kwargs)
            hprime = const * xprime ** expon
        h3, = ax.plot(xprime, hprime, c=c, ls='--', lw=0.5, zorder=100,
                      label='{:.2e} x^{:.3f}'.format(const, expon))
        if legend:
            handles, labels = ax.get_legend_handles_labels()
            handles.append(h3)
            labels.append('{:.2e} x^{:.3f}'.format(const, expon))
            leg = ax.legend(fontsize=legsize, handles=handles, labels=labels, loc=legloc)
            ax.add_artist(leg)
    else:
        print('    Not enough points to fit')
    if showallscatter:
        ax.scatter(flatx, flaty, c=c, alpha=0.05, s=10)
    return ax


#
# def plot_h_vs_Ra(Ra=None, eta=None, t1=None, fig_path=fig_path_bullard,
#                  showallscatter=False,
#                  save=True, fname='h_vs_Ra', fig_fmt='.png',
#                  labelsize=16, xlabel='', ylabel='dynamic topography', title='',
#                  c_peak='xkcd:forest green', c_rms='xkcd:periwinkle',
#                  fit=False, cases=None, x_var=None, logx=True, logy=True, legend=True,
#                  fig=None, ax=None, ylim=(6e-3, 7e-2), xlim=None, **kwargs):
#     # Ra or eta is list of strings, t1 is a list of numbers the same length
#     if cases is None:
#         cases, x_var = get_cases_list(Ra, eta)
#     if t1 is None:
#         t1 = [0] * len(x_var)
#
#     quants_h_peak = np.zeros((len(x_var), 3))
#     quants_h_rms = np.zeros((len(x_var), 3))
#     x = np.zeros((len(x_var), 3))
#     peak_all = []
#     rms_all = []
#
#     for ii, case in enumerate(cases):
#         df = pickleio(case, suffix='_h_all', postprocess_functions=['h_at_ts'], t1=t1[ii], dat_new=None,
#                       at_sol=False, **kwargs)
#         x[ii, :] = float(x_var[ii])
#         h_peak = df['h_peak']
#         h_rms = df['h_rms']
#         peak_all.append((h_peak, x[ii, 1]))
#         rms_all.append((h_rms, x[ii, 1]))
#
#         qdict = parameter_percentiles(case, df=df, keys=['h_peak', 'h_rms'], plot=False)
#         quants_h_peak[ii, :] = qdict['h_peak']
#         quants_h_rms[ii, :] = qdict['h_rms']
#
#     yerr_peak = [quants_h_peak[:, 1] - quants_h_peak[:, 0], quants_h_peak[:, 2] - quants_h_peak[:, 1]]
#     yerr_rms = [quants_h_rms[:, 1] - quants_h_rms[:, 0], quants_h_rms[:, 2] - quants_h_rms[:, 1]]
#     xerr = None
#
#     if ax is None:
#         fig = plt.figure()
#         ax = plt.gca()
#
#     if fit:
#         if len(x_var) > 1:
#             fitx = [[a[1]] * len(a[0]) for a in rms_all]
#             fith = [a[0] for a in rms_all]
#             flatfitx = [item for sublist in fitx for item in sublist]
#             flatfith = [item for sublist in fith for item in sublist]
#             expon, const = fit_log(flatfitx, flatfith)
#             xprime = [a[1] for a in rms_all]
#             hprime = const * np.array(xprime) ** expon
#             h3, = ax.plot(xprime, hprime, c=c_rms, ls='--', lw=1, zorder=100,
#                           label='{:.2e} Ra^{:.3f}'.format(const, expon))
#             if legend:
#                 ax.legend(
#                     # handles=[h3], labels=[],
#                     loc='lower left')  # always show what fit is
#         else:
#             print('    Not enough points to fit')
#
#     ax.errorbar(x[:, 1], quants_h_peak[:, 1], yerr=yerr_peak, xerr=xerr,
#                 fmt='^', c=c_peak, alpha=0.9, capsize=5)
#     ax.errorbar(x[:, 1], quants_h_rms[:, 1], yerr=yerr_rms, xerr=xerr,
#                 fmt='o', c=c_rms, alpha=0.9, capsize=5)
#
#     if logx:
#         ax.set_xscale('log')
#     if logy:
#         ax.set_yscale('log')
#     if ylim is not None:
#         ax.set_ylim(ylim[0], ylim[1])  # for fair comparison
#     if xlim is not None:
#         ax.set_xlim(xlim)
#     ax.set_ylabel(ylabel, fontsize=labelsize)
#     ax.set_xlabel(xlabel, fontsize=labelsize)
#     ax.set_title(title, fontsize=labelsize)
#
#     if save:
#         savefig(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt)
#     return fig, ax

#
# def plot_h_vs_components(Ra=None, eta=None, t1=None, data_path=data_path_bullard, fig_path=fig_path_bullard,
#                          fig_fmt='.png',
#                          save=True, fname='h_T', showallscatter=False,
#                          labelsize=16, xlabel=r'$\delta_rh \Delta T_{rh}$', ylabel='dynamic topography', title='',
#                          c_peak='xkcd:forest green', c_rms='xkcd:periwinkle', legend=True,
#                          fit=False, cases=None, x_var=None, logx=True, logy=True,
#                          fig=None, ax=None, ylim=(6e-3, 7e-2), xlim=None, **kwargs):
#     # Ra or eta is list of strings, t1 is a list of numbers the same length
#     # instead of plotting vs Ra or eta, plot vs theoretical components of scaling relationship
#
#     if cases is None:
#         cases, x_var = get_cases_list(Ra, eta)
#     if t1 is None:
#         t1 = [0] * len(cases)
#
#     quants_h_peak = np.zeros((len(x_var), 3))
#     quants_h_rms = np.zeros((len(x_var), 3))
#     quants_h_components = np.zeros((len(x_var), 3))
#     peak_all = []
#     rms_all = []
#
#     for ii, case in enumerate(cases):
#
#         pickle_remove_duplicate(case, suffix='_T', which='sol', data_path=data_path)
#
#         dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False, read_statistics=True)
#
#         # load h and T
#         df1 = pickleio(case, suffix='_h', postprocess_functions=[h_at_ts], t1=t1[ii], dat_new=dat,
#                        data_path=data_path, at_sol=True, **kwargs)
#         df2 = pickleio(case, suffix='_T', postprocess_functions=[T_parameters_at_sol], t1=t1[ii],
#                        dat_new=dat, data_path=data_path, at_sol=True, **kwargs)
#
#
#         try:
#             df = pd.concat([df1, df2], axis=1)  # concatenate along columns
#         except ValueError as e:
#             # print('\n\ndf1', df1)
#             # print('\n\ndf2', df2)
#             print('\ndf1 sol range', np.min(df1['sol']), np.max(df1['sol']))
#             print('\ndf2 sol range', np.min(df2['sol']), np.max(df2['sol']))
#             print(case)
#             raise e
#
#         try:
#             h_components = df['h_components']
#         except KeyError:
#             h_components = T_components_of_h(case, df=df, dat=dat, data_path=data_path, t1=t1[ii], update=True,
#                                              **kwargs)
#             df['h_components'] = h_components
#
#         h_peak = df['h_peak']
#         h_rms = df['h_rms']
#         # old way of accounting for loading all h instead of just at sols
#         # t1_idx = np.argmax(dat.stats_time > t1[ii])  # timestep corresponding to t1
#         # sol_idx = dat.find_time_at_sol()
#         # sol_idx = np.array(sol_idx[sol_idx >= t1_idx])  # cut any values of idx below t1_idx
#         # # account for the stored peak_list and rms_list starting at t1 for indexing
#         # if np.shape(sol_idx) != np.shape(h_components):
#         #     print('inconsistent times: sol_idx', np.shape(sol_idx), 'delta', np.shape(df['delta_rh']))
#         # try:
#         #     peak_list = [peak_list[j - t1_idx] for j in sol_idx]
#         #     rms_list = [rms_list[j - t1_idx] for j in sol_idx]
#         # except IndexError:
#         #     print('sol_idx - t1_idx', [j - t1_idx for j in sol_idx])
#         peak_all.append((h_peak, h_components))
#         rms_all.append((h_rms, h_components))
#
#         qdict = parameter_percentiles(case, df=df, keys=['h_peak', 'h_rms', 'h_components'], plot=False)
#         quants_h_peak[ii, :] = qdict['h_peak']
#         quants_h_rms[ii, :] = qdict['h_rms']
#         quants_h_components[ii, :] = qdict['h_components']
#
#     yerr_peak = [quants_h_peak[:, 1] - quants_h_peak[:, 0], quants_h_peak[:, 2] - quants_h_peak[:, 1]]
#     yerr_rms = [quants_h_rms[:, 1] - quants_h_rms[:, 0], quants_h_rms[:, 2] - quants_h_rms[:, 1]]
#     xerr = [quants_h_components[:, 1] - quants_h_components[:, 0],
#             quants_h_components[:, 2] - quants_h_components[:, 1]]
#
#     if ax is None:
#         fig = plt.figure()
#         ax = plt.gca()
#
#     fitx = [a[1] for a in rms_all]
#     fith_rms = [a[0] for a in rms_all]
#     flatfitx = [item for sublist in fitx for item in sublist]
#     flatfith_rms = [item for sublist in fith_rms for item in sublist]
#     fith_peak = [a[0] for a in peak_all]
#     flatfith_peak = [item for sublist in fith_peak for item in sublist]
#
#     if fit:
#         #         fitx = x_var
#         #         fitidx = np.where(np.intersect1d(x[:,1], fitx))[0] # only fit subset?
#         if len(x_var) > 1:  # can only fit if at least 2 data
#             expon, const = fit_log(flatfitx, flatfith_rms)
#             xprime = np.linspace(np.min(flatfitx), np.max(flatfitx))
#             hprime = const * xprime ** expon
#             h3, = ax.plot(xprime, hprime, c=c_rms, ls='--', lw=1, zorder=100,
#                           label='{:.2e} x^{:.3f}'.format(const, expon))
#             if legend:
#                 ax.legend(
#                     # handles=[h3], labels=[],
#                     loc='lower left')
#         else:
#             print('    Not enough points to fit')
#
#     ax.errorbar(quants_h_components[:, 1], quants_h_peak[:, 1], yerr=yerr_peak, xerr=xerr,
#                 fmt='^', c=c_peak, alpha=0.9, capsize=5)
#     ax.errorbar(quants_h_components[:, 1], quants_h_rms[:, 1], yerr=yerr_rms, xerr=xerr,
#                 fmt='o', c=c_rms, alpha=0.9, capsize=5)
#
#     if showallscatter:
#         ax.scatter(flatfitx, flatfith_rms, c=c_rms, alpha=0.1, s=20)
#         ax.scatter(flatfitx, flatfith_peak, c=c_peak, alpha=0.1, s=20)
#
#     if logx:
#         ax.set_xscale('log')
#     if logy:
#         ax.set_yscale('log')
#     if ylim is not None:
#         ax.set_ylim(ylim[0], ylim[1])  # for fair comparison
#     if xlim is not None:
#         ax.set_xlim(xlim)
#     ax.set_ylabel(ylabel, fontsize=labelsize)
#     ax.set_xlabel(xlabel, fontsize=labelsize)
#     ax.set_title(title, fontsize=labelsize)
#
#     if save:
#         savefig(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt)
#     return fig, ax


def subplots_topo_regimes(Ra_ls, eta_ls, regime_grid, regime_names, c_regimes=None, save=True, t1_grid=None, nrows=2,
                          ncols=2, which_x='Ra', leftleg_bbox=(-0.05, 1), p_dimensionals=None,
                          load_grid='auto', fig_path=fig_path_bullard, fname='h_Ra_all', fig_fmt='.png', end_grid=None,
                          show_bounds=False, regimes_title='', Ra_i=False, show_isoviscous=False, y2label='',
                          labelsize=14, xlabel='Ra', include_regimes=None, ylabel='dynamic topography', xlabelpad=12, ylabelpad=2, **kwargs):
    Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid) = reshape_inputs(Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid))
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
                                    load_grid=load_ii[Ra_regime_idx], regime_grid=regime_grid[ii,Ra_regime_idx],
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
                        fig_path=fig_path_bullard, load_grid='auto', regime_grid=None, include_regimes=None, Ra_i=False, compare_exponent=None,
                        save=True, fname='Ra_scalings', labelsize=16, ylabels=None, psuffixes='', title='',
                        postprocess_functions=[], xlim=None, ylim=None, legloc=None, averagescheme=None, regime_names=None,
                        cmap='magma', compare_pub=None, compare_label=None, vmin=None, vmax=None,
                        fig=None, axes=None, fig_fmt='.png', postprocess_kwargs={}, **kwargs):
    # Ra or eta is list of strings, t1 is a list of numbers the same length
    # instead of plotting vs Ra or eta, plot vs theoretical components of scaling relationship

    Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid) = reshape_inputs(Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid))
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
        cases, Ra_var = get_cases_list(Ra_ls, eta_str, end_grid[jj])
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
                    dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False,
                                           read_statistics=True)
                else:
                    dat = None
                dfs = []
                for ip, suffix in enumerate(psuffixes):
                    df1 = pickleio(case, suffix=suffix, postprocess_functions=postprocess_functions[ip], t1=t1_ii,
                                   dat_new=dat, data_path=data_path, load=load_ii, postprocess_kwargs=postprocess_kwargs, **kwargs)
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
                    T_av, y = time_averaged_profile_from_df(df, 'T_av')
                    uv_mag_av, y = time_averaged_profile_from_df(df, 'uv_mag_av')
                    df_av = T_parameters_at_sol(case, n=None, T_av=T_av, uv_mag_av=uv_mag_av, **postprocess_kwargs,
                                                **kwargs)
                    # h_components = T_components_of_h(case, df=df_av, data_path=data_path, t1=t1_ii, load=load_ii,
                    #                                  update=False, **postprocess_kwargs, **kwargs)

                    if 'Nu' in keys:
                        df_av['Nu'] = np.mean(df['Nu'])  # using average flux anyways so scheme doesn't matter
                    df = df_av

                for key in keys:
                    med = np.median(df[key])
                    if np.isnan(med):
                        raise Exception('NaN in median, key:', key, '\n', df[key])
                    plot_data[key].append(np.median(df[key]))
                if Ra_i == 'eff':
                    plot_data['Ra'].append(Ra_i_eff(Ra_1=Ra_ii, d_eta=float(eta_str), T_i=np.median(df['T_i']),
                                                    T_l=np.median(df['T_l']), delta_L=np.median(df['delta_L'])))
                elif Ra_i:
                    plot_data['Ra'].append(Ra_interior(Ra_1=Ra_ii, d_eta=float(eta_str), T_i=np.median(df['T_i'])))
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


def moresi95(Ra=None, d_eta=None, df=None, dat=None, case=None, T1=1, dT=1,
             data_path=data_path_bullard, load='auto', **kwargs):
    if df is None:
        df = pickleio(case, postprocess_functions=[T_parameters_at_sol], suffix='_T', dat_new=dat, data_path=data_path,
                      load=load)

    T_i = np.median(df['T_i'])
    gamma = np.log(d_eta)  # gamma for this delta eta
    p = gamma * dT
    eta_1 = np.exp(-gamma * T1)
    eta_i = np.exp(-gamma * T_i)
    Ra_i = np.array(Ra) * eta_1 / eta_i
    delta_1 = 0.58 * p ** 0.29 * Ra_i ** -0.24
    Nu = 1.89 * p ** -1.02 * Ra_i ** 0.2
    T_i_scaling = 1 - (1.1 * p ** -0.73 * Ra_i ** -0.04)
    delta_0 = T_i_scaling / Nu
    return {'Ra_i': Ra_i, 'delta_0': delta_0, 'Nu': Nu, 'T_i': T_i_scaling, 'delta_1': delta_1}


def Ra_interior(Ra_1=None, d_eta=None, T_i=None, T1=1, T0=0):
    theta = np.log(d_eta)  # gamma for this delta eta
    gamma = theta / (np.array(T1) - np.array(T0))
    eta_1 = np.exp(-gamma * np.array(T1))  # bottom viscosity
    eta_i = np.exp(-gamma * np.array(T_i))
    Ra_i = np.array(Ra_1) * eta_1 / eta_i
    return Ra_i


def Ra_eff(Ra=None, T1=1, T0=0, T_l=None, Z=1, delta_L=None):
    return Ra * (T1 - T_l) / (T1 - T0) * ((Z - np.array(delta_L)) / Z) ** 3


def Ra_i_eff(Ra_1=None, d_eta=None, T_i=None, T1=1, T0=0, T_l=None, Z=1, delta_L=None):
    # accounting for effective depth of convection with large lid
    # order of Ra_i and Ra_eff matters if updating T1 etc
    Ra_i = Ra_interior(Ra_1, d_eta=d_eta, T_i=T_i, T1=T1, T0=T0)
    return Ra_eff(Ra_i, T1=T1, T_l=T_l, Z=Z, delta_L=delta_L)


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
        if os.path.exists(data_path + 'output-' + case) and not (regime_grid[ii] == 'no convection') and not (regime_grid[ii] == 'sluggish'):
            print('Plotting summary for', case, 'using t1 =', t1[ii])
            t1_ii = t1[ii]
            if iterable_not_string(load):
                load_ii = load[ii]
            else:
                load_ii = load
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
                if t1[ii] < 1:
                    sol_df = pickleio(case, suffix='_T', postprocess_functions=[T_parameters_at_sol], t1=t1_ii,
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
                        sol_df = pickleio(case, suffix='_T', postprocess_functions=[T_parameters_at_sol], t1=t1_ii,
                                          dat_new=dat, load=load_ii, data_path=data_path, fig_path=fig_path, **kwargs)

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
                    ts_df = pickleio(case, suffix='_h_all', postprocess_functions=[h_at_ts], t1=t1_ii,
                                     dat_new=dat, load=load_ii, data_path=data_path, fig_path=fig_path,
                                     at_sol=False, **kwargs)

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


def plot_parameter_grid(Ra, eta, function, data_path=data_path_bullard, fig_path=fig_path_bullard, load='auto',
                        vmin=None, vmax=None, set_under=None, set_over=None,
                        save=True, fname='grid', labelsize=16, fig_fmt='.png', t1=None, end=None, cticklabels=None,
                        cticks=None, title='', lognorm=False, log=False, clabel=None, discrete=False,
                        overplot_h=False, nlevels_contour=10, cmap='jet', clist=None, cmap_contours='spring', **kwargs):
    # plot output of any (scalar-returning) function (with optional topo contours?)

    if t1 is None:
        t1 = np.zeros((len(eta), len(Ra)))
    if not_iterable(load):  #
        load = np.array([[load] * len(Ra)] * len(eta))
    fig, ax = plt.subplots(1, 1)
    plot_grid = np.zeros((len(eta), len(Ra)))

    for jj, eta_str in enumerate(eta):
        cases, _ = get_cases_list(Ra, eta_str, end[jj])
        for ii, Ra_str in enumerate(Ra):
            # calculate value at this parameter-space coordinate
            if os.path.exists(data_path + 'output-' + cases[ii] + '/'):
                plot_grid[jj, ii] = function(Ra=Ra_str, eta=eta_str, ii=ii, jj=jj, case=cases[ii], load=load[jj][ii],
                                             t1=t1[jj][ii], data_path=data_path, **kwargs)
            else:
                plot_grid[jj, ii] = np.nan

    if log:
        plot_grid = np.log10(plot_grid)
    m = np.ma.masked_where(np.isnan(plot_grid), plot_grid)

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

    im = ax.imshow(m, origin='bottom', aspect='equal', interpolation='None', cmap=cmap, vmin=vmin, vmax=vmax,
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

    if overplot_h:
        if t1 is None:
            t1 = np.zeros((len(eta), len(Ra)))
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
        CS = ax.contour(h_grid, nlevels_contour, cmap=cmap_contours)
        ax.clabel(CS, inline=1, fontsize=labelsize)

    ax.set_title(title, fontsize=labelsize)
    if save:
        plot_save(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt, tight_layout=False)


def regime_to_digital(ii=None, jj=None, regime_grid=None, regime_names=None, **kwargs):
    label = regime_grid[jj, ii]
    digi = np.nonzero(np.array(regime_names) == label)[0]
    if not list(digi):
        return np.nan
    else:
        return digi[0] + 1


def surf_mobility_at_sol(case=None, dat=None, n=None, data_path=data_path_bullard, **kwargs):
    if dat is None:
        dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False,
                               read_statistics=False, read_parameters=False)
    if n is None:
        n = dat.final_step()

    df_T = pickleio(case, '_T', [T_parameters_at_sol], dat_new=dat, at_sol=True,
                    data_path=data_path, **kwargs)
    if not df_T.empty:
        df_sol = df_T.set_index('sol')
        S = dat.surface_mobility(n=n, delta_0=df_sol.loc[n, 'delta_0'], delta_rh=df_sol.loc[n, 'delta_rh'],
                                 delta_l=df_sol.loc[n, 'delta_L'])
        return S
    else:
        return np.nan


def plot_velocity_profile(case, dat=None, n=None, xlabel='rms velocity', ylabel='depth', fig=None, ax=None,
                          data_path=data_path_bullard,
                          labelsize=16, fig_path=fig_path_bullard, fname='velocity', save=True, fig_fmt='.png',
                          **kwargs):
    if fig is None:
        fig, ax = plt.subplots(figsize=(4, 4))

    if dat is None:
        dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False,
                               read_statistics=False, read_parameters=False)
    if n is None:
        n = dat.final_step()

    x, y, _, u, v, _, mag = dat.read_velocity(n=n, verbose=False)
    mag_av = post.horizontal_mean(mag, x)
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
        T_params = pickleio(case, suffix='_T', postprocess_functions=[T_parameters_at_sol], t1=t1,
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


def subplots_evol_at_sol(Ra_ls, eta_ls, regime_grid=None, save=True, t1_grid=None,
                         load_grid='auto', psuffixes=['_T'], postprocess_functions=[T_parameters_at_sol],
                         fig_path=fig_path_bullard, fname='evol', fig_fmt='.png', end_grid=None, normtime=True,
                         labelsize=14, xlabel=r'Time', ylabels=None, keys=None, title='', legsize=10,
                         xlabelpad=8, ylabelpad=-2, markers=None, markersize=20, alpha=0.5,
                         fig=None, cmap='magma', vmin=None, vmax=None, include_regimes=None, regime_names=None,
                         data_path=data_path_bullard, **kwargs):
    # plot time-evolution of list of keys for all cases in given regime

    Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid) = reshape_inputs(Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid))
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
    c_list = colorize(logeta_fl, cmap=cmap, vmin=vmin, vmax=vmax)[0]

    for jj, eta_str in enumerate(eta_ls):
        cases, Ra_var = get_cases_list(Ra_ls, eta_str, end_grid[jj])
        c_jj = c_list[jj]

        for ii, case in enumerate(cases):
            t1_ii = t1_grid[jj][ii]
            load_ii = load_grid[jj][ii]
            marker_ii = markers[ii]

            if (t1_ii != 1) and (os.path.exists(data_path + 'output-' + case)) and (
                    regime_grid[jj][ii] in include_regimes):
                Ra_ii = float(Ra_var[ii])

                # load data
                if load_ii == 'auto':
                    dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', verbose=False,
                                           read_statistics=True)
                else:
                    dat = None
                dfs = []
                for ip, suffix in enumerate(psuffixes):
                    df1 = pickleio(case, suffix=suffix, postprocess_functions=postprocess_functions[ip], t1=t1_ii,
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
                    ax.scatter(x_data, y_data, c=c_jj, s=markersize, marker=marker_ii, alpha=alpha)
                    ax.plot(x_data, y_data, c=c_jj, lw=0.8, alpha=alpha)

    # legend proxy artist
    ax = axes[0]
    # scat = ax.scatter(logeta_fl, logeta_fl, visible=False, c=np.array(logeta_fl), cmap=cmap, s=markersize,
    #                   vmin=vmin, vmax=vmax)  # dummy - neeeds matplotlib 3.1.1
    # legend1 = ax.legend(*scat.legend_elements(num=len(logeta_fl)),
    #                     loc="upper left", title=r"log $\Delta \eta$", fontsize=legsize)
    lines = []
    for jj, leta in enumerate(logeta_fl):
        p = mlines.Line2D([], [], lw=0, color=c_list[jj], marker='o', markersize=markersize/4, label=leta, alpha=alpha)
        lines.append(p)
    legend1 = ax.legend(lines, [l.get_label() for l in lines], fontsize=legsize, frameon=True, loc="upper left",
                        title=r"log $\Delta \eta$",)
    ax.add_artist(legend1)

    ax = axes[-1]
    lines = []
    for ii, Ra in enumerate(Ra_ls):
        p = mlines.Line2D([], [], lw=0, color='k', marker=markers[ii], markersize=markersize/4, label=Ra, alpha=alpha)
        lines.append(p)
    legend2 = ax.legend(lines, [l.get_label() for l in lines], fontsize=legsize, frameon=True, loc="lower right",
                        title="Ra", )
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
    time, y = read_evol(case, col, dat=dat)
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
        plot_save(fig, case + '_h_' + '{:05}'.format(ts), fig_path=fig_path, fig_fmt=fig_fmt)


def read_JFR(fname, path='/raid1/cmg76/aspect/benchmarks/JFR/'):
    df = pd.read_csv(path + fname, header=0, index_col=False)
    print('Loaded', fname, df.columns)
    Ra = np.array(df.Ra)
    h_peak = np.array(df.peak_topo)
    h_rms = np.array(df.RMS_topo)
    Nu = np.array(df.Nu)
    return df  # Ra, h_peak, h_rms, Nu


def Nu_eff(gamma=None, d_m=None, delta_L=None, alpha_m=None, g=None, b=None, kappa=None, T_i=None, n=1, k=1, **kwargs):
    return None


def reprocess_all_at_sol(Ra_ls, eta_ls, psuffixes, postprocess_functions, t1_grid=None, end_grid=None,
                         data_path=data_path_bullard, redo=True, load_grid=None, regime_grid=None, include_regimes=None, **kwargs):
    Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid) = reshape_inputs(Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid))

    for jj, eta_str in enumerate(eta_ls):
        cases, Ra_var = get_cases_list(Ra_ls, eta_str, end_grid[jj])
        for ii, case in enumerate(cases):
            t1_ii = t1_grid[jj][ii]
            if (t1_ii != 1) and (os.path.exists(data_path + 'output-' + case)):
                if include_regimes is not None and (regime_grid[jj][ii] in include_regimes):
                    print(case)
                    if redo:
                        # for recalculating everything if you fucked up e.g.
                        load = False
                    elif load_grid is not None:
                        load = load_grid[jj][ii]
                    else:
                        load = 'auto'
                    for ip, suffix in enumerate(psuffixes):
                        pickleio(case, suffix=suffix, postprocess_functions=postprocess_functions[ip], t1=t1_ii,
                                 data_path=data_path, at_sol=True, load=load, **kwargs)


def pickleio_average(case, postprocess_fn=None, t1=0, load=True, suffix='', data_path=data_path_bullard, fend='.pkl', **kwargs):
    # for these time-average ones 'auto' counts as reprocess
    case_path = data_path + 'output-' + case + '/'
    fname = case + suffix + fend
    if load is not True:
        dat = post.Aspect_Data(directory=case_path, verbose=False,
                               read_statistics=True, read_parameters=False)
        time = dat.stats_time
        i_time = np.argmax(time >= t1)  # index of first timestep to process
        ts0 = i_time
        tsf = len(time) - 1
        df = postprocess_fn(case, ts0, tsf, **kwargs)
        print('Processed', postprocess_fn, 'for time steps', ts0, 'to', tsf)
        pkl.dump(df, open(case_path + 'pickle/' + fname, "wb"))
    else:
        df = pkl.load(open(case_path + 'pickle/' + fname, "rb"))
    return df


def reprocess_all_average(Ra_ls, eta_ls, t1_grid=None, end_grid=None,
                          data_path=data_path_bullard, redo=True, load_grid=None, regime_grid=None, include_regimes=None, **kwargs):
    Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid) = reshape_inputs(Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid))

    for jj, eta_str in enumerate(eta_ls):
        cases, Ra_var = get_cases_list(Ra_ls, eta_str, end_grid[jj])
        for ii, case in enumerate(cases):
            t1_ii = t1_grid[jj][ii]
            if (t1_ii != 1) and (os.path.exists(data_path + 'output-' + case)):
                if include_regimes is not None and (regime_grid[jj][ii] in include_regimes):
                    print('Found', case)
                    if redo:
                        # for recalculating everything if you fucked up e.g.
                        load = False
                    elif load_grid is not None:
                        load = load_grid[jj][ii]
                    else:
                        load = 'auto'
                    df = pickleio_average(case, t1=t1_ii, load=load, data_path=data_path, **kwargs)
                    print(df)


# def reprocess_time_averages(Ra_ls, eta_ls, psuffixes, postprocess_functions, t1_grid=None, end_grid=None,
#                             data_path=data_path_bullard, redo=True, load_grid=None, **kwargs):
#     Ra_ls, eta_ls, t1_grid, load_grid, end_grid = reshape_inputs(Ra_ls, eta_ls, t1_grid, load_grid, end_grid)
#
#     for jj, eta_str in enumerate(eta_ls):
#         cases, Ra_var = get_cases_list(Ra_ls, eta_str, end_grid[jj])
#         for ii, case in enumerate(cases):
#             t1_ii = t1_grid[jj][ii]
#             if (t1_ii != 1) and (os.path.exists(data_path + 'output-' + case)):
#                 print(case)
#                 if redo:
#                     # for recalculating everything if you fucked up e.g.
#                     load = False
#                 elif load_grid is not None:
#                     load = load_grid[jj][ii]
#                 else:
#                     load = 'auto'
#                 for ip, suffix in enumerate(psuffixes):
#                     pickleio(case, suffix=suffix, postprocess_functions=postprocess_functions[ip], t1=t1_ii,
#                              data_path=data_path, at_sol=True, load=load, **kwargs)


def plot_heuristic_scalings(Ra_ls, eta_ls, regime_grid=None, t1_grid=None, load_grid=None, end_grid=None,
                            literature_file=None,
                            legend=True, postprocess_kwargs={}, regime_names=None,
                            c='k', averagescheme=None, ylim=None, which_h='rms', data_path=data_path_bullard,
                            save=True, fname='model-data', labelsize=16, clist=None,
                            cmap='magma', cbar='eta', include_regimes=None, **kwargs):

    Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid) = reshape_inputs(Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid))
    if include_regimes is None:
        include_regimes = regime_names
    if averagescheme == 'timefirst':
        psuffixes = ['_T']
        postprocess_functions = [T_parameters_at_sol]
    else:
        psuffixes = ['_T', '_h']
        postprocess_functions = [T_parameters_at_sol, h_at_ts]

    if cbar == 'eta':
        clabel = r'$\Delta \eta$'
        cticklabels = None
        vmin, vmax = 0.9e5, 2e8
        crot = 0
        cnorm = LogNorm()
        discrete = False
    elif cbar == 'regime':
        clabel = ''  # 'Stationarity'
        cticklabels = ['steady', 'transitional', 'chaotic']
        vmin, vmax = 1, 3
        cnorm = None
        crot = 0  # 70
        cmap = cmap_from_list(clist, cmap_name='regimes')
        discrete = True

    h_data_all = []
    x_data_all = []
    c_data_all = []
    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    for jj, eta_str in enumerate(eta_ls):
        cases, Ra_var = get_cases_list(Ra_ls, eta_str, end_grid[jj])
        for ii, case in enumerate(cases):
            t1_ii = t1_grid[jj][ii]
            load_ii = load_grid[jj][ii]
            if (t1_ii != 1) and (os.path.exists(data_path + 'output-' + case)) and (
                    regime_grid[jj][ii] in include_regimes):
                # load outputs
                dfs = []
                for ip, ps in enumerate(psuffixes):
                    df1 = pickleio(case, suffix=ps, postprocess_functions=postprocess_functions[ip], t1=t1_ii,
                                   load=load_ii, data_path=data_path, at_sol=True, **kwargs)
                    dfs.append(df1)
                df = pd.concat(dfs, axis=1)
                df = df.loc[:, ~df.columns.duplicated()]

                if averagescheme == 'timelast':
                    print('    plot_heuristic_scalings(): Averaging T components calculated at each timestep')
                    df = df.dropna(axis=0, how='any',
                                   subset=['h_peak', 'h_rms', 'h_components']).mean(axis=0)  # remove any rows with nans
                    h_components = T_components_of_h(case, df=df, data_path=data_path, t1=t1_ii,
                                                     load=load_ii, update=False, **postprocess_kwargs, **kwargs)

                elif averagescheme == 'timefirst':
                    print('    plot_heuristic_scalings(): Calculating T components using time-averaged profiles')
                    T_av, y = time_averaged_profile_from_df(df, 'T_av')
                    uv_mag_av, y = time_averaged_profile_from_df(df, 'uv_mag_av')
                    df_av = T_parameters_at_sol(case, n=None, T_av=T_av, uv_mag_av=uv_mag_av, **postprocess_kwargs,
                                                **kwargs)
                    h_components = T_components_of_h(case, df=df_av, data_path=data_path, t1=t1_ii, load=load_ii,
                                                     update=False, **postprocess_kwargs, **kwargs)
                    df = pickleio_average(case, suffix='_h_mean', postprocess_fn=h_timeaverage, t1=t1_ii, load=True,
                                            data_path=data_path, postprocess_kwargs=postprocess_kwargs, **kwargs)

                if which_h == 'rms':
                    h = np.array(df['h_rms'])
                elif which_h == 'peak':
                    h = np.array(df['h_peak'])

                h_data_all.append((np.mean(h)))
                x_data_all.append(np.mean(h_components))
                if cbar == 'eta':
                    c_data_all.append(float(eta_str))
                elif cbar == 'regime':
                    c_data_all.append(regime_to_digital(ii, jj, regime_grid=regime_grid, regime_names=regime_names))
    x_data, h_data = [list(tup) for tup in zip(*sorted(zip(x_data_all, h_data_all)))]  # sort according to x
    print('x_data: min', np.min(x_data), 'max', np.max(x_data), '(n =', len(x_data), ')')
    expon, const = fit_log(x_data, h_data, weights=None, **kwargs)
    xprime = np.linspace(np.min(x_data), np.max(x_data))
    h_fit = const * np.array(x_data) ** expon
    print('(data, fit)', [(da, fi) for (da, fi) in zip(h_data, h_fit)])
    if not cbar:
        c = 'k'
        vmin, vmax = None, None
    else:
        c = [q for _, q in sorted(zip(x_data_all, c_data_all))]
    ax.set_ylabel('Model', fontsize=labelsize)
    ax.set_xlabel('Data', fontsize=labelsize)
    title = 'Fit to h = ({:.2e}'.format(const) + r') $\alpha \Delta T_{rh} \delta_{rh}$' + '^{:.3f}'.format(expon)
    ax.set_title(title, fontsize=labelsize)
    ax.set_xscale('log')
    ax.set_yscale('log')
    if ylim is not None:
        ax.set_ylim(ylim)
        ax.set_xlim(ylim)
    else:
        ax.axis('equal')

    scat = ax.scatter(h_data, h_fit, s=30, zorder=100, c=c, cmap=cmap, norm=cnorm, vmin=vmin, vmax=vmax)

    if not (not cbar):
        cbar = colourbar(scat, label=clabel, ticklabels=cticklabels, labelsize=labelsize, discrete=discrete,
                         vmin=vmin, vmax=vmax, rot=crot)

    fig, ax = plot_error_contours(fig, ax)

    if save:
        plot_save(fig, fname, **kwargs)
    return const, expon


def plot_error_contours(fig, ax, errs=None, c='k', labels=True):
    if errs is None:
        errs = [0.5, 0.2, 0.1]
    x0 = np.array(ax.get_xlim())
    y0 = np.array(ax.get_ylim())
    # set 1:1 line
    ax.plot(x0, y0, c=c, lw=2)
    for err in errs:
        y_plus = y0 + err * y0
        y_minus = y0 - err * y0
        for y in [y_plus, y_minus]:
            l, = ax.plot(x0, y, c=c, lw=1, ls='--')
            if labels:
                pos = [(x0[-2] + x0[-1]) / 3., (y[-2] + y[-1]) / 3.]
                # transform data points to screen space
                xscreen = ax.transData.transform(np.array((x0[-2::], y[-2::])))
                rot = np.rad2deg(np.arctan2(*np.abs(np.gradient(xscreen)[0][0][::-1])))
                ltex = ax.text(pos[0], pos[1], '{0:.0f}%'.format(err * 100), size=9, rotation=rot, color=l.get_color(),
                               ha="center", va="center", bbox=dict(boxstyle='square,pad=-0.0', ec='1', fc='1'))
    return fig, ax

# def subplots_h_fit_2D(Ra_ls=None, eta_ls=None, t1=None, end=None, load='auto', data_path=data_path_bullard,
#               fig_path=fig_path_bullard, alpha_m=1.35e-5,
#               fig_fmt='.png', regime_grid=None, nrows=2, ncols=2,
#               save=True, fname='h-2d', xlabelpad=10, ylabelpad=10,
#               labelsize=16, xlabel=r'$\delta_{rh}$', ylabel='$h_{rms}$',  title='',
#               c_rms='xkcd:periwinkle',
#                logx=True, logy=True,
#               ax=None, ylim=None, xlim=None, **kwargs):
#
#     if t1 is None:
#         t1 = [[0] * len(Ra_ls)] * len(eta_ls)
#     if not_iterable(load):  #
#         load = np.array([[load] * len(Ra_ls)] * len(eta_ls))
#
#     fig = plt.figure(figsize=(7, 7))
#     bigax = fig.add_subplot(111)  # The big subplot
#     bigax.spines['top'].set_color('none')
#     bigax.spines['bottom'].set_color('none')
#     bigax.spines['left'].set_color('none')
#     bigax.spines['right'].set_color('none')
#     bigax.tick_params(labelcolor='w', which='both', bottom=False, left=False, right=False, top=False)
#     bigax.set_xlabel(xlabel, fontsize=labelsize, labelpad=xlabelpad)
#     bigax.set_ylabel(ylabel, fontsize=labelsize, labelpad=ylabelpad)
#
#     psuffixes = ['_T', '_h']
#     postprocess_functions = [T_parameters_at_sol, h_at_ts]
#
#     h_data = np.zeros_like(load)
#     x_data = np.zeros_like(load)
#     y_data = np.zeros_like(load)
#
#     for jj, eta_str in enumerate(eta_ls):
#         cases, Ra = get_cases_list(Ra_ls, eta_str, end[jj])
#         for ii, case in enumerate(cases):
#             t1_ii = t1[jj][ii]
#             load_ii = load[jj][ii]
#             if (t1_ii != 1) and (os.path.exists(data_path + 'output-' + case)) and (regime_grid[jj][ii] != 'sluggish'):
#                 # load outputs
#                 dfs = []
#                 for ip, ps in enumerate(psuffixes):
#                     df1 = pickleio(case, suffix=ps, postprocess_functions=postprocess_functions[ip], t1=t1_ii,
#                                    load=load_ii, data_path=data_path, at_sol=True, **kwargs)
#                     dfs.append(df1)
#                 df = pd.concat(dfs, axis=1)
#                 df = df.loc[:, ~df.columns.duplicated()]
#                 df = df.dropna(axis=0, how='any', thresh=None, subset=['h_rms', 'dT_rh', 'delta_rh'])
#                 h_data[jj,ii] = df.h_rms.mean()
#                 x_data[jj,ii] = df.dT_rh.mean()
#                 y_data[jj,ii] = df.delta_rh.mean()
#
#     coef_x, coef_y, C = fit_2d_log(x_data.flatten(), y_data.flatten(), h_data.flatten())
#
#     y_sample = np.mean(y_data, axis=0)
#     print('y_sample', np.shape(y_sample))
#     for jj, eta_ii in enumerate(eta_ls):
#         z = int(str(nrows) + str(ncols) + str(jj + 1))
#         ax = fig.add_subplot(z)
#
#         x_data = x_data[jj]
#         print('jj', jj)
#         print('x_data', np.shape(x_data))
#         y_data = y_sample[jj]
#         h_data = h_data[jj]
#         h_model = C * x_data**coef_x * y_data**coef_y
#
#         ax.plot(x_data, h_model, label='model', c='k')
#         ax.scatter(x_data, h_data, label='data', c=c_rms)
#         ax.legend(frameon=False)
#         ax.text(0.01, 0.98, r'$\Delta T_{rh}$=' + '{:.2f}'.format(y_data), fontsize=labelsize-4, ha='left', va='top',
#                 transform=ax.transAxes)  # label
#         if logx:
#             ax.set_xscale('log')
#         if logy:
#             ax.set_yscale('log')
#         if ylim is not None:
#             ax.set_ylim(ylim[0], ylim[1])  # for fair comparison
#         if xlim is not None:
#             ax.set_xlim(xlim)
#         # ax.set_ylabel(ylabel, fontsize=labelsize)
#         # ax.set_xlabel(xlabel, fontsize=labelsize)
#         # ax.set_title(title, fontsize=labelsize)
#
#     if save:
#         plot_save(fig, fname, fig_path=fig_path, fig_fmt=fig_fmt)
#     return fig, ax


#
# def fit_2d_log(x, y, h):
#     try:
#         x1 = np.log10(np.array(x))  # this should work for time-series of all x corresponding to h
#         y1 = np.log10(np.array(y))
#         h1 = np.log10(np.array(h))
#     except Exception as e:
#         print('h', h, type(h))
#         print('x', x, type(x))
#         print('y', y, type(y))
#         raise e
#
#     # define the data/predictors as the pre-set feature names
#     df = pd.DataFrame({'x1':x1, 'y1':y1})
#
#     # Put the target (housing value -- MEDV) in another DataFrame
#     target = pd.DataFrame({'h1':h1})
#
#     X = df
#     y = target['h1']
#     lm = linear_model.LinearRegression()
#     model = lm.fit(X, y)
#
#     coefs = lm.coef_
#     intercept = lm.intercept_
#
#     print('coefs', coefs)
#     print('intercept', intercept)
#
#     return coefs[0], coefs[1], 10**intercept
#
#     # def func(xdata_tuple, a, b, c):
#     #     (x, y) = xdata_tuple
#     #     g = a * x**b * y**c
#     #     return g.ravel()
#     # popt, pcov = curve_fit(func, (x1, x2), h.ravel())
