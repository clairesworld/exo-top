""" ASPECT runs: functions for bulk postprocessing, storing/loading data, extracting scaling relationships, etc. """

import numpy as np
import pandas as pd
import pickle as pkl
from scipy.optimize import curve_fit
from scipy import stats
from scipy import odr
import os
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
# from scipy.spatial import distance
# import statsmodels.api as sm
from postaspect import aspectdata as ad
from postaspect.setup_postprocessing import data_path_bullard
from useful_and_bespoke import colorize, iterable_not_string, not_iterable, find_nearest_idx, \
    not_string, minmaxnorm, reduced_chisq, mahalanobis


def find_ts(case, t, dat=None, data_path=data_path_bullard, **kwargs):
    if dat is None:
        dat = ad.Aspect_Data(directory=data_path + 'output-' + case + '/', read_statistics=False,
                             read_parameters=False, **kwargs)

    dat.read_times(**kwargs)
    time = dat.stats_time
    ts = find_nearest_idx(time, t)
    return ts


def read_topo_stats(case, ts, data_path=data_path_bullard, fast=False, **kwargs):
    if fast:
        pass

    df = pd.read_csv(data_path + 'output-' + case + '/dynamic_topography_surface.' + '{:05}'.format(ts), header=None,
                     names=['x', 'y', 'h'], skiprows=1, index_col=False, delimiter=r"\s+", engine='python')
    return df['x'].to_numpy(), df['h'].to_numpy()


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
        # try:
        grids_new.append(reshape_one_input(grid, proper, default))
        # except ValueError:
        #     raise Exception('Problem reshaping input grid', grid, 'check indexing in first function call')

    # make sure Ra and eta are iterable the right way
    if not not_string(Ra_ls):  # if a string
        Ra_ls = [Ra_ls]
    if not not_string(eta_ls):
        eta_ls = [eta_ls]

    return Ra_ls, eta_ls, grids_new


def pickleio(case, suffix, t1=0, load='auto', dat_new=None,
             data_path=data_path_bullard, fend='.pkl', **kwargs):
    # do pickling strategy. only saving processed data for runs in quasi-steady state (past their t1 value)
    case_path = data_path + 'output-' + case + '/'
    fname = case + suffix + fend
    df = pd.DataFrame()
    dump_flag = False
    reprocess_flag = False
    t1_new = t1

    # retrieve eta
    eta_str = case[9:12]
    if 'col_vis' in kwargs:
        col_vis = kwargs.pop('col_vis')
    elif eta_str == '1e9':  # not sure what happened here but...
        col_vis = 23
    else:
        col_vis = 20
    # print('using col vis', col_vis)

    # auto-determine postprocessing based on pickle name
    if suffix == '_T':
        postprocess_functions = T_parameters_at_sol
        at_sol = True
    elif suffix == '_h':
        postprocess_functions = h_at_ts
        at_sol = True
    elif suffix == '_h_all':
        postprocess_functions = h_at_ts
        at_sol = False
    elif suffix == '_Nu':
        postprocess_functions = Nu_at_ts
        at_sol = True

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
                        dat_new = ad.Aspect_Data(directory=case_path,
                                                 read_statistics=False, read_parameters=False, **kwargs)
                        dat_new.read_times(**kwargs)
                    try:
                        if at_sol:
                            print('      Checking for new solutions...')
                            sol_f_old = df.sol.iat[-1]
                            sol_new = dat_new.read_stats_sol_files(col_vis=col_vis)
                            sol1_new = sol_new[np.argmax(sol_new > sol_f_old)]  # first solution after latest saved
                            t1_new = dat_new.find_time_at_sol(n=sol1_new, sol_files=sol_new, return_indices=False,
                                                              col_vis=col_vis)
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
                return None
            else:
                if not load:  # load is False so automatically calculate shit
                    print('    Declined to load', fname, ', processing afresh...')
                elif load == 'auto':
                    print('    File', fname, 'not found, processing...')
                else:
                    raise Exception('load value not understood:', load, type(load))
                reprocess_flag = True
                dat_new = ad.Aspect_Data(directory=case_path,
                                         read_statistics=False, read_parameters=False, **kwargs)
                dat_new.read_times(**kwargs)

            if reprocess_flag:
                if not hasattr(dat_new, 'stats_time'):
                    dat_new.read_times(**kwargs)
                dat_new.read_stats_heatflux(**kwargs)
                if at_sol:
                    if not hasattr(dat_new, 'sol_files'):
                        dat_new.read_stats_sol_files(col_vis=col_vis)
                    sol_new = dat_new.sol_files
                    df = process_at_solutions(case, postprocess_functions=postprocess_functions, dat=dat_new,
                                              t1=np.maximum(t1, t1_new),  # whichever comes later in time
                                              data_path=data_path, sol_files=sol_new, df_to_extend=df, col_vis=col_vis,
                                              **kwargs)
                    dump_flag = True  # always save if you did something

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


def pickleio_multi(case, psuffixes=None, t1=None, load=None, data_path=data_path_bullard,
                   postprocess_kwargs=None, **kwargs):
    dfs = []
    for ip, ps in enumerate(psuffixes):
        df1 = pickleio(case, suffix=ps, t1=t1,
                       load=load, data_path=data_path, postprocess_kwargs=postprocess_kwargs, **kwargs)
        dfs.append(df1)
    df = pd.concat(dfs, axis=1)
    df = df.loc[:, ~df.columns.duplicated()]
    return df


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
                index_key=None, test_run=False,
                **kwargs):
    case_path = data_path + 'output-' + case + '/'
    fname = case + suffix + fend
    df = pkl.load(open(case_path + 'pickle/' + fname, "rb"))  # open pickled file
    if index_key is not None:
        df.set_index(index_key, drop=False, inplace=False)
    print('df loaded and indexed\n', df.head(10))
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

    df2.reset_index(drop=False, inplace=False)
    print('dumping new df2\n', df2.head(10), '\n', df2.keys())
    if not test_run:
        pkl.dump(df2, open(case_path + 'pickle/' + fname, "wb"))
    else:
        print('  (just a test)')
    return df2


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


def peak_and_rms(h, trap=False):
    peak = np.max(h)
    if trap:
        rms = np.sqrt(trapzmean(h ** 2))
    else:
        rms = np.sqrt(np.mean(h ** 2))
    return peak, rms


def test_h_avg(case, data_path=data_path_bullard):
    dat = ad.Aspect_Data(directory=data_path + 'output-' + case + '/', read_statistics=True)
    ts = dat.stats_timestep[-1]  # final
    x, h = read_topo_stats(case, ts)
    h_mean = np.mean(h)
    h_peak = np.max(h)
    print('h mean raw', h_mean)
    print('h max raw', h_peak)

    # processed verj
    dic = h_at_ts(case, ts=ts)
    print('h max saved', dic['h_peak'])
    # --> h_avg = 0 --> h_peak - h_avg = h_peak


def time_averaged_profile(case, n0, nf, which='temperature', dat=None, data_path=data_path_bullard,
                          **kwargs):
    " get time-averaged profile for T, u, v etc. "
    profs_time = []
    if dat is None:
        dat = ad.Aspect_Data(directory=data_path + 'output-' + case + '/', read_statistics=False, **kwargs)
    for n in range(n0, nf + 1):
        print('Loading', which, 'profile at n =', n)
        if which == 'temperature':
            x, y, _, A = dat.read_temperature(n, **kwargs)
        elif which == 'velocity':
            x, y, _, _, _, _, A = dat.read_velocity(n, **kwargs)
        else:
            raise Exception('Profile type not implemented')
        prof = ad.horizontal_mean(A, x)
        profs_time.append(prof)
    profs_time = np.array(profs_time)
    profs_mean = np.mean(profs_time, axis=0)
    return profs_mean, y


def time_averaged_profile_from_df(df, col):
    nsols = len(df)
    try:
        srs = df[col].to_numpy()
    except KeyError:
        print(df.head(5))
    y = df.y.to_numpy()[0]
    profs = np.zeros((nsols, len(srs[-1])))
    for ii in range(nsols):
        try:
            profs[ii, :] = srs[ii]
        except:
            # this is probably because the mesh size changed lol
            raise Exception('error at ii', ii, 'with size', len(srs[ii]), '--> probably due to adaptive mesh size')
    av = np.mean(profs, axis=0)
    return av, y


def read_evol(case, col, dat=None, data_path=data_path_bullard, **kwargs):
    # return time, column i, (from statistics file)
    if dat is None:
        dat = ad.Aspect_Data(directory=data_path + 'output-' + case + '/', read_statistics=True, **kwargs)
    return dat.stats_time, eval('dat.stats_' + col)


def process_at_solutions(case, postprocess_functions, dat=None, t1=0, data_path=data_path_bullard,
                         df_to_extend=None, sol_files=None, postprocess_kwargs=None, **kwargs):
    if postprocess_kwargs is None:
        postprocess_kwargs = {}
    if df_to_extend is None:
        df_to_extend = pd.DataFrame()
    if dat is None:
        dat = ad.Aspect_Data(directory=data_path + 'output-' + case + '/', read_statistics=False,
                             read_parameters=False, **kwargs)
    try:
        time = dat.stats_time
    except AttributeError:
        dat.read_times(**kwargs)
        time = dat.stats_time
    i_time = np.argmax(time >= t1)  # index of first timestep to process

    if i_time > 0:
        if sol_files is None:
            try:
                sol_files = dat.sol_files
            except AttributeError:
                sol_files = dat.read_stats_sol_files(**kwargs)
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


def process_steadystate(case, postprocess_functions, dat=None, t1=0, data_path=data_path_bullard,
                        df_to_extend=None, postprocess_kwargs=None, **kwargs):
    if postprocess_kwargs is None:
        postprocess_kwargs = {}
    if df_to_extend is None:
        df_to_extend = pd.DataFrame()
    if dat is None:
        dat = ad.Aspect_Data(directory=data_path + 'output-' + case + '/', read_statistics=False,
                             read_parameters=False, **kwargs)
    try:
        time = dat.stats_time
    except AttributeError:
        dat.read_times(**kwargs)
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
        x, h = read_topo_stats(case, ts, **kwargs)
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


def h_timeaverage(case, ts0, tsf=1e50, **kwargs):
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
        raise Exception(case + 'ts0 does not catch topography')
    h_all = np.vstack(h_all)
    peak, rms = peak_and_rms(np.mean(h_all, axis=0))
    h_params['h_peak'] = [peak]
    h_params['h_rms'] = [rms]
    h_params['n'] = [np.shape(h_all)[0]]
    return pd.DataFrame.from_dict(h_params)


def Nu_at_ts(case, ts=None, dat=None, data_path=data_path_bullard, **kwargs):
    if dat is None:
        dat = ad.Aspect_Data(directory=data_path + 'output-' + case + '/', read_statistics=False,
                             read_parameters=True, **kwargs)
    N = dat.nusselt(k=1, **kwargs)
    return {'Nu': N[ts]}


def T_components_of_h(case, df=None, dat=None, psuffix='_T', data_path=data_path_bullard, update=False,
                      fend='.pkl', alpha_m=None, postprocess_kwargs=None, **kwargs):
    # calculate T components in h heuristic for all processed solutions, df can be dict
    if postprocess_kwargs is None:
        postprocess_kwargs = {}
    if 'alpha_m' in postprocess_kwargs:
        if (alpha_m is not None) and (alpha_m != postprocess_kwargs['alpha_m']):
            raise Exception('Error: competing alpha_m values passed to T_componenents_of_h. What the heck is going on!')
        alpha_m = postprocess_kwargs['alpha_m']
    if df is None:
        df = pickleio(case, suffix=psuffix,
                      dat_new=dat, data_path=data_path, fend=fend, **kwargs)

    try:
        h_components = alpha_m * np.array(df['dT_rh']) / np.array(df['dT_m']) * np.array(df['delta_rh']) / np.array(
            df['d_m'])
    except KeyError as e:
        print(df)
        raise e

    if update:
        df['h_components'] = h_components
        pkl.dump(df, open(data_path + 'output-' + case + '/pickle/' + case + psuffix + fend, 'wb'))

    return h_components


def T_parameters_at_sol(case, n, dat=None, T_av=None, uv_mag_av=None, y=None, data_path=data_path_bullard, **kwargs):
    if dat is None:
        dat = ad.Aspect_Data(directory=data_path + 'output-' + case + '/',
                             read_statistics=False, read_parameters=False, **kwargs)
    if y is None:
        dat.read_mesh(n)
        y = dat.y
    if uv_mag_av is None:
        x, y, z, u, v, _, uv_mag = dat.read_velocity(n, **kwargs)
        uv_mag_av = ad.horizontal_mean(uv_mag, x)
    if T_av is None:
        x, y, z, T = dat.read_temperature(n, **kwargs)
        T_av = ad.horizontal_mean(T, x)
    d_n = dat.T_components(n, T_av=T_av, uv_mag_av=uv_mag_av, y=y, data_path=data_path,
                           **kwargs)  # dict of components just at solution n
    d_n['h_components'] = T_components_of_h(case, df=d_n, dat=dat, data_path=data_path, **kwargs)

    return d_n


def parameter_percentiles(case=None, df=None, keys=None, sigma=2, **kwargs):
    # probability distribution for a number of parameters with df containing time evolutions
    # df can also be a dict...
    if sigma == 2:
        qs = [2.5, 50, 97.5]
    elif sigma == 1:
        qs = [16, 50, 84]
    else:
        raise Exception('Unrecognized sigma value')
    qdict = {}
    for key in keys:
        vals = df[key]  # .values
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

    return qdict


def fit_log(x, h, intercept=False, weights=None, slope=1, **kwargs):
    def coefficient(x, b):
        return x + b

    def coefficient2(x, b):
        global slope2
        return slope2 * x + b

    def f_wrapper_for_odr(beta, x):  # parameter order for odr
        return coefficient(x, *beta)

    try:
        x1 = np.log10(np.array(x))  # this should work for time-series of all x corresponding to h
        h1 = np.log10(np.array(h))
    except Exception as e:
        print('h', h, type(h))
        print('x', x, type(x))
        raise e
    if intercept:
        try:
            # check for and remove nans
            df = pd.DataFrame({'x': x1, 'h': h1})
            df.dropna(inplace=True)
            if slope == 1:
                parameters, cov = curve_fit(coefficient, df.x.to_numpy(), df.h.to_numpy())
            else:
                global slope2
                slope2 = slope
                parameters, cov = curve_fit(coefficient2, df.x.to_numpy(), df.h.to_numpy())
            slope = slope
            intercept = parameters[0]

            model = odr.odrpack.Model(f_wrapper_for_odr)
            data = odr.odrpack.Data(df.x.to_numpy(), df.h.to_numpy())
            myodr = odr.odrpack.ODR(data, model, beta0=parameters, maxit=0)
            myodr.set_job(fit_type=2)
            parameter_stats = myodr.run()
            df_e = len(x) - len(parameters)  # degrees of freedom, error
            cov_beta = parameter_stats.cov_beta  # parameter covariance matrix from ODR
            sd_beta = parameter_stats.sd_beta * parameter_stats.sd_beta
            ci = []
            t_df = stats.t.ppf(0.975, df_e)
            ci = []
            for i in range(len(parameters)):
                ci.append([parameters[i] - t_df * parameter_stats.sd_beta[i],
                           parameters[i] + t_df * parameter_stats.sd_beta[i]])

            tstat_beta = parameters / parameter_stats.sd_beta  # coeff t-statistics
            pstat_beta = (1.0 - stats.t.cdf(np.abs(tstat_beta), df_e)) * 2.0  # coef. p-values

            for i in range(len(parameters)):
                print('parameter:', parameters[i])
                print('   conf interval:', ci[i][0], ci[i][1])
                print('   tstat:', tstat_beta[i])
                print('   pstat:', pstat_beta[i])
                print()

        except ValueError as e:
            print('x', x1, np.shape(x1))
            print('h', h1, np.shape(h1))
            raise e
    elif weights is None:
        try:
            result = stats.linregress(x1, h1)
            slope = result.slope
            intercept = result.intercept
            print('\n  slope standard error', result.stderr,  # '  intercept standard error', result.intercept_stderr,
                  '  p value',
                  result.pvalue)  # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html
        except ValueError as e:
            print('x', np.shape(x1), 'h', np.shape(h1))
            raise e
    else:
        # use np polyfit for weighted least squares
        intercept, slope = np.polynomial.polynomial.polyfit(x1, h1, deg=1, w=weights)
    return slope, 10 ** intercept


#
# def fit_2log(x1, x2, h, **kwargs):
#     # https://stackoverflow.com/questions/35041266/scipy-odr-multiple-variable-regression
#     try:
#         df = pd.DataFrame({'x': np.log10(np.array(x1)), 'y': np.log10(np.array(x2)), 'h': np.log10(np.array(h))})
#     except Exception as e:
#         print('h', h, type(h))
#         print('x', x1, type(x1))
#         print('y', x2, type(x2))
#         raise e
#
#     def linfit(beta, x):
#         return beta[0] * x[0] + beta[1] * x[1] + beta[2]  # notice changed indices for x
#
#     print('x1:', np.log10(np.array(x1)))
#     print('x2:', np.log10(np.array(x2)))
#
#     x = np.row_stack((np.log10(np.array(x1)), np.log10(np.array(x2))))  # odr doesn't seem to work with column_stack
#
#     linmod = odr.Model(linfit)
#     data = odr.Data(x, np.log10(np.array(h)))
#     odrfit = odr.ODR(data, linmod, beta0=[1., 1., 1.])
#     parameter_stats = odrfit.run()
#     parameter_stats.pprint()
#     parameters = parameter_stats.beta
#
#     df_e = len(x1) - len(parameters)  # degrees of freedom, error
#     cov_beta = parameter_stats.cov_beta  # parameter covariance matrix from ODR
#     sd_beta = parameter_stats.sd_beta * parameter_stats.sd_beta
#     ci = []
#     t_df = stats.t.ppf(0.975, df_e)
#     ci = []
#     for i in range(len(parameters)):
#         ci.append([parameters[i] - t_df * parameter_stats.sd_beta[i],
#                    parameters[i] + t_df * parameter_stats.sd_beta[i]])
#
#     tstat_beta = parameters / parameter_stats.sd_beta  # coeff t-statistics
#     pstat_beta = (1.0 - stats.t.cdf(np.abs(tstat_beta), df_e)) * 2.0  # coef. p-values
#
#     for i in range(len(parameters)):
#         print('parameter:', parameters[i])
#         print('   conf interval:', ci[i][0], ci[i][1])
#         print('   tstat:', tstat_beta[i])
#         print('   pstat:', pstat_beta[i])
#         print()
#
#     intercept = 10 ** parameters[0]
#     slope1 = parameters[1]
#     slope2 = parameters[2]
#
#     return (slope1, slope2), intercept


def fit_SE(x, h, beta, err_x=1, err_h=1, xn=None, num=20):
    # only works for linear fit
    # convert to log
    try:
        logx1 = np.log10(np.array(x))  # this should work for time-series of all x corresponding to h
        logh = np.log10(np.array(h))
        err_logx1 = 0.434 * np.array(err_x) / np.array(x)
        err_logh = 0.434 * np.array(err_h) / np.array(h)
    except Exception as e:
        print('h', h, type(h))
        print('x', x, type(x))
        raise e

    # standard error (single parameter only) - where is this from??

    if xn is None:
        logxn1 = np.linspace(np.min(logx1), np.max(logx1), num)
    else:
        logxn1 = np.log10(xn)
    loghn = beta[0] + logxn1 * beta[1]

    sigma_logh = np.std(logh)
    xbar = np.mean(logx1)
    n = len(logh)
    SE = sigma_logh * np.sqrt(1 / n + (logxn1 - xbar) ** 2 / np.sum((logxn1 - xbar) ** 2))
    print('         -> average SE_y:', np.mean(SE))

    # unlog
    SE_unlog = 2.302585 * 10 ** loghn * SE

    return SE_unlog


def fit_logerror(x1, h, x2=None, err_x=1, err_h=1, ci=0.95, slope=True, **kwargs):
    def func_lin(beta, u):
        a, b = beta
        return a + u * b

    def func_lin2(beta, u):
        a, b, c = beta
        return a + u[0] * b + u[1] * c

    def func_lin0(beta, u):
        a = beta
        return a + u

    if len(np.shape(x1)) != 1:
        x1 = np.array(x1).reshape(np.shape(h))

    if not slope:
        beta0 = [np.log10(2)]
        func = func_lin0
        x = x1
    elif x2 is None:
        beta0 = [np.log10(0.1), -0.1]
        func = func_lin
        x = x1
    else:
        beta0 = [np.log10(0.1), -0.1, -0.1]
        func = func_lin2
        if len(np.shape(x2)) != 1:
            x2 = np.array(x2).reshape(np.shape(h))
        x = np.row_stack((x1, x2))  # odr doesn't seem to work with column_stack
        err_x = 1  # doesn't work with multidimensional weights
        err_h = 1

    # convert to log
    try:
        logx = np.log10(np.array(x))  # this should work for time-series of all x corresponding to h
        logh = np.log10(np.array(h))
        if err_x == 1:
            err_logx = 1
        else:
            err_logx = 0.434 * np.array(err_x) / np.array(x)
        if err_h == 1:
            err_logh = 1
        else:
            err_logh = 0.434 * np.array(err_h) / np.array(h)
    except Exception as e:
        print('h', h, type(h))
        print('x', x, type(x))
        raise e

    data = odr.RealData(logx, logh, sx=err_logx, sy=err_logh)
    model = odr.Model(func)
    try:
        odrfit = odr.ODR(data, model, beta0)
    except Exception as e:
        print('logx', logx)
        print('logh', logh)
        raise e
    output = odrfit.run()

    if len(output.beta) == 1:
        logcon = output.beta
        s_logcon = output.sd_beta
    elif len(output.beta) == 2:
        logcon, power = output.beta
        s_logcon, s_pow = output.sd_beta
    elif len(output.beta) == 3:
        logcon, power, power2 = output.beta
        s_logcon, s_pow, s_pow2 = output.sd_beta
    # s_logcon, s_pow = np.sqrt(np.diag(output.cov_beta))

    con = 10 ** logcon
    s_con = 2.302585 * 10 ** logcon * s_logcon

    # confidence intervals

    df_e = len(x1) - len(output.beta)  # degrees of freedom, error
    conf = []
    t_df = stats.t.ppf(ci, df_e)  # 0.975
    for i in range(len(output.beta)):
        conf.append([output.beta[i] - t_df * output.sd_beta[i],
                     output.beta[i] + t_df * output.sd_beta[i]])

    # chi sqr

    expected = 10 ** func(output.beta, logx)
    chisqr = np.sum(((h - expected) ** 2) / expected)
    MSE = chisqr / df_e

    print('\n')
    print('       -> ODR RESULTS')
    print('         -> reason for halting:', output.stopreason)
    if slope:
        print('         -> slope:', power, '+/-', s_pow, '   CI:', conf[1][0], conf[1][1])
    print('         -> intercept:', logcon, '+/-', s_logcon, '   CI:', conf[0][0], conf[0][1])
    print('         -> constant = 10^intercept:', con, '+/-', s_con)
    print('         -> chi sqr:', chisqr)
    print('         -> MSE:', MSE)

    return output.beta, output.sd_beta, chisqr, MSE


def fit_powererror(x1, h, x2=None, err_x=1, err_h=1, ci=0.95, **kwargs):
    def func_pow1(beta, u):
        a, b = beta
        return a * u ** b

    def func_pow2(beta, u):
        a, b = beta
        return a * u[0] ** (b * u[1])

    if len(np.shape(x1)) != 1:
        x1 = np.array(x1).reshape(np.shape(h))
    elif x2 is None:
        beta0 = [np.log10(0.1), -0.1]
        func = func_pow1
        x = x1
    else:
        beta0 = [np.log10(0.1), -0.1]
        func = func_pow2
        if len(np.shape(x2)) != 1:
            x2 = np.array(x2).reshape(np.shape(h))
        x = np.row_stack((x1, x2))  # odr doesn't seem to work with column_stack
        err_x = 1  # doesn't work with multidimensional weights
        err_h = 1

    data = odr.RealData(x, h, sx=err_x, sy=err_h)
    model = odr.Model(func)
    try:
        odrfit = odr.ODR(data, model, beta0)
    except Exception as e:
        print('x', x)
        print('h', h)
        raise e
    output = odrfit.run()

    if len(output.beta) == 2:
        con, power = output.beta
        s_con, s_pow = output.sd_beta
    else:
        print('beta length not implemented')

    # confidence intervals
    df_e = len(x1) - len(output.beta)  # degrees of freedom, error
    conf = []
    t_df = stats.t.ppf(ci, df_e)  # 0.975
    for i in range(len(output.beta)):
        conf.append([output.beta[i] - t_df * output.sd_beta[i],
                     output.beta[i] + t_df * output.sd_beta[i]])

    # chi sqr

    expected = func(output.beta, x)
    chisqr = np.sum(((h - expected) ** 2) / expected)
    MSE = chisqr / df_e

    print('\n')
    print('       -> ODR RESULTS')
    print('         -> reason for halting:', output.stopreason)
    print('         -> constant:', con, '+/-', s_con, '   CI:', conf[0][0], conf[0][1])
    print('         -> power:', power, '+/-', s_pow, '   CI:', conf[1][0], conf[1][1])
    print('         -> chi sqr:', chisqr)
    print('         -> MSE:', MSE)

    return output.beta, output.sd_beta, chisqr, MSE


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


def solomatov95(Ra=None, d_eta=None, df=None, case=None, dat=None,
                data_path=data_path_bullard, load='auto', **kwargs):
    if df is None:
        df = pickleio(case, suffix='_T', dat_new=dat, data_path=data_path,
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
        df = pickleio(case, suffix='_T', dat_new=dat, data_path=data_path,
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


def Ra_F_eff(delta_L=None, q_sfc=None, T_i=None, d_eta=None, T1=1, T0=0, Z=1):
    theta = np.log(d_eta)  # gamma for this delta eta
    gamma = theta / (np.array(T1) - np.array(T0))
    eta_i = np.exp(-gamma * np.array(T_i))
    return (Z - delta_L) ** 4 * q_sfc / eta_i


def regime_to_digital(ii=None, jj=None, regime_grid=None, regime_names=None, **kwargs):
    label = regime_grid[jj, ii]
    digi = np.nonzero(np.array(regime_names) == label)[0]
    if not list(digi):
        return np.nan
    else:
        return digi[0] + 1


def surf_mobility_at_sol(case=None, dat=None, n=None, data_path=data_path_bullard, **kwargs):
    if dat is None:
        dat = ad.Aspect_Data(directory=data_path + 'output-' + case + '/',
                             read_statistics=False, read_parameters=False, **kwargs)
    df_T = pickleio(case, '_T', dat_new=dat,
                    data_path=data_path, **kwargs)

    if not df_T.empty:
        df_sol = df_T.set_index('sol')
        if n is None:
            delta_0 = np.mean(df_sol['delta_0'].to_numpy())
            uv_mag_av = np.mean(df_sol['uv_mag_av'].to_numpy())
            if np.isnan(delta_0):
                print(df_sol.loc[pd.isna(df_sol["delta_L"]), :].index)
                print(df_sol.loc[pd.isna(df_sol["delta_rh"]), :].index)
                delta_0 = np.nanmean(df_sol['delta_0'].to_numpy())
            if np.isnan(uv_mag_av).any():
                print(case, df_sol['uv_mag_av'])
                uv_mag_av = np.nanmean(df_sol['uv_mag_av'].to_numpy())
        else:
            delta_0 = df_sol.loc[n, 'delta_0']
            uv_mag_av = df_sol.loc[n, 'uv_mag_av']
        S = dat.surface_mobility(n=n, delta_0=delta_0, uv_mag_av=uv_mag_av)
        return S
    else:
        return np.nan


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


def reprocess_all_at_sol(Ra_ls, eta_ls, psuffixes, t1_grid=None, end_grid=None,
                         data_path=data_path_bullard, redo=True, load_grid=None, regime_grid=None, include_regimes=None,
                         regime_names=None, check_t0=False, test_run=False, **kwargs):
    Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid) = reshape_inputs(Ra_ls, eta_ls, (
        t1_grid, load_grid, end_grid, regime_grid))
    if include_regimes is None:
        include_regimes = regime_names
    print('load_grid for reprocessing:', load_grid)
    for jj, eta_str in enumerate(eta_ls):
        cases, Ra_var = get_cases_list(Ra_ls, eta_str, end_grid[jj])
        for ii, case in enumerate(cases):
            t1_ii = t1_grid[jj][ii]
            if (t1_ii != 1) and (os.path.exists(data_path + 'output-' + case)):
                if regime_grid[jj][ii] in include_regimes:
                    print('(Re)processing', case)
                    if redo:
                        # for recalculating everything if you fucked up e.g.
                        load = False
                    elif load_grid is not None:
                        load = load_grid[jj][ii]
                    else:
                        load = 'auto'
                    for ip, suffix in enumerate(psuffixes):
                        df = pickleio(case, suffix=suffix, t1=t1_ii,
                                 data_path=data_path, load=load, **kwargs)
                        if check_t0:
                            try:
                                if 'index' in df.columns:
                                    dat = ad.Aspect_Data(directory=data_path + 'output-' + case + '/',
                                                         read_statistics=False, read_parameters=False, **kwargs)
                                    try:
                                        time = dat.stats_time
                                    except AttributeError:
                                        dat.read_times(**kwargs)
                                        time = dat.stats_time
                                    i_time = np.argmax(time >= t1_ii)  # index of first timestep to process
                                    ts_save = np.arange(i_time + 1, len(time))
                                    print('           ts range', ts_save, 'given t1', t1_ii, 'time[i_time]',
                                          time[i_time])
                                    idx_stored = df.index.values
                                    if idx_stored[0] < ts_save[0]:
                                        droppy = np.isin(idx_stored,
                                                         ts_save)  # these are what u want to keep but im attache to the name droppy
                                        print('           dropping idx (ts?)', idx_stored[~droppy])
                                        df = pickle_drop(case, suffix, keys=None, index=idx_stored[~droppy],
                                                         errors='raise', data_path=data_path,
                                                         test_run=test_run, **kwargs)
                            except AttributeError as e:
                                print(e)
                        print('\n', case, suffix)
                        print(df.head())
    print('>>>>>>>  done reprocessing!')


# def pickleio_average(case, postprocess_fn=None, t1=0, load=True, suffix='', data_path=data_path_bullard,
#                      fend='.pkl', **kwargs):
#     case_path = data_path + 'output-' + case + '/'
#     fname = case + suffix + fend
#     if not load:
#         dat = ad.Aspect_Data(directory=case_path,
#                              read_statistics=False, read_parameters=False, **kwargs)
#         dat.read_times(**kwargs)
#         time = dat.stats_time
#         i_time = np.argmax(time >= t1)  # index of first timestep to process
#         ts0 = i_time
#         tsf = len(time) - 1
#         df = postprocess_fn(case, ts0, tsf, **kwargs)
#         print('Processed', postprocess_fn, 'for time steps', ts0, 'to', tsf)
#         pkl.dump(df, open(case_path + 'pickle/' + fname, "wb"))
#     else:
#         df = pkl.load(open(case_path + 'pickle/' + fname, "rb"))
#     return df


# def reprocess_all_average(Ra_ls, eta_ls, t1_grid=None, end_grid=None,
#                           data_path=data_path_bullard, redo=True, load_grid=None, regime_grid=None,
#                           include_regimes=None, **kwargs):
#     Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid) = reshape_inputs(Ra_ls, eta_ls, (
#         t1_grid, load_grid, end_grid, regime_grid))
#
#     for jj, eta_str in enumerate(eta_ls):
#         cases, Ra_var = get_cases_list(Ra_ls, eta_str, end_grid[jj])
#         for ii, case in enumerate(cases):
#             t1_ii = t1_grid[jj][ii]
#             if (t1_ii != 1) and (os.path.exists(data_path + 'output-' + case)):
#                 if include_regimes is not None and (regime_grid[jj][ii] in include_regimes):
#                     print('Found', case)
#                     if redo:
#                         # for recalculating everything if you fucked up e.g.
#                         load = False
#                     elif load_grid is not None:
#                         load = load_grid[jj][ii]
#                     else:
#                         load = 'auto'
#                     df = pickleio_average(case, t1=t1_ii, load=load, data_path=data_path, **kwargs)
#                     print(df)


def fit_wrapper(x, h, yerr=1, xerr=1, n_fitted=2, fit_linear=True, **kwargs):
    if len(x) > 1:  # can only fit if at least 2 data
        slope = True
        x1 = x
        x2 = None
        if n_fitted == 3:
            x1 = x[0]
            x2 = x[1]
        elif n_fitted == 1:
            slope = False
        if fit_linear:
            beta, sd_beta, chisqr, MSE = fit_logerror(x1=x1, h=h, x2=x2, err_x=xerr, err_h=yerr, slope=slope, **kwargs)
            const = 10 ** beta[0]
            const_err = 2.302585 * 10 ** beta[0] * sd_beta[0]
        else:
            # fit power law directly
            print('x1 =', x1)
            print('x2 =', x2)
            print('h =', h)
            beta, sd_beta, chisqr, MSE = fit_powererror(x1=x1, h=h, x2=x2, err_x=xerr, err_h=yerr, **kwargs)
            const = beta[0]
            const_err = sd_beta[0]
        if len(beta) > 1:
            expon = [beta[1]]
            expon_err = [sd_beta[1]]
        else:
            expon = [1]
            expon_err = [0]
        if len(beta) > 2:
            expon.append(beta[2])
            expon_err.append(sd_beta[2])

    else:
        print('    Not enough points to fit')
        return [None] * 6

    return const, expon, const_err, expon_err, chisqr, MSE


def check_convergence(case, window=100, t1=0, fig_path='', plot=True, **kwargs):
    def savitzky_golay(y, window_size, order, deriv=0, rate=1):
        r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
        The Savitzky-Golay filter removes high frequency noise from data.
        It has the advantage of preserving the original shape and
        features of the signal better than other types of filtering
        approaches, such as moving averages techniques.
        Parameters
        ----------
        y : array_like, shape (N,)
            the values of the time history of the signal.
        window_size : int
            the length of the window. Must be an odd integer number.
        order : int
            the order of the polynomial used in the filtering.
            Must be less then `window_size` - 1.
        deriv: int
            the order of the derivative to compute (default = 0 means only smoothing)
        Returns
        -------
        ys : ndarray, shape (N)
            the smoothed signal (or it's n-th derivative).
        Notes
        -----
        The Savitzky-Golay is a type of low-pass filter, particularly
        suited for smoothing noisy data. The main idea behind this
        approach is to make for each point a least-square fit with a
        polynomial of high order over a odd-sized window centered at
        the point.
        Examples
        --------
        t = np.linspace(-4, 4, 500)
        y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
        ysg = savitzky_golay(y, window_size=31, order=4)
        import matplotlib.pyplot as plt
        plt.plot(t, y, label='Noisy signal')
        plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
        plt.plot(t, ysg, 'r', label='Filtered signal')
        plt.legend()
        plt.show()
        References
        ----------
        .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
           Data by Simplified Least Squares Procedures. Analytical
           Chemistry, 1964, 36 (8), pp 1627-1639.
        .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
           W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
           Cambridge University Press ISBN-13: 9780521880688
        """
        import numpy as np
        from math import factorial

        try:
            window_size = np.abs(np.int(window_size))
            order = np.abs(np.int(order))
        except ValueError as msg:
            raise ValueError("window_size and order have to be of type int")
        if window_size % 2 != 1 or window_size < 1:
            raise TypeError("window_size size must be a positive odd number")
        if window_size < order + 2:
            raise TypeError("window_size is too small for the polynomials order")
        order_range = range(order + 1)
        half_window = (window_size - 1) // 2
        # precompute coefficients
        b = np.mat([[k ** i for i in order_range] for k in range(-half_window, half_window + 1)])
        m = np.linalg.pinv(b).A[deriv] * rate ** deriv * factorial(deriv)
        # pad the signal at the extremes with
        # values taken from the signal itself
        firstvals = y[0] - np.abs(y[1:half_window + 1][::-1] - y[0])
        lastvals = y[-1] + np.abs(y[-half_window - 1:-1][::-1] - y[-1])
        y = np.concatenate((firstvals, y, lastvals))
        return np.convolve(m[::-1], y, mode='valid')


    t, q_top = read_evol(case, 'heatflux_top', **kwargs)
    _, q_bot = read_evol(case, 'heatflux_bottom', **kwargs)
    _, vel = read_evol(case, 'rms_velocity', **kwargs)

    try:
        q_diff = (-q_bot - q_top)/q_top
    except ValueError:
        q_top = np.insert(q_top, 0, 0)
        q_diff = (-q_bot - q_top) / q_top
    # q_diff_smooth = savitzky_golay(q_diff, 10, 2)

    # smooth using moving avgs
    try:
        df = pd.DataFrame({'t': t, 'vel':vel, 'qdiff':q_diff})
    except ValueError:
        print('t', np.shape(t), 'vel', np.shape(vel), 'q', np.shape(q_diff))
    # df.set_index('t', inplace=True)
    df['vel_rolling'] = df['vel'].rolling(window).mean()
    df['qdiff_rolling'] = df['qdiff'].rolling(window).mean()
    df.dropna()

    # for rms velocity, want to check percent difference and see where it is small and stable
    # note for q diff, just want to see where this gets to near 0
    df['vel_change'] = df['vel_rolling'].pct_change(periods=1)

    # df.drop(df[df.t < t1].index, inplace=True)

    # print(df.head())
    print(df.tail())

    if plot:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 1)
        for ax in axes:
            ax.set_xlabel('time')
            ax.axvline(x=t1, lw=1, c='k')
        axes[0].plot(df['t'], df['vel_change'], lw=0.5)
        axes[1].plot(df['t'], df['qdiff_rolling'], lw=0.5)

        axes[0].set_ylabel('percent change rms velocity')
        axes[1].set_ylabel('(q bottom - q top) / q top')
        fig.savefig(fig_path + 'converge_test_'+case+'.png', bbox_inches='tight')
    return None
