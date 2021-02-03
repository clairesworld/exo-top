""" ASPECT runs: functions for bulk postprocessing, storing/loading data, extracting scaling relationships, etc. """

import numpy as np
import pandas as pd
import pickle as pkl
from scipy.optimize import curve_fit
from scipy import stats
from scipy import odr
import os
import matplotlib.pyplot as plt
from sklearn import linear_model
from scipy.interpolate import interp1d
from scipy.spatial import distance
# import statsmodels.api as sm
import sys

sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.postaspect import aspectdata as ad  # noqa: E402
from exotop.postaspect.plt_aspect import plot_save  # noqa: E402
from exotop.postaspect.setup_postprocessing import data_path_bullard  # noqa: E402
from exotop.useful_and_bespoke import colorize, iterable_not_string, not_iterable, \
    not_string, minmaxnorm, reduced_chisq, mahalanobis  # noqa: E402


def read_topo_stats(case, ts, data_path=data_path_bullard, **kwargs):
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


def pickleio(case, suffix, t1=0, load='auto', dat_new=None,
             data_path=data_path_bullard, fend='.pkl', **kwargs):
    # do pickling strategy. only saving processed data for runs in quasi-steady state (past their t1 value)
    case_path = data_path + 'output-' + case + '/'
    fname = case + suffix + fend
    df = pd.DataFrame()
    dump_flag = False
    reprocess_flag = False
    t1_new = t1

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
                dat_new = ad.Aspect_Data(directory=case_path,
                                         read_statistics=False, read_parameters=False, **kwargs)
                dat_new.read_times(**kwargs)

            if reprocess_flag:
                if not hasattr(dat_new, 'stats_time'):
                    dat_new.read_times(**kwargs)
                dat_new.read_stats_heatflux(**kwargs)
                if at_sol:
                    if not hasattr(dat_new, 'sol_files'):
                        dat_new.read_stats_sol_files()
                    sol_new = dat_new.sol_files
                    df = process_at_solutions(case, postprocess_functions=postprocess_functions, dat=dat_new,
                                              t1=np.maximum(t1, t1_new),  # whichever comes later in time
                                              data_path=data_path, sol_files=sol_new, df_to_extend=df, **kwargs)
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
    srs = df[col].to_numpy()
    y = df.y.to_numpy()[0]
    profs = np.zeros((nsols, len(srs[0])))
    for ii in range(nsols):
        profs[ii, :] = srs[ii]
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
        return -2*x + b

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
                popt, pcov = curve_fit(coefficient, df.x.to_numpy(), df.h.to_numpy())
            elif slope == -2:
                popt, pcov = curve_fit(coefficient2, df.x.to_numpy(), df.h.to_numpy())
            slope = slope
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
        df = pd.DataFrame({'x': np.log10(np.array(x)), 'y': np.log10(np.array(y)), 'h': np.log10(np.array(h))})
    except Exception as e:
        # print('h', h, type(h))
        # print('x', x, type(x))
        # print('y', y, type(y))
        raise e

    X = df[['x', 'y']]
    Y = df['h']

    # with sklearn
    regr = linear_model.LinearRegression()
    regr.fit(X, Y)

    # print('Intercept: \n', regr.intercept_)
    # print('Coefficients: \n', regr.coef_)

    return regr.coef_, 10 ** regr.intercept_


def fit_logerror(x, h, err_x, err_h, beta0=[0.1, -0.15], sigma=2, plot=True, **kwargs):
    print('x', x)
    print('h', h)
    print('err_x', err_x)
    print('err_h', err_h)
    try:
        logx = np.log10(np.array(x))  # this should work for time-series of all x corresponding to h
        logh = np.log10(np.array(h))
        logerr_x = np.log10(np.array(err_x))
        logerr_h = np.log10(np.array(err_h))
    except Exception as e:
        print('h', h, type(h))
        print('x', x, type(x))
        raise e

    def func_2param(p, u):
        a, b = p
        return a * u**b

    def func_lin2param(p, u):
        a, b = p
        return a + u*b

    def func_3param(p, u):
        a, b, c = p
        return a * u[0]**b * u[1]**c

    # Model object
    model = odr.Model(func_2param)
    # model = odr.Model(func_lin2param)

    # Create a RealData object
    data = odr.RealData(x, h, sx=err_x, sy=err_h)
    # data = odr.RealData(logx, logh, sx=logerr_x, sy=logerr_h)

    # Set up ODR with the model and data.
    odr_ = odr.ODR(data, model, beta0=beta0)

    # Run the regression.
    out = odr_.run()

    # print fit parameters and 1-sigma estimates
    popt = out.beta
    perr = out.sd_beta
    # popt[0] = 10**popt[0]  # unlog
    # perr[0] = 10**perr[0]

    # prepare confidence level curves
    nstd = sigma  # to draw 2-sigma intervals e.g.
    popt_up = popt + nstd * perr
    popt_dw = popt - nstd * perr

    # print('\nfit parameter', sigma, 'sigma error')
    # print('———————————–')
    # for i in range(len(popt)):
    #     print(str(popt[i]) + ' +- ' + str(perr[i]))
    #     print('min:', popt_dw[i])
    #     print('max:', popt_up[i])
    # print('y min @ x min', func_2param(popt_dw, x[0]))
    # print('y max @ x max', func_2param(popt_up, x[-1]))
    out.pprint()

    if plot:
        # plot
        x_fit = np.linspace(min(x), max(x), 100)
        fit = func_2param(popt, x_fit)
        fit_up = func_2param(popt_up, x_fit)
        fit_dw = func_2param(popt_dw, x_fit)
        fig, ax = plt.subplots(1)
        ax.errorbar(x, h, yerr=err_h, xerr=err_x, hold=True, ecolor='k', fmt ='none', label ='data')
        ax.set_xlabel('x', fontsize=18)
        ax.set_ylabel('h', fontsize=18)
        ax.set_title('fit with error on both axes', fontsize=18)
        plt.plot(x_fit, fit, 'r', lw=2, label='best fit curve')
        ax.fill_between(x_fit, fit_up, fit_dw, alpha=.25, label=str(sigma)+'-sigma interval')
        ax.legend(loc='lower right', fontsize=12)
        # ax.set_xscale('log')
        # ax.set_yscale('log')
        plot_save(fig, 'fit_test', **kwargs)
    return (*popt, *perr)


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
    if n is None:
        n = dat.final_step()

    df_T = pickleio(case, '_T', dat_new=dat,
                    data_path=data_path, **kwargs)
    if not df_T.empty:
        df_sol = df_T.set_index('sol')
        S = dat.surface_mobility(n=n, delta_0=df_sol.loc[n, 'delta_0'], delta_rh=df_sol.loc[n, 'delta_rh'],
                                 delta_l=df_sol.loc[n, 'delta_L'])
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
                         regime_names=None, **kwargs):
    Ra_ls, eta_ls, (t1_grid, load_grid, end_grid, regime_grid) = reshape_inputs(Ra_ls, eta_ls, (
        t1_grid, load_grid, end_grid, regime_grid))
    if include_regimes is None:
        include_regimes = regime_names

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
                        pickleio(case, suffix=suffix, t1=t1_ii,
                                 data_path=data_path, load=load, **kwargs)


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


def fit_cases_on_plot(yx_all, ax, yerr=None, xerr=None, legend=True, showallscatter=False, n_fitted=2, c_list=None,
                      c='xkcd:periwinkle', sigma=1, legsize=8, lw=1, legloc='lower left', showchisq=False, **kwargs):
    fiterror = (yerr is not None) and (xerr is not None)
    x = [a[1] for a in yx_all]
    y = [a[0] for a in yx_all]
    try:
        weights = [len(a[0]) for a in yx_all]
    except TypeError:
        weights = [1] * len(x)
    if np.array(x[0]).ndim > 0 and np.array(y[0]).ndim > 0:
        flatx = [item for sublist in x for item in sublist]
        flaty = [item for sublist in y for item in sublist]
    else:
        flatx, flaty = x, y
    if len(x) > 1:  # can only fit if at least 2 data

        if n_fitted > 2:  # fit to 3 parameter power law
            flatx0 = [a[0] for a in flatx]
            flatx1 = [a[1] for a in flatx]
            x0prime = np.linspace(np.min(flatx0), np.max(flatx0), num=len(flatx0))
            expon, const = fit_2log(x=flatx0, y=flatx1, h=flaty)
            z_vec = np.unique(flatx1)
            if c_list is None:
                c_list = colorize(np.log10(z_vec), cmap='winter')[0]
            for ind, z in enumerate(z_vec):
                hprime = const * x0prime ** expon[0] * z ** expon[1]
                h2, = ax.plot(x0prime, hprime, c=c_list[ind], ls='--', lw=lw, zorder=100, label='dum')

        else:
            xprime = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1])
            # xprime = np.linspace(np.min(flatx), np.max(flatx), num=len(flatx))
            if fiterror:
                const, expon, const_err, expon_err = fit_logerror(flatx, flaty, xerr, yerr, sigma=sigma, **kwargs)
                const_up = const + sigma * const_err
                const_dw = const - sigma * const_err
                expon_up = expon + sigma * expon_err
                expon_dw = expon - sigma * expon_err
                hprime_up = const_up * xprime ** expon_up
                hprime_dw = const_dw * xprime ** expon_dw
                hprime_dw[hprime_dw < 0] = 0  # cannot be negative
                ax.fill_between(xprime, hprime_up, hprime_dw, alpha=.25#, label=str(sigma) + '-sigma interval'
                                )
            else:
                expon, const = fit_log(flatx, flaty, weights=weights, **kwargs)
            hprime = const * xprime ** expon
            h3, = ax.plot(xprime, hprime, c=c, ls='--', lw=lw, zorder=100, label='dum')

        print('fit: {:.2e} x^{:.3f}'.format(const, expon))

        if legend:
            handles, labels = ax.get_legend_handles_labels()
            # handles.append(h3)
            # labels.append('{:.2e} x^{:.3f}'.format(const, expon))
            if n_fitted == 2:
                if fiterror:
                    newlabel = '({:.2e} +- {:.2e}) x^({:.3f} +- {:.3f})'.format(const, const_err, expon, expon_err)
                else:
                    newlabel = '{:.2e} x^{:.3f}'.format(const, expon)
            elif n_fitted == 3:
                newlabel = '{:.3e} x0^{:.3f} x1^{:.3f}'.format(const, expon[0], expon[1])
            else:
                raise Exception('Legend labelling for this n fitted parameters not implemented')
            if showchisq:
                chisq = reduced_chisq(O_y=np.log10(flaty), C_y=np.log10(hprime), n_fitted=n_fitted, **kwargs)
                newlabel = newlabel + r'; $\chi^2_\nu$ = ' + '{:.4f}'.format(chisq)
            try:
                labels[-1] = newlabel
            except IndexError:
                labels = newlabel
            leg = ax.legend(fontsize=legsize, handles=handles, labels=labels, loc=legloc)
            ax.add_artist(leg)
    else:
        print('    Not enough points to fit')
    if showallscatter:
        ax.scatter(flatx, flaty, c=c, alpha=0.05, s=10)
    return ax
