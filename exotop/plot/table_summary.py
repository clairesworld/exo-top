from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, load_grid, regime_grid_td, \
    data_path_bullard, fig_path_bullard, fig_fmt, postprocess_kwargs
from postaspect import plt_aspect as plat
from postaspect import aspect_post as pro
from useful_and_bespoke import not_iterable
import os
import numpy as np
import pandas as pd


def save_table(Ra, eta, fname, fig_path=fig_path_bullard, t1_grid=None, load_grid=None, end_grid=None, regime_grid=None,
               data_path=data_path_bullard, include_regimes=None, regime_names=None, cols=None,
               postprocess_kwargs=postprocess_kwargs, **kwargs):
    if include_regimes is None:
        include_regimes = regime_names
    if cols is None:
        cols = ('Ra_1', 'delta_eta', 'Ra_i_eff', 'delta_L', 'delta_rh', 'T_l', 'dT_rh', 'h_rms', 'h_peak')

    Ra, eta, (t1_grid, load_grid, end_grid, regime_grid) = pro.reshape_inputs(Ra, eta,
                                                                              (t1_grid, load_grid, end_grid,
                                                                               regime_grid))
    i = 0
    df_print = pd.DataFrame(columns=cols)
    # loop over cases
    for jj, etastr in enumerate(eta):
        cases, cases_var = pro.get_cases_list(Ra, etastr, end_grid[jj])
        for ii, case in enumerate(cases):
            if regime_grid[jj][ii] in include_regimes:
                t1_ii = t1_grid[jj][ii]
                load_ii = load_grid[jj][ii]
                df_T = pro.pickleio_multi(case, psuffixes=['_T'], t1=t1_ii, load=load_ii,
                                data_path=data_path, postprocess_kwargs=postprocess_kwargs, **kwargs)

                # load time-averages
                T_av, y = pro.time_averaged_profile_from_df(df_T, 'T_av')
                uv_mag_av, y = pro.time_averaged_profile_from_df(df_T, 'uv_mag_av')
                dic_av = pro.T_parameters_at_sol(case, n=None, T_av=T_av, uv_mag_av=uv_mag_av, y=y,
                                                 data_path=data_path, postprocess_kwargs=postprocess_kwargs, **kwargs)  # actually a dict
                # really hacky bit
                for k in ['T_av', 'uv_mag_av', 'y']:
                    dic_av.pop(k, None)
                    df_T = df_T.drop(k, axis=1)  # drop lists you don't need
                df_av = pd.DataFrame({key: value for (key, value) in dic_av.items()}, index=[0])
                df = df_T.mean(axis=0).to_frame().transpose()  # mean of other parameters
                df.set_index(pd.Series([0]))
                df.update(df_av)  # update with properly timefirst-averaged temperature params

                df_h = pro.pickleio_multi(case, psuffixes=['_h_all'], t1=t1_ii, load=load_ii,
                                        data_path=data_path, postprocess_kwargs=postprocess_kwargs, **kwargs)
                h_rms = df_h.h_rms.mean()
                h_peak = df_h.h_peak.mean()

                row = [Ra[ii], etastr]
                for col in cols[2:]:
                    if row not in ['h_peak', 'h_rms']:
                        row.append(df[col])
                row.extend([h_rms, h_peak])
                df_print.loc[i] = row
    df_print.to_csv(fig_path + fname)
    print(df_print.head())
    return df_print


save_table(Ra_ls, eta_ls, fname='test.csv', fig_path=fig_path_bullard, t1_grid=t1_grid, load_grid=load_grid,
           regime_grid=regime_grid_td, end_grid=end_grid, data_path=data_path_bullard, include_regimes=['chaotic'])




