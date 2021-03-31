import sh_things as sh
import numpy as np
from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, \
    fig_fmt, regime_grid_td, load_grid, p_Earth, postprocess_kwargs

""" set dimensionalisation factors """
R_p = 6371
# d, dT, alpha = 600, 442, 4e-5 # Lees table 1-2: Ra=1e6
# d, dT, alpha = 2700, 3000, 2e-5  # Venus
d, dT, alpha = 2890, 3000, 3e-5  # Hoggard AGU Monograph
# d, dT, alpha = 2700, 3000, 3e-5  # test


""" """


""" get all time-averaged spectra and store """

regimes_use = ['chaotic']
for ii, eta in enumerate(eta_ls):  # across eta_ls
    cases_ii = ['Ra' + Ra + '-eta' + eta + e for Ra, e in zip(Ra_ls, end_grid[ii])]
    labels_ii = ['Ra=' + Ra for Ra in Ra_ls]
    for jj, Ra in enumerate(Ra_ls):
        if regime_grid_td[ii, jj] in regimes_use:
            case = cases_ii[jj]
            t1 = t1_grid[ii, jj]
            print('Calculating spectrum for', case)
            fig, ax = sh.dct_spectrum_avg(case, L_x=8,
                                          dim=False, R_p=d, d=d, dT=dT, alpha=alpha,
                                          t0=t1, x_res=1, t_res=1,
                                          test=False, data_path=data_path, fig_path=fig_path,
                                          check_norm=False,
                                          plot=True, load=False, dump=True, save=True, y0_guide=1e0,
                                          )
            print('    ...finished!')
