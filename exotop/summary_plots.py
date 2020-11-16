import sys
sys.path.insert(0, '/home/cmg76/Works/exo-top/')
# import numpy as np
# from exotop import aspect_postprocessing2 as asp
from exotop import aspect_scalings as sc
from setup_postprocessing import Ra_ls, eta_ls, t1, end, data_path, fig_path, c_rms, c_peak, fig_fmt, regime_grid

# plot summaries

for ii, eta in enumerate(eta_ls):  # across eta_ls
    cases_ii = ['Ra' + Ra + '-eta' + eta + e for Ra, e in zip(Ra_ls, end[ii])]
    labels_ii = ['Ra=' + Ra for Ra in Ra_ls]
    fig, ax = sc.subplots_cases(
        cases_ii,
        labels=labels_ii,
        t1=t1[ii], save=True, load='auto',
        includeTz=True,
        fname='all-eta' + eta, suptitle='$\Delta \eta$=' + eta,
        includepdf=True,  # turn on once you know where steady state starts
        includegraphic=True, data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt, regime_grid=regime_grid[ii],
    )

# compare 64 and 129 resolution for Ra=3e7
# fig, ax = sc.case_subplots(
#     ['Ra3e7-eta1e5-wide', 'Ra3e7-eta1e5-wide-128'],
#         ['64', '128'], data_path=data_path, fig_path=fig_path,
#     t1=[0.25, 0.25], loadpickle=True, dumppickle=dump, fname='all-Ra3e7-res.png', suptitle='$\Delta \eta$=1e5, Ra=3e7',
#     includepd=True, # turn on once you know where steady state starts
#    )

for ii, Ra in enumerate(Ra_ls):  # across Ra_ls
    cases_ii = ['Ra' + Ra + '-eta' + eta + e for eta, e in zip(eta_ls, end.T[ii])]
    labels_ii = [r'$\Delta \eta$=' + eta for eta in eta_ls]
    fig, ax = sc.subplots_cases(
        cases_ii, labels=labels_ii, t1=t1.T[ii], save=True, load='auto',
        includeTz=True,
        fname='all-Ra-' + Ra, suptitle='Ra = 3e8',
        includepdf=True,  # turn on once you know where steady state starts
        includegraphic=True, data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt, regime_grid=regime_grid.T[ii],
    )

print('Summary plots complete')
