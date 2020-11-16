import sys
sys.path.insert(0, '/home/cmg76/Works/exo-top/')
# import numpy as np
# from exotop import aspect_postprocessing2 as asp
from exotop import aspect_scalings as sc
from setup_postprocessing import Ra_ls, eta_ls, t1, end, data_path, fig_path, c_rms, c_peak, fig_fmt, regime_grid

# plot summaries

# i_plot = list(range(len(eta_ls))
# for ii, eta in enumerate(eta_ls[i_plot]):  # across eta_ls
#     if ii in i_plot:
#         cases_ii = ['Ra' + Ra + '-eta' + eta + e for Ra, e in zip(Ra_ls, end[i_plot, ii])]
#         labels_ii = ['Ra=' + Ra for Ra in Ra_ls]
#         fig, ax = sc.subplots_cases(
#             cases_ii,
#             labels=labels_ii,
#             t1=t1[i_plot, ii], save=True, load='auto',
#             fname='all-eta' + eta, suptitle='$\Delta \eta$ = ' + eta,
#             includepdf=True, includeTz=True, show_sols=True,  # set False for faster summary with stats only
#             includegraphic=True, data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt, regime_grid=regime_grid[i_plot, ii],
#         )

# compare 64 and 129 resolution for Ra=3e7
# fig, ax = sc.case_subplots(
#     ['Ra3e7-eta1e5-wide', 'Ra3e7-eta1e5-wide-128'],
#         ['64', '128'], data_path=data_path, fig_path=fig_path,
#     t1=[0.25, 0.25], loadpickle=True, dumppickle=dump, fname='all-Ra3e7-res.png', suptitle='$\Delta \eta$=1e5, Ra=3e7',
#     includepd=True, # turn on once you know where steady state starts
#    )

i_plot = list(range(5,6))  # range(len(Ra_ls))
for ii, Ra in enumerate(Ra_ls):  # across Ra_ls
    if ii in i_plot:
        cases_ii = ['Ra' + Ra + '-eta' + eta + e for eta, e in zip(eta_ls, end.T[ii])]
        labels_ii = [r'$\Delta \eta$=' + eta for eta in eta_ls]
        fig, ax = sc.subplots_cases(
            cases_ii, labels=labels_ii, t1=t1.T[ii], save=True, load='auto',
            fname='all-Ra-' + Ra, suptitle='Ra = 3e8',
            includepdf=True, includeTz=True, show_sols=True,  # set False for faster summary with stats only
            includegraphic=True, data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt, regime_grid=regime_grid.T[i_plot, ii],
        )

print('Summary plots complete')
