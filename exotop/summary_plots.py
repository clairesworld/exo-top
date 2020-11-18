import sys
sys.path.insert(0, '/home/cmg76/Works/exo-top/')
# import numpy as np
# from exotop import aspect_postprocessing2 as asp
from exotop import aspect_scalings as sc
from setup_postprocessing import Ra_ls, eta_ls, t1, end, data_path, fig_path, c_rms, c_peak, fig_fmt, regime_grid

# case = 'Ra3e8-eta1e5-wide'
# sc.pickle_remove_duplicate_row(case, suffix='_T', which='sol', data_path=data_path)

## plot summaries across delta eta

# i_plot = list(range(len(eta_ls)))
# for ii, eta in enumerate(eta_ls):  # across eta_ls
#     if ii in i_plot:
#         cases_ii = ['Ra' + Ra + '-eta' + eta + e for Ra, e in zip(Ra_ls, end[ii])]
#         labels_ii = ['Ra=' + Ra for Ra in Ra_ls]
#         fig, ax = sc.subplots_cases(
#             cases_ii, labels=labels_ii, t1=t1[ii], save=True, load='auto',
#             fname='all-eta' + eta, suptitle='$\Delta \eta$ = ' + eta,
#             includepdf=True, includeTz=True, show_sols=True,  # set False for faster summary with stats only
#             includegraphic=True, data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt, regime_grid=regime_grid[ii],
#         )

# compare 64 and 129 resolution for Ra=3e7
# fig, ax = sc.case_subplots(
#     ['Ra3e7-eta1e5-wide', 'Ra3e7-eta1e5-wide-128'],
#         ['64', '128'], data_path=data_path, fig_path=fig_path,
#     t1=[0.25, 0.25], loadpickle=True, dumppickle=dump, fname='all-Ra3e7-res.png', suptitle='$\Delta \eta$=1e5, Ra=3e7',
#     includepd=True, # turn on once you know where steady state starts
#    )

## plot summaries across Ra

# i_plot = list(range(len(Ra_ls)))  # range(4,5)
# for ii, Ra in enumerate(Ra_ls):  # across Ra_ls
#     if ii in i_plot:
#         cases_ii = ['Ra' + Ra + '-eta' + eta + e for eta, e in zip(eta_ls, end.T[ii])]
#         labels_ii = [r'$\Delta \eta$=' + eta for eta in eta_ls]
#         fig, ax = sc.subplots_cases(
#             cases_ii, labels=labels_ii, t1=t1.T[ii], save=True, load='auto',
#             fname='all-Ra-' + Ra, suptitle='Ra = '+Ra,
#             includepdf=True, includeTz=True, show_sols=True,  # set False for faster summary with stats only
#             includegraphic=True, data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt, regime_grid=regime_grid.T[ii],
#         )


## look at individual case data

case = 'Ra1e8-eta1e7-wide'
# sc.print_solution_data(case, suffix='_T', keys=['sol', 'time', 'delta_rh'], data_path=data_path)
sc.print_solution_data(case, suffix='_h_all', keys=['time', 'h_rms', 'h_peak'], data_path=data_path)

print('Summary plots complete')
