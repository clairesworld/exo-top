from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    fig_fmt, regime_grid_td, postprocess_kwargs, regime_names_td, \
    load_grid  # noqa: E402
from postaspect import plt_aspect as plat  # noqa: E402
from postaspect import aspect_post as ap
# import pstats

"""(re)process all"""

# ap.reprocess_all_at_sol(Ra_ls, eta_ls[4], psuffixes=['_T', '_h', '_h_all', '_Nu'],  redo=True,
#                         t1_grid=t1_grid[4,:], end_grid=end_grid[4,:], load_grid=load_grid[4,:], regime_grid=regime_grid_td[4,:],
#                         include_regimes=['chaotic'],
#                         data_path=data_path, postprocess_kwargs=postprocess_kwargs)

# ap.reprocess_all_at_sol(Ra_ls, eta_ls, psuffixes=['_T', '_h', '_h_all', '_Nu'], regime_names=regime_names_td,
#                         t1_grid=t1_grid, end_grid=end_grid, data_path=data_path, redo=False,
#                         load_grid=load_grid, regime_grid=regime_grid_td, postprocess_kwargs=postprocess_kwargs)

"""plot summaries across delta eta and/or Ra """

# i_plot = [4]  # list(range(len(eta_ls)))
# for ii, eta in enumerate(eta_ls):  # across eta_ls
#     if ii in i_plot:
#         cases_ii = ['Ra' + Ra + '-eta' + eta + e for Ra, e in zip(Ra_ls, end_grid[ii])]
#         labels_ii = ['Ra=' + Ra for Ra in Ra_ls]
#         plat.subplots_cases(
#             cases_ii, labels=labels_ii, t1=t1_grid[ii], save=True, load=load_grid[ii],
#             fname='all-eta' + eta, suptitle='$\Delta \eta$ = ' + eta, c_rms=c_rms, c_peak=c_peak,
#             includepdf=True, includeTz=True, show_sols=True,  # set False for faster summary with stats only
#             includegraphic=True, data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt,
#             regime_grid=regime_grid_td[ii],
#             postprocess_kwargs=postprocess_kwargs,
#         )

# i_plot = range(4, 6)  # list(range(len(Ra_ls)))  # range(4,5)
# for ii, Ra in enumerate(Ra_ls):  # across Ra_ls
#     if ii in i_plot:
#         cases_ii = ['Ra' + Ra + '-eta' + eta + e for eta, e in zip(eta_ls, end_grid.T[ii])]
#         labels_ii = [r'$\Delta \eta$=' + eta for eta in eta_ls]
#         plat.subplots_cases(
#             cases_ii, labels=labels_ii, t1=t1_grid.T[ii], save=True, load=True,
#             fname='all-Ra' + Ra, suptitle='Ra = ' + Ra, c_rms=c_rms, c_peak=c_peak,
#             includepdf=True, includeTz=True, show_sols=True,  # set False for faster summary with stats only
#             includegraphic=True, data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt,
#             regime_grid=regime_grid_td.T[ii],
#             postprocess_kwargs=postprocess_kwargs,
#         )

""" 2D isoviscous benchmark """
plat.subplots_cases(
    ['Lees-Ra1e6-2D'], labels=['2D isoviscous benchmark'], t1=[0.6], save=True, load='auto',
    fname='Lees-benchmark', suptitle='Ra = 1e6', c_rms=c_rms, c_peak=c_peak, dt_xlim=(0.0, 0.1),
    includepdf=True, includeTz=False, show_sols=False,  # set False for faster summary with stats only
    includegraphic=True, data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt,
    postprocess_kwargs=postprocess_kwargs,
)

# compare 64 and 129 resolution for Ra=3e7
# fig, ax = plat.case_subplots(
#     ['Ra3e7-eta1e5-wide', 'Ra3e7-eta1e5-wide-128'],
#         ['64', '128'], data_path=data_path, fig_path=fig_path,
#     t1=[0.25, 0.25], loadpickle=True, dumppickle=dump, fname='all-Ra3e7-res.png', suptitle='$\Delta \eta$=1e5, Ra=3e7',
#     includepd=True, # turn on once you know where steady state starts
#    )

# # eta 1e9 --> this has been subsumed into above
# fig, ax = plat.subplots_cases(
#     ['Ra1e8-eta1e9-wide-ascii', 'Ra3e8-eta1e9-wide-ascii'],
#     labels=['Ra 1e8', 'Ra 3e8'], data_path=data_path, fig_path=fig_path,
#     t1=[0.085, 0.07], load='auto', fname='all-eta1e9', suptitle='$\Delta \eta$=1e9',
#     includepd=False, includeTz=False,  # turn on once you know where steady state starts
# )

"""look at individual case data"""

# case = 'Ra3e8-eta1e5-wide'
# pro.print_solution_data(case, suffix='_T', keys=None, data_path=data_path)
# pro.print_solution_data(case, suffix='_h', keys=None, data_path=data_path)

print('Summary plots complete')
