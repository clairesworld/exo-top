from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    fig_fmt, regime_grid_td, postprocess_kwargs, regime_names_td, \
    load_grid  # noqa: E402
from postaspect import plt_aspect as plat  # noqa: E402
from postaspect import aspect_post as ap
# import pstats

"""(re)process all"""

# ap.reprocess_all_at_sol(Ra_ls[4:], eta_ls[1:-1], psuffixes=['_h', '_h_all'],  redo=True,
#                         t1_grid=t1_grid[1:-1, 4:], end_grid=end_grid[1:-1, 4:], load_grid=load_grid[1:-1, 4:],
#                         regime_grid=regime_grid_td[1:-1, 4:],
#                         include_regimes=['chaotic'],
#                         data_path=data_path, postprocess_kwargs=postprocess_kwargs)

# ap.reprocess_all_at_sol(Ra_ls[4:], eta_ls[-1:], psuffixes=['_h', '_h_all'],  redo=True, col_vis=23,
#                         t1_grid=t1_grid[-1:, 4:], end_grid=end_grid[-1:, 4:], load_grid=load_grid[-1:, 4:],
#                         regime_grid=regime_grid_td[-1:, 4:],
#                         include_regimes=['chaotic'],
#                         data_path=data_path, postprocess_kwargs=postprocess_kwargs)  # eta 1e9 has diff col vis

# Ra 2e8 eta 1e9 even different?? that's col_vis=20??
# ap.reprocess_all_at_sol(Ra_ls[5], eta_ls[-1:], psuffixes=['_T', '_h', '_h_all', '_Nu'],  redo=True, col_vis=20,
#                         t1_grid=t1_grid[-1:, 5], end_grid=end_grid[-1:, 5], load_grid=load_grid[-1:, 5],
#                         regime_grid=regime_grid_td[-1:, 5],
#                         include_regimes=['chaotic'],
#                         data_path=data_path, postprocess_kwargs=postprocess_kwargs)  # eta 1e9 has diff col vis

# process latest from Ra 1e8
ap.reprocess_all_at_sol(Ra_ls[4:], eta_ls[1:], psuffixes=['_T', '_h', '_h_all', '_Nu'], redo=False,
                        t1_grid=t1_grid[1:, 4:], end_grid=end_grid[1:, 4:], load_grid=load_grid[1:, 4:],
                        regime_grid=regime_grid_td[1:, 4:],
                        check_t0=False, test_run=False, data_path=data_path, regime_names=regime_names_td,
                        postprocess_kwargs=postprocess_kwargs)

# Ra 2e8
ap.reprocess_all_at_sol(Ra_ls[-2], eta_ls, psuffixes=['_T', '_h', '_h_all', '_Nu'], redo=False,
                        t1_grid=t1_grid[:, -2], end_grid=end_grid[:, -2], load_grid=load_grid[:, -2],
                        regime_grid=regime_grid_td[:, -2],
                        data_path=data_path, regime_names=regime_names_td,
                        postprocess_kwargs=postprocess_kwargs)

"""plot summaries across delta eta and/or Ra """

# i_plot = [1, 2, 3, 4]  # list(range(len(eta_ls)))
# for ii, eta in enumerate(eta_ls):  # across eta_ls
#     if ii in i_plot:
#         cases_ii = ['Ra' + Ra + '-eta' + eta + e for Ra, e in zip(Ra_ls, end_grid[ii])]
#         labels_ii = ['Ra=' + Ra for Ra in Ra_ls]
#         plat.subplots_cases(
#             cases_ii, labels=labels_ii, t1=t1_grid[ii], save=True, load=True, #load_grid[ii],
#             fname='all-eta' + eta, suptitle='$\Delta \eta$ = ' + eta, c_rms=c_rms, c_peak=c_peak,
#             includepdf=True, includeTz=True, show_sols=True,  # set False for faster summary with stats only
#             includegraphic=True, data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt,
#             regime_grid=regime_grid_td[ii],
#             postprocess_kwargs=postprocess_kwargs,
#         )


i_plot = [5, 6]  # list(range(len(Ra_ls)))  # range(4,5)
for ii, Ra in enumerate(Ra_ls):  # across Ra_ls
    if ii in i_plot:
        cases_ii = ['Ra' + Ra + '-eta' + eta + e for eta, e in zip(eta_ls, end_grid.T[ii])]
        labels_ii = [r'$\Delta \eta$=' + eta for eta in eta_ls]
        plat.subplots_cases(
            cases_ii, labels=labels_ii, t1=t1_grid.T[ii], save=True, load=True,
            fname='all-Ra' + Ra, suptitle='Ra = ' + Ra, c_rms=c_rms, c_peak=c_peak,
            includepdf=True, includeTz=False, show_sols=True,  # set False for faster summary with stats only
            includegraphic=True, data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt,
            regime_grid=regime_grid_td.T[ii],
            postprocess_kwargs=postprocess_kwargs,
        )

# """ 2D isoviscous benchmark """
# plat.subplots_cases(
#     ['Lees-Ra1e6-2D'], labels=['2D isoviscous benchmark'], t1=[0.6], save=True, load='auto',
#     fname='Lees-benchmark', suptitle='Ra = 1e6', c_rms=c_rms, c_peak=c_peak, dt_xlim=(0.0, 0.1),
#     includepdf=True, includeTz=False, show_sols=False,  # set False for faster summary with stats only
#     includegraphic=True, data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt,
#     postprocess_kwargs=postprocess_kwargs,
# )

# compare 64 and 129 resolution for Ra=3e7
# fig, ax = plat.case_subplots(
#     ['Ra3e7-eta1e5-wide', 'Ra3e7-eta1e5-wide-128'],
#         ['64', '128'], data_path=data_path, fig_path=fig_path,
#     t1=[0.25, 0.25], loadpickle=True, dumppickle=dump, fname='all-Ra3e7-res.png', suptitle='$\Delta \eta$=1e5, Ra=3e7',
#     includepd=True, # turn on once you know where steady state starts
#    )

# eta 1e9 --> this has been subsumed into above tho
# fig, ax = plat.subplots_cases(
#     ['Ra1e8-eta1e9-wide-ascii', 'Ra2e8-eta1e9-wide-ascii', 'Ra3e8-eta1e9-wide-ascii'],
#     labels=['Ra 1e8', 'Ra 2e8', 'Ra 3e8'], data_path=data_path, fig_path=fig_path, col_vis=23,  # not sure why 23 for eta 1e9 only
#     t1=[0.15, 0.1], load=True, fname='all-eta1e9', suptitle='$\Delta \eta$=1e9',
#     includepd=True, includeTz=True, includegraphic=True,  # turn on once you know where steady state starts
# )

# Ra2e8 -->
# fig, ax = plat.subplots_cases(
#     ['Ra2e8-eta1e6-wide', 'Ra2e8-eta1e7-wide-ascii', 'Ra2e8-eta1e8-wide-ascii'],
#     labels=['$\Delta \eta$=1e6', '$\Delta \eta$=1e7', '$\Delta \eta$=1e8'], data_path=data_path, fig_path=fig_path,
#     t1=[0.06, 0.065, 0.06], load=True, fname='all-Ra2e8', suptitle='Ra=2e8', postprocess_kwargs=postprocess_kwargs,
#     includepd=True, includeTz=True, includegraphic=True,  # turn on once you know where steady state starts
# )

"""look at individual case data"""

# case = 'Ra2e8-eta1e7-wide-ascii'
# ap.print_solution_data(case, suffix='_T', keys=None, data_path=data_path)
# # ap.print_solution_data(case, suffix='_h', keys=None, data_path=data_path)

print('Summary plots complete')

