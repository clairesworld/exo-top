from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    fig_fmt, regime_grid_td, postprocess_kwargs, regime_names_td, \
    load_grid    # noqa: E402
from postaspect import plt_aspect as plat  # noqa: E402
from postaspect import aspect_post as ap
# import pstats


# (re)process all

# ap.reprocess_all_at_sol(Ra_ls, eta_ls, psuffixes=['_T'],  redo=True,
#                         t1_grid=t1_grid, end_grid=end_grid, load_grid=load_grid, regime_grid=regime_grid_td,
#                         include_regimes=['steady', 'trans.'],
#                         data_path=data_path,
#                         postprocess_kwargs=postprocess_kwargs)

ap.reprocess_all_at_sol(Ra_ls, eta_ls, psuffixes=['_T', '_h', '_Nu'], regime_names=regime_names_td,
                          t1_grid=t1_grid, end_grid=end_grid, data_path=data_path, redo=False,
                          load_grid=load_grid, regime_grid=regime_grid_td, postprocess_kwargs=postprocess_kwargs)

# ## plot summaries across delta eta
#
# i_plot = list(range(len(eta_ls)))
# for ii, eta in enumerate(eta_ls):  # across eta_ls
#     if ii in i_plot:
#         cases_ii = ['Ra' + Ra + '-eta' + eta + e for Ra, e in zip(Ra_ls, end_grid[ii])]
#         labels_ii = ['Ra=' + Ra for Ra in Ra_ls]
#         plat.subplots_cases(
#             cases_ii, labels=labels_ii, t1=t1_grid[ii], save=True, load=load_grid[ii],
#             fname='all-eta' + eta, suptitle='$\Delta \eta$ = ' + eta, c_rms=c_rms, c_peak=c_peak,
#             includepdf=True, includeTz=True, show_sols=True,  # set False for faster summary with stats only
#             includegraphic=True, data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt, regime_grid=regime_grid_td[ii],
#             postprocess_kwargs=postprocess_kwargs,
#         )


## plot summaries across Ra

i_plot = list(range(len(Ra_ls)))  # range(4,5)
for ii, Ra in enumerate(Ra_ls):  # across Ra_ls
    if ii in i_plot:
        cases_ii = ['Ra' + Ra + '-eta' + eta + e for eta, e in zip(eta_ls, end_grid.T[ii])]
        labels_ii = [r'$\Delta \eta$=' + eta for eta in eta_ls]
        plat.subplots_cases(
            cases_ii, labels=labels_ii, t1=t1_grid.T[ii], save=True, load=True,
            fname='all-Ra' + Ra, suptitle='Ra = '+Ra, c_rms=c_rms, c_peak=c_peak,
            includepdf=True, includeTz=True, show_sols=True,  # set False for faster summary with stats only
            includegraphic=True, data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt, regime_grid=regime_grid_td.T[ii],
            postprocess_kwargs=postprocess_kwargs,
        )


# compare 64 and 129 resolution for Ra=3e7
# fig, ax = plat.case_subplots(
#     ['Ra3e7-eta1e5-wide', 'Ra3e7-eta1e5-wide-128'],
#         ['64', '128'], data_path=data_path, fig_path=fig_path,
#     t1=[0.25, 0.25], loadpickle=True, dumppickle=dump, fname='all-Ra3e7-res.png', suptitle='$\Delta \eta$=1e5, Ra=3e7',
#     includepd=True, # turn on once you know where steady state starts
#    )


## look at individual case data
#
# case = 'Ra3e8-eta1e5-wide'
# pro.print_solution_data(case, suffix='_T', keys=None, data_path=data_path)
# pro.print_solution_data(case, suffix='_h', keys=None, data_path=data_path)

print('Summary plots complete')
