import sys
sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    fig_fmt, regime_grid_td, postprocess_kwargs, regime_names_td, \
    load_grid    # noqa: E402
from exotop.postaspect import aspect_scalings as sc  # noqa: E402
# import pstats

# (re)process all

# sc.reprocess_all_average(Ra_ls, eta_ls, suffix='_h_mean', postprocess_fn=sc.h_timeaverage, redo=True,
#                          t1_grid=t1_grid, end_grid=end_grid, data_path=data_path, load_grid=load_grid,
#                          regime_grid=regime_grid_td, include_regimes=['steady', 'trans.'])
#
# sc.reprocess_all_at_sol(Ra_ls, eta_ls, psuffixes=['_T'], postprocess_functions=[sc.T_parameters_at_sol], redo=True,
#                         t1_grid=t1_grid, end_grid=end_grid, load_grid=load_grid, regime_grid=regime_grid_td,
#                         include_regimes=['steady', 'trans.'],
#                         data_path=data_path,
#                         postprocess_kwargs=postprocess_kwargs)

sc.reprocess_all_at_sol(Ra_ls, eta_ls, psuffixes=['_T', '_h', '_Nu'], regime_names=regime_names_td,
                        postprocess_functions=[sc.T_parameters_at_sol, sc.h_at_ts, sc.Nu_at_ts],
                        t1_grid=t1_grid, end_grid=end_grid, data_path=data_path, redo=False,
                        load_grid=load_grid, regime_grid=regime_grid_td, postprocess_kwargs=postprocess_kwargs)

# ## plot summaries across delta eta
#
# i_plot = list(range(len(eta_ls)))
# for ii, eta in enumerate(eta_ls):  # across eta_ls
#     if ii in i_plot:
#         cases_ii = ['Ra' + Ra + '-eta' + eta + e for Ra, e in zip(Ra_ls, end_grid[ii])]
#         labels_ii = ['Ra=' + Ra for Ra in Ra_ls]
#         sc.subplots_cases(
#             cases_ii, labels=labels_ii, t1=t1_grid[ii], save=True, load=load_grid[ii],
#             fname='all-eta' + eta, suptitle='$\Delta \eta$ = ' + eta, c_rms=c_rms, c_peak=c_peak,
#             includepdf=True, includeTz=True, show_sols=True,  # set False for faster summary with stats only
#             includegraphic=True, data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt, regime_grid=regime_grid_td[ii],
#             postprocess_kwargs=postprocess_kwargs,
#         )
#
## plot summaries across Ra

i_plot = list(range(len(Ra_ls)))  # range(4,5)
for ii, Ra in enumerate(Ra_ls):  # across Ra_ls
    if ii in i_plot:
        cases_ii = ['Ra' + Ra + '-eta' + eta + e for eta, e in zip(eta_ls, end_grid.T[ii])]
        labels_ii = [r'$\Delta \eta$=' + eta for eta in eta_ls]
        sc.subplots_cases(
            cases_ii, labels=labels_ii, t1=t1_grid.T[ii], save=True, load=True,
            fname='all-Ra' + Ra, suptitle='Ra = '+Ra, c_rms=c_rms, c_peak=c_peak,
            includepdf=True, includeTz=True, show_sols=True,  # set False for faster summary with stats only
            includegraphic=True, data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt, regime_grid=regime_grid_td.T[ii],
            postprocess_kwargs=postprocess_kwargs,
        )



# # plot evolutions for debugging T components
# sc.subplots_evol_at_sol(Ra_ls, eta_ls, regime_grid=regime_grid_td, save=True, t1_grid=t1_grid, load_grid=True,
#                         psuffixes=['_T'], postprocess_functions=[sc.T_parameters_at_sol], fig_path=fig_path,
#                         fname='evol', fig_fmt=fig_fmt, end_grid=end_grid, normtime=True, labelsize=14, xlabel=r'Time',
#                         keys=['y_L', 'T_l', 'T_i', 'dT_rh', 'delta_rh', 'h_components'],
#                         ylabels=['y_L', 'T_L', 'T_i', 'dT_rh', 'delta_rh', r'$\alpha \Delta T{_{rh} \delta_{rh}$'], title='', xlabelpad=8,
#                         ylabelpad=8, markers=None, markersize=24, cmap='magma', vmin=5, vmax=8.5,
#                         include_regimes=['chaotic'], data_path=data_path, postprocess_kwargs=postprocess_kwargs)

sc.subplots_evol_at_sol(Ra_ls[-1], eta_ls[1:], regime_grid=regime_grid_td[1:,-1], save=True, t1_grid=t1_grid[1:,-1], load_grid=True,
                        psuffixes=['_T'], postprocess_functions=[sc.T_parameters_at_sol], fig_path=fig_path,
                        fname='evol-Ra3e8', fig_fmt=fig_fmt, end_grid=end_grid[1:,-1], normtime=True, labelsize=14, xlabel=r'Time',
                        ylabels=None, keys=['y_L', 'T_l', 'T_i', 'dT_rh', 'delta_rh', 'h_components'], title='', xlabelpad=8,
                        ylabelpad=8, markers=None, markersize=24, cmap='magma', vmin=5, vmax=8.5,
                        data_path=data_path, regime_names=regime_names_td, postprocess_kwargs=postprocess_kwargs)

# sc.subplots_evol_at_sol(Ra_ls[-2:], eta_ls[1], regime_grid=regime_grid_td, save=True, t1_grid=t1_grid[1,-2:], load_grid=True,
#                         psuffixes=['_T'], postprocess_functions=[sc.T_parameters_at_sol], fig_path=fig_path,
#                         fname='evol-eta1e6', fig_fmt=fig_fmt, end_grid=end_grid[1,-2:], normtime=True, labelsize=14, xlabel=r'Time',
#                         ylabels=None, keys=['y_L', 'T_l', 'T_i', 'dT_rh', 'delta_rh', 'h_components'], title='', xlabelpad=8,
#                         ylabelpad=8, markers=None, markersize=24, cmap='magma', vmin=5, vmax=8.5,
#                         data_path=data_path, regime_names=regime_names_td, postprocess_kwargs=postprocess_kwargs)

# sc.subplots_evol_at_sol(Ra_ls[-2:], eta_ls[2], regime_grid=regime_grid_td, save=True, t1_grid=t1_grid[2,-2:], load_grid=True,
#                         psuffixes=['_T'], postprocess_functions=[sc.T_parameters_at_sol], fig_path=fig_path,
#                         fname='evol-eta1e7', fig_fmt=fig_fmt, end_grid=end_grid[2,-2:], normtime=True, labelsize=14, xlabel=r'Time',
#                         ylabels=None, keys=['y_L', 'T_l', 'T_i', 'dT_rh', 'delta_rh', 'h_components'], title='', xlabelpad=8,
#                         ylabelpad=8, markers=None, markersize=24, cmap='magma', vmin=5, vmax=8.5,
#                         data_path=data_path, regime_names=regime_names_td, postprocess_kwargs=postprocess_kwargs)
#
# sc.subplots_evol_at_sol(Ra_ls[-2:], eta_ls[3], regime_grid=regime_grid_td, save=True, t1_grid=t1_grid[3,-2:], load_grid=True,
#                         psuffixes=['_T'], postprocess_functions=[sc.T_parameters_at_sol], fig_path=fig_path,
#                         fname='evol-eta1e8', fig_fmt=fig_fmt, end_grid=end_grid[3,-2:], normtime=True, labelsize=14, xlabel=r'Time',
#                         ylabels=None, keys=['y_L', 'T_l', 'T_i', 'dT_rh', 'delta_rh', 'h_components'], title='', xlabelpad=8,
#                         ylabelpad=8, markers=None, markersize=24, cmap='magma', vmin=5, vmax=8.5,
#                         data_path=data_path, regime_names=regime_names_td, postprocess_kwargs=postprocess_kwargs)



# compare 64 and 129 resolution for Ra=3e7
# fig, ax = sc.case_subplots(
#     ['Ra3e7-eta1e5-wide', 'Ra3e7-eta1e5-wide-128'],
#         ['64', '128'], data_path=data_path, fig_path=fig_path,
#     t1=[0.25, 0.25], loadpickle=True, dumppickle=dump, fname='all-Ra3e7-res.png', suptitle='$\Delta \eta$=1e5, Ra=3e7',
#     includepd=True, # turn on once you know where steady state starts
#    )


# ## look at individual case data
# #
# # case = 'Ra3e8-eta1e5-wide'
# # sc.print_solution_data(case, suffix='_T', keys=None, data_path=data_path)
# # sc.print_solution_data(case, suffix='_h', keys=None, data_path=data_path)

# look closely at individual pickles
# jj = 3
# eta = eta_ls[jj]
# cases_ii = ['Ra' + Ra + '-eta' + eta + e for Ra, e in zip(Ra_ls, end_grid[jj])]
# cases_ii = ['Ra3e8-eta1e5-wide', 'Ra1e8-eta1e7-wide']
# t1_ii = [0.0590015, 0.3]
# for ii, case in enumerate(cases_ii):
#     print('\n\n', case)
#     sc.pickleio(case, '_T', [sc.T_parameters_at_sol], t1=t1_ii[ii], load=False, at_sol=True, data_path=data_path)
#     sc.print_solution_data(case, suffix='_T', keys=['h_components', 'sol'], data_path=data_path)
#     # sc.print_solution_data(case, suffix='_h', keys=['h_rms', 'sol'], data_path=data_path)
#     # sc.print_solution_data(case, suffix='_h_all', keys=['h_rms', 'time'], data_path=data_path)
# case = 'Ra3e8-eta1e6-wide'
# print('\n\n', case)
# sc.pickleio(case, '_h', [sc.h_at_ts], t1=0.0590015, load=False, at_sol=True, data_path=data_path)
# sc.print_solution_data(case, suffix='_h', keys=['h_rms', 'sol'], data_path=data_path)
# #


# sc.pickleio('Ra1e8-eta1e6-wide', suffix='_T', postprocess_functions=[sc.T_parameters_at_sol], t1=0.3,
#             load=False, at_sol=True, data_path=data_path, postprocess_kwargs=postprocess_kwargs,)
#
# sc.pickleio('Ra3e8-eta1e6-wide', suffix='_T', postprocess_functions=[sc.T_parameters_at_sol], t1=0.055083,
#             load=False, at_sol=True, data_path=data_path, postprocess_kwargs=postprocess_kwargs,)
#


print('Summary plots complete')
