from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    fig_fmt, regime_grid_td, postprocess_kwargs, regime_names_td, \
    load_grid  # noqa: E402
from postaspect import plt_aspect as plat  # noqa: E402

# import pstats
cmap = 'cool'

# plat.subplots_hist(Ra_ls[-1], eta_ls[1:], regime_grid=regime_grid_td[1:, -1], save=True, t1_grid=t1_grid[1:, -1], load_grid=True,
#                    psuffixes=['_T', '_h_all'], fig_path=fig_path,
#                    fname='hist-Ra3e8', fig_fmt=fig_fmt, end_grid=end_grid[1:,-1], labelsize=14,
#                    keys=['y_L', 'T_l', 'T_i', 'dT_rh', 'delta_rh', 'h_rms'], title='Ra 3e8', xlabelpad=8,
#                    cmap='magma', vmin=5, vmax=8.5, nbins=10,
#                    data_path=data_path, regime_names=regime_names_td, postprocess_kwargs=postprocess_kwargs)

# plat.subplots_hist(Ra_ls[-2], eta_ls[1:], regime_grid=regime_grid_td[1:, -2], save=True, t1_grid=t1_grid[1:, -2], load_grid=True,
#                    psuffixes=['_T', '_h'], fig_path=fig_path,
#                    fname='hist-Ra1e8', fig_fmt=fig_fmt, end_grid=end_grid[1:,-2], labelsize=14,
#                    keys=['y_L', 'T_l', 'T_i', 'dT_rh', 'delta_rh', 'h_rms'], title='Ra 1e8', xlabelpad=8,
#                    cmap='magma', vmin=5, vmax=8.5, nbins=10,
#                    data_path=data_path, regime_names=regime_names_td, postprocess_kwargs=postprocess_kwargs)


# plat.subplots_hist(Ra_ls[-2:], eta_ls[-1], regime_grid=regime_grid_td[-1, -2:],  t1_grid=t1_grid[-1, -2:],
#                    load_grid=True, end_grid=end_grid[-1, -2:],
#                    psuffixes=['_T', '_h', '_Nu'], fig_path=fig_path, save=True,
#                    fname='hist-eta1e8', fig_fmt=fig_fmt, labelsize=14,
#                    keys=['y_L', 'T_l', 'dT_rh', 'delta_rh', 'h_rms', 'Nu'], title='eta 1e8', xlabelpad=8,
#                    cmap=cmap, nbins=10, colour_by='Ra',
#                    data_path=data_path, regime_names=regime_names_td, postprocess_kwargs=postprocess_kwargs)
#
# plat.subplots_hist(Ra_ls[-2:], eta_ls[-2], regime_grid=regime_grid_td[-2, -2:],  t1_grid=t1_grid[-2, -2:],
#                    load_grid=True, end_grid=end_grid[-2, -2:],
#                    psuffixes=['_T', '_h', '_Nu'], fig_path=fig_path, save=True,
#                    fname='hist-eta1e7', fig_fmt=fig_fmt, labelsize=14,
#                    keys=['y_L', 'T_l', 'dT_rh', 'delta_rh', 'h_rms', 'Nu'], title='eta 1e7', xlabelpad=8,
#                    cmap=cmap, nbins=10, colour_by='Ra',
#                    data_path=data_path, regime_names=regime_names_td, postprocess_kwargs=postprocess_kwargs)
#
# plat.subplots_hist(Ra_ls[-2:], eta_ls[-3], regime_grid=regime_grid_td[-3, -2:],  t1_grid=t1_grid[-3, -2:],
#                    load_grid=True, end_grid=end_grid[-3, -2:],
#                    psuffixes=['_T', '_h', '_Nu'], fig_path=fig_path, save=True,
#                    fname='hist-eta1e6', fig_fmt=fig_fmt, labelsize=14,
#                    keys=['y_L', 'T_l', 'dT_rh', 'delta_rh', 'h_rms', 'Nu'], title='eta 1e6', xlabelpad=8,
#                    cmap=cmap, nbins=10, colour_by='Ra',
#                    data_path=data_path, regime_names=regime_names_td, postprocess_kwargs=postprocess_kwargs)


#
#
# # # plot evolutions for debugging T components

zoom = False
normtime = True
percent_change = False
fig, axes = plat.subplots_evol_at_sol(Ra_ls[-2:], eta_ls, regime_grid=regime_grid_td[:,-2:], save=False,
                                      t1_grid=t1_grid[:,-2:], load_grid=True, psuffixes=['_T', '_h_all'], fig_path=fig_path,
                                      fig_fmt=fig_fmt, end_grid=end_grid[:,-2:], normtime=normtime, alpha=0.4,
                                      labelsize=14, xlabel=r'Time', percent_change=percent_change, zoom=zoom,
                                      keys=['T_l', 'delta_rh', 'h_rms'], title='', legend=False,
                                      xlabelpad=8, ylabelpad=8, markers=None, markersize=14, colour_by=None,
                                      data_path=data_path, regime_names=regime_names_td,
                                      postprocess_kwargs=postprocess_kwargs)

fig, axes = plat.subplots_evol_at_sol(Ra_ls[4], eta_ls[1], regime_grid=regime_grid_td[1, 4], save=False,
                                      t1_grid=t1_grid[1, 4], load_grid=True, psuffixes=['_T', '_h_all'], fig_path=fig_path,
                                      fname='evol_bi', fig_fmt=fig_fmt, end_grid=end_grid[1, 4], normtime=normtime,
                                      labelsize=14, xlabel=r'Time', percent_change=percent_change, alpha=1, zoom=zoom,
                                      keys=['T_l', 'delta_rh', 'h_rms'], title='', legend=False,
                                      xlabelpad=8, ylabelpad=8, markers=['X'], markersize=24, c_const='r', colour_by=None,
                                      data_path=data_path, regime_names=regime_names_td, fig=fig, axes=axes,
                                      postprocess_kwargs=postprocess_kwargs)

fig, axes = plat.subplots_evol_at_sol(Ra_ls[5], eta_ls[3], regime_grid=regime_grid_td[3, 5], save=True,
                                      t1_grid=t1_grid[3, 5], load_grid=True, psuffixes=['_T', '_h_all'], fig_path=fig_path,
                                      fname='evol_bi', fig_fmt=fig_fmt, end_grid=end_grid[3, 5], normtime=normtime,
                                      labelsize=14, xlabel=r'Time', percent_change=percent_change, alpha=1, zoom=zoom,
                                      keys=['T_l', 'delta_rh', 'h_rms'], title='', legend=False,
                                      xlabelpad=8, ylabelpad=8, markers=['X'], markersize=24, c_const='r', colour_by=None,
                                      data_path=data_path, regime_names=regime_names_td, fig=fig, axes=axes,
                                      postprocess_kwargs=postprocess_kwargs)

# plat.subplots_evol_at_sol(Ra_ls[-2:], eta_ls[1:], regime_grid=regime_grid_td[1:, -2:], save=True, t1_grid=t1_grid[1:, -2:],
#                           load_grid=True,
#                         psuffixes=['_T'], fig_path=fig_path,
#                         fname='evol', fig_fmt=fig_fmt, end_grid=end_grid[1:, -2:], normtime=True, labelsize=14, xlabel=r'Time',
#                         keys=['y_L', 'T_l', 'dT_rh', 'delta_rh', 'h_components'],
#                         ylabels=['y_L', 'T_L', 'dT_rh', 'delta_rh', r'$\alpha \Delta T{_{rh} \delta_{rh}$'], title='', xlabelpad=8,
#                         ylabelpad=8, markers=None, markersize=24, cmap=cmap, colour_by='Ra',
#                         data_path=data_path, regime_names=regime_names_td, postprocess_kwargs=postprocess_kwargs)

# plat.subplots_evol_at_sol(Ra_ls[-2], eta_ls[1:], regime_grid=regime_grid_td[1:, -2], save=True, t1_grid=t1_grid[1:, -2], load_grid=True,
#                           psuffixes=['_T'], fig_path=fig_path,
#                           fname='evol-Ra1e8', fig_fmt=fig_fmt, end_grid=end_grid[1:,-1], normtime=True, labelsize=14, xlabel=r'Time',
#                           ylabels=None, keys=['y_L', 'T_l', 'dT_rh', 'delta_rh', 'h_components'], title='', xlabelpad=8,
#                           ylabelpad=8, markers=None, markersize=24, cmap=cmap, colour_by='eta',
#                           data_path=data_path, regime_names=regime_names_td, postprocess_kwargs=postprocess_kwargs)
#
#
# plat.subplots_evol_at_sol(Ra_ls[-1], eta_ls[1:], regime_grid=regime_grid_td[1:, -1], save=True, t1_grid=t1_grid[1:, -1], load_grid=True,
#                           psuffixes=['_T'], fig_path=fig_path,
#                           fname='evol-Ra3e8', fig_fmt=fig_fmt, end_grid=end_grid[1:,-1], normtime=True, labelsize=14, xlabel=r'Time',
#                           ylabels=None, keys=['y_L', 'T_l', 'dT_rh', 'delta_rh', 'h_components'], title='', xlabelpad=8,
#                           ylabelpad=8, markers=None, markersize=24, cmap=cmap, colour_by='eta',
#                           data_path=data_path, regime_names=regime_names_td, postprocess_kwargs=postprocess_kwargs)
#
# plat.subplots_evol_at_sol(Ra_ls[-2:], eta_ls[1], regime_grid=regime_grid_td[1, -2:], save=True, t1_grid=t1_grid[1,-2:],
#                         psuffixes=['_T'], fig_path=fig_path, end_grid=end_grid[1,-2:],  load_grid=True,
#                         fname='evol-eta1e6', fig_fmt=fig_fmt, normtime=True, labelsize=14, xlabel=r'Time',
#                         ylabels=None, keys=['y_L', 'T_l', 'T_i', 'dT_rh', 'delta_rh', 'h_components'], title='', xlabelpad=8,
#                         ylabelpad=8, markers=None, markersize=24, cmap=cmap, colour_by='Ra',
#                         data_path=data_path, regime_names=regime_names_td, postprocess_kwargs=postprocess_kwargs)
#
# plat.subplots_evol_at_sol(Ra_ls[-2:], eta_ls[2], regime_grid=regime_grid_td[2,-2:], save=True, t1_grid=t1_grid[2,-2:], load_grid=True,
#                         psuffixes=['_T'], fig_path=fig_path,
#                         fname='evol-eta1e7', fig_fmt=fig_fmt, end_grid=end_grid[2,-2:], normtime=True, labelsize=14, xlabel=r'Time',
#                         ylabels=None, keys=['y_L', 'T_l', 'dT_rh', 'delta_rh', 'h_components'], title='', xlabelpad=8,
#                         ylabelpad=8, markers=None, markersize=24, cmap=cmap, colour_by='Ra',
#                         data_path=data_path, regime_names=regime_names_td, postprocess_kwargs=postprocess_kwargs)
#
# plat.subplots_evol_at_sol(Ra_ls[-2:], eta_ls[3], regime_grid=regime_grid_td[3,-2:], save=True, t1_grid=t1_grid[3,-2:], load_grid=True,
#                         psuffixes=['_T'], fig_path=fig_path,
#                         fname='evol-eta1e8', fig_fmt=fig_fmt, end_grid=end_grid[3,-2:], normtime=True, labelsize=14, xlabel=r'Time',
#                         ylabels=None, keys=['y_L', 'T_l',  'dT_rh', 'delta_rh', 'h_components'], title='', xlabelpad=8,
#                         ylabelpad=8, markers=None, markersize=24, cmap=cmap, colour_by='Ra',
#                         data_path=data_path, regime_names=regime_names_td, postprocess_kwargs=postprocess_kwargs)
