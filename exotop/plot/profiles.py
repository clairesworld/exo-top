from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, load_grid, data_path, fig_path, fig_fmt
from postaspect import plt_aspect as plat
from postaspect.aspect_post import T_parameters_at_sol, get_cases_list, pickleio, h_at_ts, Nu_at_ts  # noqa: E402
from postaspect import ani_aspect as anims
import os

# fig, ax = anims.static_uv_prof('Ra3e8-eta1e8-wide-ascii', data_path=data_path, fig_path=fig_path, labelsize=30,
#                          ticksize=16, c='k', save=True, avg=True, fig=None, ax=None,
#                                xlim=(0, 3000), xlabel='RMS Velocity')
#
# fig, ax = anims.animate_uv_prof('Ra3e8-eta1e8-wide-ascii', data_path=data_path, fig_path=fig_path, labelsize=30,
#                          ticksize=16, c='k', fig=None, ax=None,
#                                xlim=(0, 3000), xlabel='RMS Velocity')

# fig, axes = anims.static_h(case='Lees-Ra1e6-2D', data_path=data_path, fig_path=fig_path, labelsize=30, ticksize=16, dark=False,
#              xlabel='x', ylabel='h', return_artists=False, c='k', save=True, i_ts=43650, avg=False, #ylim=(-9e-2, 9e-2)
#  )

""" manuscript gridspec """
labelsize = 30
ticksize = 20
wspace = 0.2
hspace = 0.2

# fig, axes = anims.T_h_gridspec(case='Lees-Ra1e6-2D', data_path=data_path, fig_path=fig_path, labelsize=30, ticksize=20,
#                                cmap='gist_heat', save=True, c='k')
fig, axes = anims.T_h_gridspec(case='Ra1e8-eta1e8-wide-ascii', data_path=data_path, fig_path=fig_path,
                               labelsize=labelsize, ticksize=ticksize, wspace=wspace, hspace=hspace,
                               cmap='gist_heat', save=True, c='k')
fig, axes = anims.T_h_gridspec(case='Ra3e8-eta1e8-wide-ascii', data_path=data_path, fig_path=fig_path,
                               labelsize=labelsize, ticksize=ticksize, wspace=wspace, hspace=hspace,
                               cmap='gist_heat', save=True, c='k')

# cases = ['Ra1e7-eta1e5-wide', 'Ra3e7-eta1e5-wide', 'Ra1e8-eta1e5-wide', 'Ra3e8-eta1e5-wide', 'Ra3e8-eta1e6-wide']
# for case in cases:
#     plat.plot_velocity_profile(case, n=None, data_path=data_path, labelsize=16, fname=case+'_u-mean',
#                              save=True, fig_path=fig_path+'case_profiles/', fig_fmt=fig_fmt)
#
# i_plot = range(3,4) # list(range(len(eta_ls)))  #
# for ii, eta in enumerate(eta_ls):  # eta_ls
#     if ii in i_plot:
#         for jj, Ra in enumerate(Ra_ls):
#             case = 'Ra' + Ra + '-eta' + eta + end_grid[ii, jj]
#
#             ### plot just T(z) at final timestep (or mean?), and extract all temperature params
#             # T_params = pickleio(case, suffix='_T', postprocess_functions=[T_parameters_at_sol],
#             #                          t1=t1_grid[ii, jj], load=False, data_path=data_path)
#             #
#             # fig, ax = plat.plot_T_profile(case, T_params=T_params, n=-1, data_path=data_path, t1=t1_grid[ii, jj],
#             #                               setylabel=True,
#             #                               setxlabel=True, save=True, fig_path=fig_path+'case_profiles/', fend='_T-z', fig_fmt=fig_fmt,
#             #                               legend=True,
#             #                               load=load_grid[ii, jj])
#
#             ### look at distribution of T parameters over time
#
#             # plat.plot_pdf(case, df=T_params, keys=['h_components'], fig_path=fig_path, save=True,
#             #             legend=False, c_list=['xkcd:forest green'], fig_fmt=fig_fmt,
#             #             xlabel=r'$\alpha\Delta T_{rh} \delta_{rh}$', fend='T_hist')
#
# =
#             # fig, ax = anims.static_h(case, data_path=data_path, fig_path=fig_path, labelsize=30,
#             #                          ticksize=16, c='k', save=True, i_ts=-1, fig=None, ax=None)
#             # fig, ax = anims.static_T_prof(case, data_path=data_path, fig_path=fig_path, labelsize=30,
#             #                          ticksize=16, c='k', save=True, avg=True, fig=None, ax=None)
#             # fig, ax = anims.static_T_field(case, data_path=data_path, fig_path=fig_path, labelsize=30,
#             #                                ticksize=16, cmap='gist_heat',
#             #                               shading='nearest', save=True, i_n=0, avg=True, c='k', dark=True,
#             #                               fig=None, ax=None)
#
#             fig, axes = anims.T_h_gridspec(case, data_path=data_path, fig_path=fig_path, labelsize=30, ticksize=16,
#                                             cmap='gist_heat', save=True, c='k')
#
#             # anims.animate_T_field(case, data_path=data_path_bullard, fig_path=fig_path_bullard+'animations/', labelsize=30, ticksize=16,
#             #                 shading='nearest',#'gouraud',
#             #                 cmap='gist_heat')
#             # anims.animate_T_prof(case, data_path=data_path_bullard, fig_path=fig_path_bullard + 'animations/', labelsize=30, ticksize=16)
#             # anims.animate_h(case, data_path=data_path_bullard, fig_path=fig_path_bullard + 'animations/', labelsize=30, ticksize=16)
#
#             print('finished case')