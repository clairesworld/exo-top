import sys
# sys.path.insert(0, '/home/cmg76/Works/exo-top/')
# sys.path.insert(0, '//')
from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path_bullard, data_path_home, fig_path_bullard, fig_path_home, c_rms, c_peak, \
    c_regimes_td, regime_grid_td, regime_names_td, load_grid, alpha_m, p_Earth  # noqa: E402
from model_1D.parameters import M_E  # noqa: E402
from postaspect.aspect_post import get_cases_list  # noqa: E402
from postaspect import plt_1Dscaling as plt1  # noqa: E402
from postaspect import ani_aspect as anims  # noqa: E402

figsize = (12, 4)
labelsize = 30
ticksize = 18
xlabel = r'Planet mass ($M_E$)'
fname = ['ani_h_time', 'ani_drh_time', 'ani_Trh_time']
ylabel = [r'$\Delta h_{rms}^{\prime}$', r'$\delta_{rh}^{\prime}$', r'$\Delta T_{rh}^{\prime}$']
y_param = ['dyn_top_rms', 'delta_rh', 'dT_rh']
textb = [True, False, False]
yticks = [[0.001, 0.01],[0.002, 0.01, 0.04],[0.1, 0.2]]

for name, y, yt, label, tb in zip(fname, y_param, yticks, ylabel, textb):
    plt1.animate_Ra_time(default='Earthbaseline', fig_path=fig_path_home, figsize=figsize, labelsize=labelsize,
                         ylabelpad=20,
                         xlabelpad=20, ticksize=ticksize, fname=name, ms_scat=100, c_scat='xkcd:coral',
                         secy=False, yticks=yt, i_static=-1,
                         data_path=data_path_home,
                         x_min=0.1 * M_E,

                         # x_min=0.03 * M_E,
                         x_max=6 * M_E, y_param=y, x_param='M_p', x2_param='Ra_i_eff', xscale=M_E ** -1, x2scale=1, yscale=1, fps=15,
                         c='xkcd:pale green', anim_tpref='', tf=9.9, ti=2, ylabel=label,
                         x_min_planets=0.1 * M_E, ini_dict=None, logx=True, dark=True, lw=3, nplanets=20,
                         alpha_lines=1, text=tb, x_test2=0.107 * M_E,
                         xlabel=xlabel, x2label=r'Ra$_{i, eff}$', marker_scat='s')


# dimensionalised secondary y axis
# plt1.animate_Ra_time(default='Earthbaseline', fig_path=fig_path_home, figsize=figsize, labelsize=labelsize, ylabelpad=20,
#                       xlabelpad=20, ticksize=ticksize, fname='ani_h_time_dim', ms_scat=100, c_scat='xkcd:coral', secy=True,
#                       data_path=data_path_home, x_min=0.02745715027114846 * M_E, x_max=6 * M_E, y_param='dyn_top_rms',
#                       x_param='M_p', x2_param='Ra_i_eff', xscale=M_E ** -1, x2scale=1, yscale=1, fps=15,
#                       c='xkcd:pale green', anim_tpref='', tf=9.9, ti=2, ylabel=r'$\Delta h_{rms}^{\prime}$',
#                       x_min_planets=0.1 * M_E, ini_dict=None, logx=True, dark=True, lw=3, nplanets=20, alpha_lines=1,
#                       xlabel=r'$M_p$ ($M_E$)', x2label=r'Ra$_{i, eff}$', marker_scat='s')

#
# for jj, etastr in enumerate(eta_ls):
#     if jj <= 20:
#         cases, cases_var = get_cases_list(Ra_ls, etastr, end_grid[jj])
#         for ii, case in enumerate(cases):
#             if (os.path.exists(data_path_bullard + 'output-' + case)) and (ii >= 4):
#                 anims.animate_T_field(case, data_path=data_path_bullard, fig_path=fig_path_bullard+'animations/', labelsize=30, ticksize=16,
#                                 shading='nearest',#'gouraud',
#                                 cmap='gist_heat')
#                 # anims.animate_T_prof(case, data_path=data_path_bullard, fig_path=fig_path_bullard + 'animations/', labelsize=30, ticksize=16)
#                 # anims.animate_h(case, data_path=data_path_bullard, fig_path=fig_path_bullard + 'animations/', labelsize=30, ticksize=16)
#                 print('finished case')
