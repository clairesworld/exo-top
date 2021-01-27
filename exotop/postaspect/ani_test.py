import sys
sys.path.insert(0, '/home/cmg76/Works/exo-top/')
sys.path.insert(0, '/home/claire/Works/exo-top/')
from exotop.postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, c_rms, c_peak, \
    c_regimes_td, fig_fmt, \
    regime_grid_td, regime_names_td, load_grid, alpha_m, p_Earth  # noqa: E402
from exotop.model_1D.parameters import M_E  # noqa: E402
from exotop.postaspect import aspect_combined as combo  # noqa: E402

fig_path_home = '/home/claire/Works/exo-top/exotop/figs_scratch/'
data_path_home = '/home/claire/Works/aspect/runs/model-output/'

figsize = (12, 4)
labelsize = 30
ticksize = 18

combo.animate_Ra_time(default='Earthbaseline', fig_path=fig_path_home, figsize=figsize, labelsize=labelsize, ylabelpad=20,
                      xlabelpad=20, ticksize=ticksize, fname='ani_h_time_dim', ms_scat=100, c_scat='xkcd:coral', secy=True,
                      data_path=data_path_home, x_min=0.02745715027114846 * M_E, x_max=6 * M_E, y_param='dyn_top_rms',
                      x_param='M_p', x2_param='Ra_i_eff', xscale=M_E ** -1, x2scale=1, yscale=1, fps=15,
                      c='xkcd:pale green', anim_tpref='', tf=9.9, ti=2, ylabel=r'$\Delta h_{rms}^{\prime}$',
                      x_min_planets=0.1 * M_E, ini_dict=None, logx=True, dark=True, lw=3, nplanets=20, alpha_lines=1,
                      xlabel=r'$M_p$ ($M_E$)', x2label=r'Ra$_{i, eff}$', marker_scat='s')

combo.animate_Ra_time(default='Earthbaseline', fig_path=fig_path_home, figsize=figsize, labelsize=labelsize, ylabelpad=20,
                      xlabelpad=20, ticksize=ticksize, fname='ani_h_time', ms_scat=100, c_scat='xkcd:coral', secy=False,
                      data_path=data_path_home, x_min=0.02745715027114846 * M_E, x_max=6 * M_E, y_param='dyn_top_rms',
                      x_param='M_p', x2_param='Ra_i_eff', xscale=M_E ** -1, x2scale=1, yscale=1, fps=15,
                      c='xkcd:pale green', anim_tpref='', tf=9.9, ti=2, ylabel=r'$\Delta h_{rms}^{\prime}$',
                      x_min_planets=0.1 * M_E, ini_dict=None, logx=True, dark=True, lw=3, nplanets=20, alpha_lines=1,
                      xlabel=r'$M_p$ ($M_E$)', x2label=r'Ra$_{i, eff}$', marker_scat='s')

combo.animate_Ra_time(default='Earthbaseline', fig_path=fig_path_home, figsize=figsize, labelsize=labelsize, ylabelpad=20,
                      xlabelpad=20, ticksize=ticksize, fname='ani_drh_time', ms_scat=100, c_scat='xkcd:coral', secy=False,
                      data_path=data_path_home, x_min=0.02745715027114846 * M_E, x_max=6 * M_E, y_param='delta_rh',
                      x_param='M_p', x2_param='Ra_i_eff', xscale=M_E ** -1, x2scale=1, yscale=1, fps=15,
                      c='xkcd:pale green', anim_tpref='', tf=9.9, ti=2, ylabel=r'$\delta_{rh}^{\prime}$',
                      x_min_planets=0.1 * M_E, ini_dict=None, logx=True, dark=True, lw=3, nplanets=20, alpha_lines=1,
                      xlabel=r'$M_p$ ($M_E$)', x2label=r'Ra$_{i, eff}$', marker_scat='s', text=False)

combo.animate_Ra_time(default='Earthbaseline', fig_path=fig_path_home, figsize=figsize, labelsize=labelsize, ylabelpad=20,
                      xlabelpad=20, ticksize=ticksize, fname='ani_Trh_time', ms_scat=100, c_scat='xkcd:coral', secy=False,
                      data_path=data_path_home, x_min=0.02745715027114846 * M_E, x_max=6 * M_E, y_param='dT_rh',
                      x_param='M_p', x2_param='Ra_i_eff', xscale=M_E ** -1, x2scale=1, yscale=1, fps=15,
                      c='xkcd:pale green', anim_tpref='', tf=9.9, ti=2, ylabel=r'$\Delta T_{rh}^{\prime}$',
                      x_min_planets=0.1 * M_E, ini_dict=None, logx=True, dark=True, lw=3, nplanets=20, alpha_lines=1,
                      xlabel=r'$M_p$ ($M_E$)', x2label=r'Ra$_{i, eff}$', marker_scat='s', text=False)