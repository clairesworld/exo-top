""" Plot new 1D thermal model topography after processing of ASPECT runs """

from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path_home, fig_path_home, c_rms, c_peak, \
    c_regimes_td, fig_fmt, \
    regime_grid_td, regime_names_td, load_grid, alpha_m, p_Earth  # noqa: E402
from model_1D.parameters import M_E  # noqa: E402
from postaspect import plt_1Dscaling as plt1  # noqa: E402

M_p_min = 0.03 * M_E
M_p_max = 6 * M_E
M_p_min_display = 0.1 * M_E
tf = 10  # time in Gyr

#
# plt1.plot_1D_evolutions(default='Earthbaseline', nplanets=45, labelsize=23, ticksize=16, clabelpad=4, save=True,
#                          fname='evol_summary_mass', zmin=M_p_min_display, zmax=M_p_max, zname='M_p', zlabel=r'$M_p$ ($M_E$)',
#                          zscale=M_E**-1, backwards_cooling=False,
#                          ylabels_dict=None, fig_path=fig_path_home, fig_fmt=fig_fmt,
#                          age=tf, ini_dict=None, cmap='rainbow')
#
#
# plt1.plot_1D_evolutions(default='Earthbaseline', nplanets=45, labelsize=23, ticksize=16, clabelpad=4, save=True,
#                          fname='evol_summary_Ea_hi', zmin=300e3, zmax=540e3, zname='Ea', zlabel=r'$E_a$ (kJ/mol)',
#                          zscale=1e-3, backwards_cooling=False,
#                          ylabels_dict=None, fig_path=fig_path_home, fig_fmt=fig_fmt,
#                          age=tf, ini_dict=None, cmap='rainbow')
#
#
# plt1.plot_1D_evolutions(default='Earthbaseline', nplanets=45, labelsize=23, ticksize=16, clabelpad=4, save=True,
#                          fname='evol_summary_Ea_lo', zmin=200e3, zmax=350e3, zname='Ea', zlabel=r'$E_a$ (kJ/mol)',
#                          zscale=1e-3, backwards_cooling=False,
#                          ylabels_dict=None, fig_path=fig_path_home, fig_fmt=fig_fmt,
#                          age=tf, ini_dict=None, cmap='rainbow')

plt1.plot_h_heuristic_variations(default='Earthbaseline', x_param='M_p', x_min=M_p_min_display, x_max=M_p_max,
                                 x_min_planets=None, verbose=False,
                                 xlabel=r'$M_p$ ($M_E$)', x2label=r'Ra$_{i, eff}$', save=True,
                                 xleglabel=r'$M_p$, 1D model', x2leglabel=r'Ra$_{i, eff}$, 1D model',
                                 fig_path=fig_path_home, fig_fmt=fig_fmt,
                                 legend=True, nplanets=20, age=tf,  #ini_dict={'D_l0':10e3},
                                 fname='heuristic_mass_dependence', nondimensional=False,
                                 c='xkcd:drab green', lw=3, labelsize=24, ticksize=12, legsize=12)

plt1.plot_h_heuristic_variations(default='Earthbaseline', x_param='M_p', x_min=M_p_min, x_max=M_p_max,
                                 x_min_planets=M_p_min_display,
                                 xlabel=r'$M_p$ ($M_E$)', x2label=r'Ra$_{i, eff}$', save=True,
                                 xleglabel=r'$M_p$, 1D model', x2leglabel=r'Ra$_{i, eff}$, 1D model',
                                 aspectleglabel=r'Ra$_{i, eff}$, 2D model',
                                 fig_path=fig_path_home, fig_fmt=fig_fmt,
                                 legend=True, nplanets=20, age=tf,
                                 fname='heuristic_mass_dependence_nondim', nondimensional=True,
                                 c='xkcd:drab green', lw=3, labelsize=24, ticksize=12, legsize=8,
                                 data_path=data_path_home,
                                 overplot_aspect_cases_x_param='Ra_i_eff',
                                 overplot_aspect_cases=['Ra3e8-eta1e8-wide-ascii', 'Ra3e8-eta1e7-wide-ascii',
                                                         'Ra3e8-eta1e6-wide', 'Ra1e8-eta1e8-wide-ascii',
                                                         'Ra1e8-eta1e7-wide', 'Ra1e8-eta1e6-wide']
                                 )



