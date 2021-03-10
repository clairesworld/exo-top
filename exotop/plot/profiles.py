from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, load_grid, data_path, fig_path, fig_fmt
from postaspect import plt_aspect as plat
from postaspect.aspect_post import T_parameters_at_sol, get_cases_list, pickleio, h_at_ts, Nu_at_ts  # noqa: E402
from postaspect import ani_aspect as anims
import os

# cases = ['Ra1e7-eta1e5-wide', 'Ra3e7-eta1e5-wide', 'Ra1e8-eta1e5-wide', 'Ra3e8-eta1e5-wide', 'Ra3e8-eta1e6-wide']
# for case in cases:
#     plat.plot_velocity_profile(case, n=None, data_path=data_path, labelsize=16, fname=case+'_u-mean',
#                              save=True, fig_path=fig_path+'case_profiles/', fig_fmt=fig_fmt)

i_plot = range(3,4) # list(range(len(eta_ls)))  #
for ii, eta in enumerate(eta_ls):  # eta_ls
    if ii in i_plot:
        for jj, Ra in enumerate(Ra_ls):
            ### plot just T(z) at final timestep (or mean?), and extract all temperature params
            case = 'Ra' + Ra + '-eta' + eta + end_grid[ii, jj]
            T_params = pickleio(case, suffix='_T', postprocess_functions=[T_parameters_at_sol],
                                     t1=t1_grid[ii, jj], load=False, data_path=data_path)

            fig, ax = plat.plot_T_profile(case, T_params=T_params, n=-1, data_path=data_path, t1=t1_grid[ii, jj],
                                          setylabel=True,
                                          setxlabel=True, save=True, fig_path=fig_path+'case_profiles/', fend='_T-z', fig_fmt=fig_fmt,
                                          legend=True,
                                          load=load_grid[ii, jj])

            ### look at distribution of T parameters over time

            # plat.plot_pdf(case, df=T_params, keys=['h_components'], fig_path=fig_path, save=True,
            #             legend=False, c_list=['xkcd:forest green'], fig_fmt=fig_fmt,
            #             xlabel=r'$\alpha\Delta T_{rh} \delta_{rh}$', fend='T_hist')


            for jj, etastr in enumerate(eta_ls):
                if jj <= 20:
                    cases, cases_var = get_cases_list(Ra_ls, etastr, end_grid[jj])
                    for ii, case in enumerate(cases):
                        if (os.path.exists(data_path + 'output-' + case)) and (ii >= 4):
                            fig, ax = anims.static_h(case, data_path=data_path, fig_path=fig_path, labelsize=30,
                                                     ticksize=16, c='k', save=True, i_ts=-1, fig=None, ax=None)
                            fig, ax = anims.static_T_prof(case, data_path=data_path, fig_path=fig_path, labelsize=30,
                                                     ticksize=16, c='k', save=True, avg=True, fig=None, ax=None)
                            fig, ax = anims.static_T_field(case, data_path=data_path, fig_path=fig_path, labelsize=30,
                                                           ticksize=16, cmap='gist_heat',
                                                          shading='nearest', save=True, i_n=0, avg=False, c='k', dark=True,
                                                          fig=None, ax=None)

                            # anims.animate_T_field(case, data_path=data_path_bullard, fig_path=fig_path_bullard+'animations/', labelsize=30, ticksize=16,
                            #                 shading='nearest',#'gouraud',
                            #                 cmap='gist_heat')
                            # anims.animate_T_prof(case, data_path=data_path_bullard, fig_path=fig_path_bullard + 'animations/', labelsize=30, ticksize=16)
                            # anims.animate_h(case, data_path=data_path_bullard, fig_path=fig_path_bullard + 'animations/', labelsize=30, ticksize=16)
                            print('finished case')