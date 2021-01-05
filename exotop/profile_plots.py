import sys

sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, load_grid, data_path, fig_path, fig_fmt
from exotop.postaspect import aspect_scalings as sc

# cases = ['Ra1e7-eta1e5-wide', 'Ra3e7-eta1e5-wide', 'Ra1e8-eta1e5-wide', 'Ra3e8-eta1e5-wide', 'Ra3e8-eta1e6-wide']
# for case in cases:
#     sc.plot_velocity_profile(case, n=None, data_path=data_path, labelsize=16, fname=case+'_u-mean',
#                              save=True, fig_path=fig_path+'case_profiles/', fig_fmt=fig_fmt)

i_plot = range(3,4) # list(range(len(eta_ls)))  #
for ii, eta in enumerate(eta_ls):  # eta_ls
    if ii in i_plot:
        for jj, Ra in enumerate(Ra_ls):
            ### plot just T(z) at final timestep (or mean?), and extract all temperature params
            case = 'Ra' + Ra + '-eta' + eta + end_grid[ii, jj]
            T_params = sc.pickleio(case, suffix='_T', postprocess_functions=[sc.T_parameters_at_sol],
                                   t1=t1_grid[ii, jj], load=False, data_path=data_path)

            fig, ax = sc.plot_T_profile(case, T_params=T_params, n=-1, data_path=data_path, t1=t1_grid[ii, jj],
                                        setylabel=True,
                                        setxlabel=True, save=True, fig_path=fig_path+'case_profiles/', fend='_T-z', fig_fmt=fig_fmt,
                                        legend=True,
                                        load=load_grid[ii, jj])

            ### look at distribution of T parameters over time

            # sc.plot_pdf(case, df=T_params, keys=['h_components'], fig_path=fig_path, save=True,
            #             legend=False, c_list=['xkcd:forest green'], fig_fmt=fig_fmt,
            #             xlabel=r'$\alpha\Delta T_{rh} \delta_{rh}$', fend='T_hist')
