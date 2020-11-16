import sys

sys.path.insert(0, '/home/cmg76/Works/exo-top/')
# import numpy as np
# from exotop import aspect_postprocessing2 as asp
from exotop import aspect_scalings as sc
from setup_postprocessing import Ra_ls, eta_ls, t1, end, data_path, fig_path, fig_fmt


for ii, eta in enumerate(eta_ls):  # eta_ls
    for jj, Ra in enumerate(Ra_ls):
        ### plot just T(z) at final timestep (or mean?), and extract all temperature params
        case = 'Ra' + Ra + '-eta' + eta + '-wide'
        T_params = sc.pickleio(case, suffix='_T', postprocess_functions=[sc.T_parameters_at_sol],t1=t1[ii, jj],
                               load='auto', data_path=data_path, fig_path=fig_path, **kwargs)

        sc.plot_T_params(case, T_params=T_params, n=-1,  # n='mean'
                         data_path=data_path, t1=t1[ii, jj],
                         setylabel=True, setxlabel=True, save=True,
                         fig_path=fig_path, fend='_T-z', fig_fmt=fig_fmt,
                         legend=True)

        ### look at distribution of T parameters over time

        sc.plot_pdf(case, df=T_params, keys=['h_components'], fig_path=fig_path, save=True,
                    legend=False, c_list=['xkcd:forest green'], fig_fmt=fig_fmt,
                    xlabel=r'$\alpha\Delta T_{rh} \delta_{rh}$', fend='T_hist')