import sys
sys.path.insert(0,'/home/cmg76/Works/exo-top/')
import numpy as np
from exotop import aspect_scalings as sc
from setup_postprocessing import Ra_ls, eta_ls, t1, end

# Ra from 1e6 to 3e8
regime_grid = np.array([['steady', 'steady' , 'steady', 'chaotic', 'sluggish', 'not converged'], # eta 1e5
               ['steady', 'steady', 'steady', 'transitional', 'chaotic', 'not converged'],  # eta 1e6 - 
               ['no convection', 'steady', 'steady', 'transitional', 'chaotic', 'not converged'], # eta 1e7
               ['no convection', 'no convection', 'steady', 'transitional', 'not converged', 'not converged']])  # eta 1e8
#c_regimes = [(0.2208, 0.5455, 0.2338), (0.984, 0.925, 0.365), (0.6882, 0.1059, 0.2059)] # ['#112a12',  '#fbec5d', '#EA2446']
c_regimes = ['xkcd:sage green', 'xkcd:blood red' , 'xkcd:dark violet']

# plot h scalings - with Ra

# print('assuming you have properly saved h values for the right time period')

# fig, axes = sc.subplots_h_vs(Ra_ls, eta_ls, regime_grid, c_regimes, fit=True, t1=t1,
#                              loadpickle=True, dumppickle=False, xlim=(8e5,1.3e8), ylim=None,
#                              save=True, fname='h_Ra_all.png', xlabel='Ra', hscale=2e-5*2700*2890,
#                             ylabel='dynamic topography (km)')

# print('Ra scaling complete')


# plot h scalings - with dT_m*delta*d_m   
# case = 'Ra3e7-eta1e5-wide'
# _ = sc.pd_h_components(case, t1=0, data_path='model-output/', fig_path='/raid1/cmg76/aspect/figs/', sigma=2,
#                     picklefrom=case+'_pdx.pkl', plotTz=False, savefig=True, 
#                     settitle=True, setxlabel=True, c='xkcd:pale purple', params_list=None,
#                     legend=True, plotpd=True, labelsize=16)

# for ii, eta in enumerate(eta_ls):
#     sc.plot_h_vs_Td(Ra=['1e6', '3e6', '1e7', '3e7', '1e8'], eta=eta, t1=t1[ii], 
#                   loadpickle=True, loadpicklex=True,
#                   save=True, fname='h_T_'+eta+'.png', sigma=2, 
#                   labelsize=16, xlabel=r'$\delta_u \Delta T_{rh}$', title='',
#                   c_peak='xkcd:forest green', c_rms='xkcd:periwinkle',legend=True,
#                   fit=False, logx=True, logy=True,
#                   xlim=None,)

# fig, ax = sc.subplots_h_vs(Ra_ls, eta_ls, regime_grid=regime_grid, c_regimes=c_regimes, x_components=True,
#                            t1=t1, path='model-output/', fig_path='/raid1/cmg76/aspect/figs/', 
#                            loadpickle=True, dumppickle=False, loadpicklex=True,
#                            save=True, fname='h_T_all.png', showallscatter=True,
#                            labelsize=14, xlabel=r'$\alpha \delta_u \Delta T_{rh}$', ylabel='dynamic topography', 
#                            fit=True, logx=True, logy=True, ylim=(3e-3,7e-2), xlim=(3e-8, 7e-7),
#                           xlabelpad=8, ylabelpad=0)
    
        
# for ii, eta in enumerate(eta_ls):
#     try:
#         fig, ax = sc.plot_h_vs(Ra=includeRa[ii], eta=eta, t1=t1[ii], sigma=2, title='$\Delta \eta$='+eta, 
#                                xlabel='$\Delta T_{rh} \delta_u d_m$',
#                                loadpickle=True, dumppickle=dump,  plotpd=False, fname='h-x_eta'+eta+'.png', x_components=True,
#                                fit=False, loadpicklex=False, plotTz=True)
#         print(eta, 'dT_m*delta*d_m scaling complete')
#     except Exception as e:
#         exc_type, exc_obj, exc_tb = sys.exc_info()
#         fename = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
#         print(exc_type, fename, exc_tb.tb_lineno, e)


# say what convective regimes are after looking at temperature fields
# example of transitional is Ra3e7 eta1e6 - still even cells
# sc.plot_convection_regimes(Ra_ls, eta_ls, regime_grid, fname='regimes.png', 
#                            loadpickle=True, dumppickle=False, save=True,
#                            overploth=True, nlevels=14, clist=c_regimes, cmap_contours='autumn')

# make bl scaling plot

sc.plot_bl_Nu_scaling(Ra_ls, eta_ls, loadpicklex=True, save=True, compare_pub='S95',
                      fname='bl-Nu.png', t1=t1, fitdelta=False, fitNu=True)
