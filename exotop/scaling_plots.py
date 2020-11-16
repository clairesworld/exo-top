import sys
sys.path.insert(0,'/home/cmg76/Works/exo-top/')
# import numpy as np
from exotop import aspect_scalings as sc
from setup_postprocessing import Ra_ls, eta_ls, t1, end, data_path, fig_path, c_rms, c_peak, fig_fmt, regime_grid

#c_regimes = [(0.2208, 0.5455, 0.2338), (0.984, 0.925, 0.365), (0.6882, 0.1059, 0.2059)] # ['#112a12',  '#fbec5d', '#EA2446']
c_regimes = ['xkcd:sage green', 'xkcd:blood red' , 'xkcd:dark violet']

### plot h scalings with Ra

print('\nAssuming you have properly saved h values for the right time period\n')
fig, axes = sc.subplots_h_vs(Ra_ls, eta_ls, regime_grid, c_regimes, fit=True, t1=t1,
                             load=True, xlim=(8e5,1.3e8), ylim=None, fig_fmt=fig_fmt,
                             save=True, fname='h_Ra_all', xlabel='Ra', hscale=2e-5*2700*2890,
                             ylabel='dynamic topography (km)', data_path=data_path, fig_path=fig_path)


### plot h scalings - with dT_m*delta*d_m

fig, ax = sc.subplots_h_vs(Ra_ls, eta_ls, regime_grid=regime_grid, c_regimes=c_regimes, save=True, t1=t1,
                           fit=True, load='auto', T_components=True, data_path=data_path,
                           fig_path=fig_path, fname='h_T_all', fig_fmt=fig_fmt,
                           ylim=(3e-3, 7e-2), labelsize=14, xlim=(3e-8, 7e-7),
                           xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography', logx=True, logy=True,
                           showallscatter=True, xlabelpad=8, ylabelpad=0)


# say what convective regimes are after looking at temperature fields
# example of transitional is Ra3e7 eta1e6 - still even cells
# sc.plot_convection_regimes(Ra_ls, eta_ls, regime_grid, fname='regimes', data_path=data_path, fig_path=fig_path,
#                            loadpickle=True, dumppickle=False, save=True,
#                            overploth=True, nlevels=14, clist=c_regimes, cmap_contours='autumn')

### make bl scaling plot

sc.plot_multi_Ra_scaling(Ra_ls, eta_ls, keys=['Nu', 'delta_0', 'T_i'], save=True, compare_pub=sc.solomatov95,
                         fname='bl-Nu', t1=t1, fit=True, fig_fmt=fig_fmt, data_path=data_path, fig_path=fig_path,)

