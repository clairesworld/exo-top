import sys
sys.path.insert(0, '/home/cmg76/Works/exo-top/')
# import numpy as np
from exotop import aspect_scalings as sc
from setup_postprocessing import Ra_ls, eta_ls, t1, end, data_path, fig_path, c_rms, c_peak, fig_fmt, regime_grid

# c_regimes = [(0.2208, 0.5455, 0.2338), (0.984, 0.925, 0.365), (0.6882, 0.1059, 0.2059)] # ['#112a12',  '#fbec5d', '#EA2446']
c_regimes = ['xkcd:sage green', 'xkcd:blood red', 'xkcd:dark violet']

print('\nAssuming you have properly saved h values for the right time period\n')

### plot h scalings with Ra

sc.subplots_h_vs(Ra_ls, eta_ls, regime_grid, t1=t1, c_regimes=c_regimes, load=True, fit=True,
                 labelsize=14, xlim=(1e6, 3e8), ylim=(1, 9), xlabelpad=8, ylabelpad=10, logx=True, logy=True,
                 fname='h_Ra_all', xlabel='Ra', hscale=2e-5 * 2700 * 2890, ylabel='dynamic topography (km)',
                 data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt, save=True)

### plot h scalings - with dT_m*delta*alpha

sc.subplots_h_vs(Ra_ls, eta_ls, regime_grid, t1=t1, c_regimes=c_regimes, load=True, fit=True,
                 labelsize=14, xlim=(3e-8, 7e-7), ylim=(8e-3, 6e-2), xlabelpad=8, ylabelpad=12, logx=True, logy=True,
                 T_components=True,
                 fname='h_T_all', xlabel=r'$\alpha \delta_{rh} \Delta T_{rh}$', ylabel='dynamic topography',
                 showallscatter=True, data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt, save=True)

### plot scalings of other output parameters with Ra

# sc.subplots_vs_Ra(Ra_ls, eta_ls, t1=t1, keys=['Nu', 'delta_0', 'T_i'], load=True, fit=True,
#                   fname='delta-Nu-Ti', compare_pub=None, #sc.moresi95
#                   data_path=data_path, fig_path=fig_path, fig_fmt=fig_fmt, save=True)

## Show convective regimes in parameter space
## example of transitional is Ra3e7 eta1e6 - still has regular cells

# sc.plot_convection_regimes(Ra_ls, eta_ls, regime_grid, fname='regimes', data_path=data_path, fig_path=fig_path,
#                            load='auto', save=True, fig_fmt=fig_fmt,
#                            overploth=True, nlevels=14, clist=c_regimes, cmap_contours='autumn')
