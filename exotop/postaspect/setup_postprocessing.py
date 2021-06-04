""" ASPECT runs: setup for postprocessing, defining variables  """

import numpy as np
from matplotlib import use as matplotlibuse
from matplotlib import rc
from matplotlib.pyplot import rcParams

matplotlibuse('Agg')  # turn on for running over ssh
# rc('text', usetex=True)  # turn off for running over ssh
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'CMU Serif'

default_load_value = True
check_new = True  # always check for new timesteps with ongoing runs

# model input parameters
alpha_m = 1
X_extent = 8
Y_extent = 1
postprocess_kwargs = {'X_extent': X_extent, 'Y_extent': Y_extent, 'alpha_m': alpha_m}

Ra_ls = np.array(['1e6', '3e6', '1e7', '3e7', '1e8', '2e8', '3e8'])  # 3e5 no convection
eta_ls = np.array(['1e5', '1e6', '1e7', '1e8', '1e9'])

# Ra from 1e6 to 3e8
t1_grid = np.array([[0.5, 0.3, 0.25, 0.7, 0.2, 0, 0.0590015],  # eta 1e5
                    [0.9, 0.3, 0.3, 0.4, 0.3, 0, 0.07],  # eta 1e6
                    [0, 0.7, 0.35, 0.4, 0.3, 0, 0.075],  # eta 1e7
                    [0, 0, 0.6, 0.55, 0.09, 0, 0.085],  # eta 1e8
                    [0, 0, 0, 0, 0.08, 0, 0.07]])  # eta 1e9 ???

# Ra from 1e6 to 3e8
regime_grid_td = np.array([['steady', 'steady', 'steady', 'trans.', 'sluggish', 'not ready', 'sluggish'],  # eta 1e5
                           ['steady', 'steady', 'steady', 'trans.', 'chaotic', 'not ready', 'chaotic'],  # eta 1e6
                           ['no convection', 'steady', 'steady', 'trans.', 'chaotic', 'not ready', 'chaotic'],  # eta 1e7
                           ['no convection', 'no convection', 'steady', 'trans.', 'chaotic', 'not ready', 'chaotic'],  # eta 1e8
                           ['not ran', 'not ran', 'not ran', 'not ran', 'not ready', 'not ran', 'not ready']  # eta 1e9???
                           ])
regime_names_td = ['steady', 'trans.', 'chaotic', 'not ready']
c_regimes_td = ['xkcd:sage green', 'xkcd:blood red', 'xkcd:azure']

load_grid = np.empty_like(t1_grid, dtype=object)
load_grid[:] = default_load_value
if check_new:
    load_grid[1:, 5:] = 'auto'   # ongoing runs: Ra3e8
    load_grid[3, 4] = 'auto'  # ongoing runs: Ra1e8 eta1e8
    load_grid[4, 4:] = 'auto'  # ongoing runs: new eta1e9 runs
    load_grid[:, -2] = False  # 2e8 initialisation

# string endings for case names
end_grid = np.empty_like(t1_grid, dtype=object)
end_grid[:] = "-wide"
end_grid[2:, 5:] = "-wide-ascii"  # which runs use ascii initialisation
end_grid[3, 4] = "-wide-ascii"  # which runs use ascii initialisation
end_grid[4, 4:] = "-wide-ascii"  # ongoing runs: new eta1e9 runs

data_path_bullard = '/raid1/cmg76/aspect/model-output/'
fig_path_bullard = '/raid1/cmg76/aspect/figs/'
fig_path_home = '/home/claire/Works/exo-top/exotop/figs_scratch/'
data_path_home = '/home/claire/Works/aspect/runs/model-output/'
data_path = data_path_bullard
fig_path = fig_path_bullard
fig_fmt = '.png'
cmap_path = 'exotop/plot/cmaps/'
benchmark_path = '/home/cmg76/Works/exo-top/benchmarks/'

c_rms = 'xkcd:forest green'
c_peak = 'xkcd:periwinkle'
highlight_colour = 'xkcd:coral'

p_Earth = {'alpha_m': 1.35e-5, 'dT_m': 2700, 'd_m': 2800e3}   # dimensionalisation factors
p_Venus = {'alpha_m': 2e-5, 'dT_m': 3000, 'd_m': 2723e3}   # dimensionalisation factors Weller & Kiefer 2020
