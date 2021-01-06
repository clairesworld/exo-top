""" ASPECT runs: setup for postprocessing, defining variables  """

import numpy as np
from matplotlib import use as matplotlibuse
matplotlibuse('Agg')

default_load_value = True
check_new = True  # always check for new timesteps with ongoing runs

alpha_m = 1
Ra_ls = np.array(['1e6', '3e6', '1e7', '3e7', '1e8', '3e8'])  # 3e5 no convection
eta_ls = np.array(['1e5', '1e6', '1e7', '1e8'])

# Ra from 1e6 to 3e8
t1_grid = np.array([[0.5, 0.3, 0.25, 0.7, 0.2, 0.0590015],  # eta 1e5
                    [0.9, 0.3, 0.3, 0.4, 0.3, 0.055083],  # eta 1e6
                    [0, 0.7, 0.35, 0.4, 0.3, 0.06],  # eta 1e7
                    [0, 0, 0.6, 0.55, 0.065, 0.06636]])  # eta 1e8

# Ra from 1e6 to 3e8
regime_grid_td = np.array([['steady', 'steady', 'steady', 'trans.', 'sluggish', 'sluggish'],  # eta 1e5
                           ['steady', 'steady', 'steady', 'trans.', 'chaotic', 'chaotic'],  # eta 1e6
                           ['no convection', 'steady', 'steady', 'trans.', 'chaotic', 'chaotic'],  # eta 1e7
                           ['no convection', 'no convection', 'steady', 'trans.', 'chaotic', 'chaotic']])  # eta 1e8
regime_names_td = ['steady', 'trans.', 'chaotic']
c_regimes_td = ['xkcd:sage green', 'xkcd:blood red', 'xkcd:azure']

load_grid = np.empty_like(t1_grid, dtype=object)
load_grid[:] = default_load_value
if check_new:
    load_grid[:, 5] = 'auto'   # ongoing runs: Ra3e8
    load_grid[3, 4] = 'auto'  # ongoing runs: Ra1e8 eta1e8

# string endings for case names
end_grid = np.empty_like(t1_grid, dtype=object)
end_grid[:] = "-wide"
end_grid[2:, 5] = "-wide-ascii"  # which runs use ascii initialisation
end_grid[3, 4] = "-wide-ascii"  # which runs use ascii initialisation

data_path_bullard = '/raid1/cmg76/aspect/model-output/'
fig_path_bullard = '/raid1/cmg76/aspect/figs/'
data_path = data_path_bullard
fig_path = fig_path_bullard
fig_fmt = '.png'

c_rms = 'xkcd:forest green'
c_peak = 'xkcd:periwinkle'

p_Earth = {'alpha_m': 1.35e-5, 'dT_m': 2700, 'd_m': 2800}   # dimensionalisation factors