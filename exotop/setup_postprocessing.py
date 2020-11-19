import numpy as np

default_load_value = True
check_new = True  # always check for new timesteps with ongoing runs

alpha_m = 1.35e-5
Ra_ls = np.array(['1e6', '3e6', '1e7', '3e7', '1e8', '3e8'])  # 3e5 no convection
eta_ls = np.array(['1e5', '1e6', '1e7', '1e8'])

# Ra from 1e6 to 3e8
t1_grid = np.array([[0.5, 0.3, 0.25, 0.7, 0.2, 0.055],  # eta 1e5
                    [0.9, 0.3, 0.3, 0.4, 0.3, 0.0559],  # eta 1e6
                    [0, 0.7, 0.35, 0.4, 0.3, 1],  # eta 1e7
                    [0, 0, 0.6, 0.55, 0.065, 1]])  # eta 1e8

# Ra from 1e6 to 3e8
regime_grid = np.array([['steady', 'steady', 'steady', 'chaotic', 'chaotic', 'chaotic'],  # eta 1e5
               ['steady', 'steady', 'steady', 'trans.', 'chaotic', 'chaotic'],  # eta 1e6
               ['no convection', 'steady', 'steady', 'trans.', 'chaotic', 'not converged'],  # eta 1e7
               ['no convection', 'no convection', 'steady', 'trans.', 'not converged', 'not converged']])  # eta 1e8
regime_names = ['steady', 'trans.', 'chaotic']

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
# c_regimes = [(0.2208, 0.5455, 0.2338), (0.984, 0.925, 0.365), (0.6882, 0.1059, 0.2059)] # ['#112a12',  '#fbec5d', '#EA2446']
c_regimes = ['xkcd:sage green', 'xkcd:blood red', 'xkcd:dark violet']