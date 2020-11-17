import numpy as np

Ra_ls = np.array(['1e6', '3e6', '1e7', '3e7', '1e8', '3e8'])  # 3e5 no convection
eta_ls = np.array(['1e5', '1e6', '1e7', '1e8'])

# Ra from 1e6 to 3e8
t1 = np.array([[0.5, 0.3, 0.25, 0.7, 0.2, 0.055],  # eta 1e5
               [0.9, 0.3, 0.3, 0.4, 0.3, 0.055],  # eta 1e6
               [0, 0.7, 0.35, 0.4, 0.3, 1],  # eta 1e7
               [0, 0, 0.6, 0.55, 1, 1]])  # eta 1e8

# Ra from 1e6 to 3e8
regime_grid = np.array([['steady', 'steady', 'steady', 'chaotic', 'chaotic', 'chaotic'],  # eta 1e5
               ['steady', 'steady', 'steady', 'transitional', 'chaotic', 'chaotic'],  # eta 1e6 -
               ['no convection', 'steady', 'steady', 'transitional', 'chaotic', 'not converged'],  # eta 1e7
               ['no convection', 'no convection', 'steady', 'transitional', 'not converged', 'not converged']])  # eta 1e8

# case names
end = np.empty_like(t1, dtype=object)
end[:] = "-wide"
end[2:, 5] = "-wide-ascii"  # which runs use ascii initialisation
end[3, 4] = "-wide-ascii"  # which runs use ascii initialisation

data_path_bullard = '/raid1/cmg76/aspect/model-output/'
fig_path_bullard = '/raid1/cmg76/aspect/figs/'
data_path = data_path_bullard
fig_path = fig_path_bullard
fig_fmt = '.png'

c_rms = 'xkcd:forest green'
c_peak = 'xkcd:periwinkle'
