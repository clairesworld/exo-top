import numpy as np

Ra_ls = ['1e6', '3e6', '1e7', '3e7', '1e8', '3e8']  # 3e5 no convection
eta_ls = ['1e5', '1e6', '1e7', '1e8']
t1 = np.array([[0.5, 0.3, 0.25, 0.7, 0.2, 0.05],  # eta 1e5
               [0.9, 0.3, 0.3, 0.4, 0.3, 0.05],  # eta 1e6
               [0, 0.7, 0.35, 0.4, 0.3, 1],  # eta 1e7
               [0, 0, 0.6, 0.55, 1, 1]])  # eta 1e8
end = np.empty_like(t1, dtype=object)
end[:] = "-wide"
end[2:3, 5] = "-wide-ascii"  # which runs use ascii initialisation
print('end', end)
