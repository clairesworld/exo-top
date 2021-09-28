from model_1D import the_results as plt1
import numpy as np
from model_1D.parameters import M_E, sec2Gyr

fig_path = ''
labelsize = 16
ticksize = 12
legsize = 16

fig, ax = plt1.plot_h_v_obvs(age=4.5, labelsize=labelsize, legsize=legsize, ticksize=ticksize, xlabelpad=20,
                             fig_path=fig_path, leg=True,
                  ylabel='$\Delta h_{rms}$ (m)', log=False, nplanets=20, save=False, x_vars=['t', 'M_p', 'H_0', 'CMF'],
                                      units=['Gyr', '$M_E$', 'pW kg$^{-1}$', 'CMF'],
                  x_range=[(1.5, 4.5), (0.1 * M_E, 6 * M_E), (10e-12, 40e-12),
                           (0.01, 0.4)], xticks=None,
                  xscales=[sec2Gyr, M_E ** -1, 1e12, 1], ylim=(500, 1000), yticks=[500, 600, 700, 800, 900, 1000],
                  xlabels=['Age\n(Gyr)', 'Planet mass\n($M_E$)',
                           'Rad. heating $t_0$\n(pW kg$^{-1}$)',
                           'Core mass fraction', 'Activation energy\n(kJ mol$^{-1}$)'],)