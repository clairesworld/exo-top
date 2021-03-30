import fractals as fract
import sh_things as sh
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

data_path = '/home/claire/Works/aspect/runs/model-output/'
case = 'Ra3e8-eta1e8-wide-ascii'
ts0 = 137000
tsf = 137900
R_p = 6371
# d, dT, alpha = 600, 442, 4e-5 # Lees table 1-2: Ra=1e6
# d, dT, alpha = 2700, 3000, 2e-5  # Venus
d, dT, alpha = 2890, 3000, 3e-5  # Hoggard AGU Monograph
# d, dT, alpha = 2700, 3000, 3e-5  # test

fig, ax = sh.dct_spectrum_avg(case, L_x=8,
                              dim=False, R_p=d, d=d, dT=dT, alpha=alpha,
                              ts0=ts0, tsf=tsf, x_res=1, t_res=100,
                              test=False, data_path=data_path, check_norm=False,
                              plot=True, load=True, dump=False, save=False, y0_guide=1e0,
                              #l_min=30, l_max=300,
                              )

path = '/home/claire/Works/exo-top/benchmarks/lees_topo_grids/'
# f = 'psd_freefree.csv'
# P, k = sh.load_spectrum_wavenumber(fpath=path, fname=f, has_header=True, wl=True, two_d=True)
# ax.plot(k, P, label='Lees Ra 1e6')
#
f = 'psd_hoggard.csv'
P, k = sh.load_spectrum_wavenumber(fpath=path, fname=f, has_header=True, wl=True, two_d=True)
ax.plot(k, P, label='Hoggard')

ax.legend(fontsize=16, frameon=False, bbox_to_anchor=(1.05, 1), loc='upper left', )
plt.xlim([1e-4, 2e-2])
plt.ylim([1e-2, 1e8])  # Hoggard: [1e-2, 1e3]
# xticks = ax.get_xticks()
# ax.set_xticks(xticks)
# ax.set_xticklabels([1 / tick for tick in xticks])
# ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
# ax.xaxis.set_minor_formatter(ticker.ScalarFormatter())
# plt.show()
fig.savefig('psd_d'+str(d)+'.png', bbox_inches='tight')
