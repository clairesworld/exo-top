import sh_things as sh
import matplotlib.pyplot as plt
from useful_and_bespoke import dark_background, cmap_from_ascii, get_continuous_cmap, minmaxnorm, cornertext
# import postaspect.plt_aspect as plat
import model_1D.the_results as results
import model_1D.evolve as evol
import model_1D.parameters as p
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
import numpy as np
from postaspect.plt_aspect import plot_error_contours
from datetime import date

today = date.today().strftime("%b-%d-%Y")
labelsize = 16
fig_path = '/home/claire/Works/exo-top/exotop/figs_scratch/'
data_path = '/home/claire/Works/aspect/runs/model-output/'
# d, dT, alpha = 1, 1, 1
d, dT, alpha = 2890, 3000, 3e-5  # Hoggard AGU Monograph dim factors

case = 'Ra1e8-eta1e7-wide'
sh.make_model_spectrum(case, R=2, data_path=data_path, fig_path='', newfname='base_spectrum',
                       bl_fudge=2 * np.pi, max_dscale=1, plot=False, verbose=False)
l, S = sh.load_model_spectrum_pkl(path='')

# original
h_ratio = 1
clm = sh.random_harms_from_psd(S, l, R=2, h_ratio=h_ratio, plot=False, verbose=False)

h_rms, h_peak = sh.coeffs_to_grid(clm, plot_grid=False, plot_spectrum=True, d=d, alpha_m=alpha, dT=dT, R=2 * d,
                                  verbose=False, cbar=False, labelsize=labelsize, cmap='nipy_spectral',
                                  # cline='xkcd:off white'
                                  )

fig = plt.gcf()
ax = plt.gca()
ax.set_xlabel("Spherical harmonic degree", fontsize=labelsize, labelpad=16)
ax.set_ylabel("Power (km$^2$ km$^2$)", fontsize=labelsize, labelpad=16)
ax.set_xticks([2, 5, 10, 20, 50])
# ax.xaxis.set_minor_formatter(ticker.ScalarFormatter())
ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
# fig, ax = dark_background(fig, ax)
# fig.savefig('psd_example.png', bbox_inches='tight', transparent=True)
plt.show()
