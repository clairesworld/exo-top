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
d, dT, alpha = 1, 1, 1
# d, dT, alpha = 2890, 3000, 3e-5  # Hoggard AGU Monograph dim factors
R = 2*d

# l, S = sh.make_any_reference(slope=-2, l_rolloff=None, newfname='spectrum', fig_path='', lmin=1, lmax=100,
#                              plot=True, R=R, S0=1, max_dscale=2, d=1)
# case = 'Ra1e8-eta1e7-wide'
# sh.make_model_spectrum(case, R_b=2, data_path=data_path, fig_path='', newfname='base_spectrum',
#                        bl_fudge=2 * np.pi, max_dscale=1, plot=False, verbose=False)
# l, S = sh.load_model_spectrum_pkl(path='', fname='spectrum_-1.pkl')


# test effect of R = N*d
d = 1
CRF = 0.2
R = (1/CRF)*d
print('R = ', 1/CRF, 'd')
l, S = sh.make_any_reference(slope=-2, l_rolloff=5, newfname='spectest', fig_path='', lmin=1, lmax=100,
                             plot=False, R=R, S0=1, max_dscale=2, d=d)

S_scaled = sh.scale_psd_to_rms(phi0=S, k=None, rms1=1, R=R, l=l)
print('rms of scaled', sh.parseval_rms(S_scaled, sh.l_to_k(l, R)))

# get mean peak

tol = 0.00001
mean_old, mean = 0, 100
h_rms = []
h_peak = []
i = 0
while abs(mean_old - mean) > tol:
    mean_old = mean
    clm = sh.random_harms_from_psd(S_scaled, l, R=R, h_ratio=1, plot=False, verbose=False)

    grid = sh.coeffs_to_grid(clm, plot_grid=False, plot_spectrum=False, R=R,
                                      verbose=False, cbar=False, labelsize=labelsize, cmap='nipy_spectral',
                                      # cline='xkcd:off white'
                                      )

    # plt.show()
    h_rms.append(np.sqrt(np.mean(grid.data ** 2)))
    h_peak.append(np.max(abs(grid.data)))
    mean = np.mean(h_rms)
    # print(i, mean)
    i += 1
print('R', R, 'rms', np.mean(h_rms), 'peak', np.mean(h_peak))


# plot example
# fig = plt.gcf()
# ax = plt.gca()
# ax.set_xlabel("Spherical harmonic degree", fontsize=labelsize, labelpad=16)
# ax.set_ylabel("Power (km$^2$ km$^2$)", fontsize=labelsize, labelpad=16)
# ax.set_xticks([2, 5, 10, 20, 50])
# # ax.xaxis.set_minor_formatter(ticker.ScalarFormatter())
# ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
# # fig, ax = dark_background(fig, ax)
# # fig.savefig('psd_example.png', bbox_inches='tight', transparent=True)
# plt.show()

