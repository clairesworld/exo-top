# from model_1D import topography as topo
import numpy as np
import matplotlib.pyplot as plt

Ra = np.logspace(6, 11)

h_rms = 0.14 * Ra**-0.18  # fit to chaotic regime
h_peak = 3.19 * Ra**-0.33

plt.loglog(Ra, h_peak, 'kv', label='peak')
plt.loglog(Ra, h_rms, 'r^', label='RMS')
plt.fill_betweenx(plt.gca().get_ylim(), 3e6, 2e7, label='approx ASPECT range')
plt.xlabel('Ra$_{i, eff}$')
plt.ylabel('nondimensional dynamic topography')
plt.legend()
plt.show()
