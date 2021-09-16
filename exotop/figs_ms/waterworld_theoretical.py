import numpy as np
import matplotlib.pyplot as plt
from useful_and_bespoke import dark_background
import matplotlib.ticker as ticker
from matplotlib import rc
from matplotlib.pyplot import rcParams
from datetime import date

today = date.today().strftime("%b-%d-%Y")
fig_path = '/home/claire/Works/exo-top/exotop/figs_ms/'
rc('text', usetex=True)  # turn off for running over ssh
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'CMU Serif'
rho_c = 2700
rho_w = 1000
M_E = 5.972e24  # earth mass in kg
R_E = 6371e3
g_E = 9.807   # earth grav
Y = 100e6

def radius_zeng(M_p, CMF=0.3):
    # input mass in kg, output radius in m
    # applicable to M_E <= 8 and CMF <= 0.4
#     print('using Zeng radius model')
    return (1.07 - 0.21*CMF)*(M_p/M_E)**(1/3.7) * R_E

def grav(M, R):
    """Calculate acceleration due to gravity on a point mass in m s^-2"""
    return 6.674e-11*M/R**2

def vol(R):
    return 4/3*np.pi*R**3

def h_peak_rock(M=None, R=None, Y=100e6, rho_c=2700, **kwargs):
    g = grav(M, R)
    return (2*Y)/(rho_c*g)

def h_peak_dt(M=None, peak_scale=3.5, **kwargs):  # M in kg
    h_rms = (M/M_E)**-0.44 * 10**2.26
    return h_rms*peak_scale

def wmf_max(M_p, R_p, rho_m=3500, rho_w=1000, h_peak_fn=None, **kwargs):
    # M_p R_p in SI units
    V_p = vol(R_p)
    h = h_peak_fn(M=M_p, R=R_p, **kwargs)  # in m
    # m_cap = rho_w * (V_p - 4/3*np.pi*(R_p - h)**3)
    m_cap = rho_w*rho_m/(rho_m - rho_w) * 4 * np.pi * R_p**2 * h  # accounts for water loading
    return m_cap / M_p  # as mass fraction

# give an M, what is the max water content that would make it a waterworld?


masses = np.linspace(0.1, 5, num=10)
gravs = np.ones_like(masses)
wmfs_rock = np.ones_like(masses)
wmfs_dt = np.ones_like(masses)
for ii, m in enumerate(masses):
    M_p = m*M_E
    R_p = radius_zeng(M_p)  # in km
    gravs[ii] = grav(M_p, R_p) / g_E

    wmf_rock = wmf_max(M_p, R_p, h_peak_fn=h_peak_rock) * 1e2
    # print(M_p/M_E, 'M_E:', wmf_rock, 'wt %')
    wmfs_rock[ii] = wmf_rock

    wmf_dt = wmf_max(M_p, R_p, h_peak_fn=h_peak_dt) * 1e2
    # print(M_p/M_E, 'M_E:', wmf_rock, 'wt %')
    wmfs_dt[ii] = wmf_dt

pr = np.polyfit(np.log10(masses), np.log10(wmfs_rock), deg=1)
print('rock:', pr)
pd = np.polyfit(np.log10(masses), np.log10(wmfs_dt), deg=1)
print('dt:', pd)

# print(wmfs)

c = '#b13c02ff'
labelsize = 32
ticksize = 25
legsize = 20
lw = 5
fig, ax = plt.subplots(1, 1, figsize=(8,8))
plt.plot(gravs, wmfs_rock, c='xkcd:saffron', lw=lw, label='Crustal strength limit')
plt.plot(gravs, wmfs_dt, c='xkcd:amethyst', lw=lw, label='Peak dynamic topography limit')
plt.xlabel(r'Surface gravity ($g_{\oplus}$)', fontsize=labelsize, labelpad=20)
plt.ylabel(r'Surface water mass fraction (wt.\%)', fontsize=labelsize, labelpad=20)
plt.loglog()

yticks = [1e-3, 1e-1, 1e0]
xticks = [0.3, 1, 2]
ax.set_yticks(yticks)
ax.set_ylim(yticks[0], yticks[-1])
# ax.set_xticks(xticks)
# ax.set_xlim(xticks[0], xticks[-1])
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
ax.xaxis.set_minor_formatter(ticker.NullFormatter())
ax.tick_params(axis='x', labelsize=ticksize)
ax.tick_params(axis='y', labelsize=ticksize)
# fig, ax = dark_background(fig, ax)

w_Tre = 0.3
M_Tre = 0.692
ax.scatter(M_Tre, w_Tre, s=30, marker='o', c='k')
ax.annotate('TRAPPIST-1e', (M_Tre + 0.05, w_Tre + 0.01), fontsize=legsize)

ax.legend(frameon=False, fontsize=legsize, loc='lower left')

plt.savefig(fig_path + 'ww_lim.png', bbox_inches='tight',
            # transparent=True
            )
plt.show()
