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
g_E = 9.807  # earth grav
Y = 100e6

labelsize = 32
ticksize = 25
legsize = 20
lw = 5
fig, ax = plt.subplots(1, 1, figsize=(8, 8))

w_Tre = [1e-5, 0.3e-2, 2.9e-2]  # CMF 18, 25, 33%
e_Tre = [3e-6, [[0], [1.8e-2]], [[1.5e-2], [1.7e-2]]]  # wmf error
M_Tre = 0.692

w_Trf = [0.7e-2, 1.9e-2, 4.5e-2]
e_Trf = [w_Trf[0] - 1e-5 + 3e-6, [[1.3e-2], [1.9e-2]], [[1.2e-2], [1.8e-2]]]
M_Trf = 1.039

w_Trg = [0.72e-2, 3.5e-2, 6.4e-2]
e_Trg = [[[0], [1.3e-2]], [[1.3e-2], [1.6e-2]], [[1.6e-2], [2e-2]]]
M_Trg = 1.321

err_c = ['0.1', '0.35', '0.6']
err_m = ['o', '^', 's']
err_kwargs = {'elinewidth': 1, 'capsize': 5, 'ms': 7, 'lw': 0}
labels = ['CMF = 18\%', 'CMF = 25\%', 'CMF = 33\%']
annosize = legsize  # - 5
limkwargs = err_kwargs.copy()
limkwargs.update({'marker': '_'})
off = [0, 0.02, 0.04]  # offset

for i_cmf in range(3):
    c = err_c[i_cmf]
    anno_c = 'k'
    err_kwargs.update({'mec': c, 'mfc': c, 'ecolor': c, 'marker': err_m[i_cmf]})
    limkwargs.update({'mec': c, 'mfc': c, 'ecolor': c, 'marker': err_m[i_cmf]})

    # TRAPPIST-1e
    if i_cmf == 0:
        uplims = True
    else:
        uplims = False
    ax.errorbar(M_Tre + off[i_cmf], w_Tre[i_cmf], yerr=e_Tre[i_cmf], xerr=0.022, uplims=uplims, **err_kwargs)
    if i_cmf == 1:
        ax.errorbar(M_Tre + off[i_cmf], w_Tre[i_cmf], yerr=w_Tre[i_cmf] - 1e-3 + 3e-4, uplims=True,
                    **err_kwargs)  # 25% CMF
    if i_cmf == 2:
        ax.annotate('e', (M_Tre + off[i_cmf] + 0.05, w_Tre[i_cmf] + e_Tre[i_cmf][0][-1] + 0.002), fontsize=annosize,
                    color=anno_c, ha='center', va='bottom')

    # TRAPPIST-1f
    off[i_cmf] = off[i_cmf] * 2
    if i_cmf == 0:
        kw = limkwargs
        uplims = True
    else:
        kw = err_kwargs
        uplims = False
    ax.errorbar(M_Trf + off[i_cmf], w_Trf[i_cmf], yerr=e_Trf[i_cmf], xerr=0.031, uplims=uplims, **kw)
    if i_cmf == 2:
        ax.annotate('f', (M_Trf + off[i_cmf] + 0.07, w_Trf[i_cmf] + e_Trf[i_cmf][0][-1] + 0.005), fontsize=annosize,
                    color=anno_c, ha='center', va='bottom')

    # TRAPPIST-1g
    off[i_cmf] = off[i_cmf] * 1.5
    ax.errorbar(M_Trg + off[i_cmf], w_Trg[i_cmf], yerr=e_Trg[i_cmf], xerr=0.038, label=labels[i_cmf], **err_kwargs)
    if i_cmf == 0:
        ax.errorbar(M_Trg + off[i_cmf], w_Trg[i_cmf], yerr=w_Trg[i_cmf] - 1e-5 + 3e-6, uplims=True, **err_kwargs)
    if i_cmf == 2:
        ax.annotate('g', (M_Trg + off[i_cmf] + 0.15, w_Trg[i_cmf] + e_Trf[i_cmf][0][-1] + 0.01), fontsize=annosize,
                    color=anno_c, ha='center', va='bottom')


def radius_zeng(M_p, CMF=0.3):
    # input mass in kg, output radius in m
    # applicable to M_E <= 8 and CMF <= 0.4
    #     print('using Zeng radius model')
    return (1.07 - 0.21 * CMF) * (M_p / M_E) ** (1 / 3.7) * R_E


def grav(M, R):
    """Calculate acceleration due to gravity on a point mass in m s^-2"""
    return 6.674e-11 * M / R ** 2


def vol(R):
    return 4 / 3 * np.pi * R ** 3


def h_peak_rock(M=None, R=None, Y=100e6, rho_c=2700, **kwargs):
    g = grav(M, R)
    return (2 * Y) / (rho_c * g)


def h_peak_dt(M=None, peak_scale=3.5, **kwargs):  # M in kg
    h_rms = (M / M_E) ** -0.44 * 10 ** 2.26
    return h_rms * peak_scale


def h_peak_dt_cold(M=None, peak_scale=3.5, **kwargs):  # M in kg
    h_rms = (M / M_E) ** -0.3434899 * (10 ** 3.1794847)
    return h_rms * peak_scale


def h_peak_dt_hot(M=None, peak_scale=3.5, **kwargs):  # M in kg
    h_rms = (M / M_E) ** -0.74629234 * (10 ** 2.45476118)
    return h_rms * peak_scale


def wmf_max(M_p, R_p, rho_m=3500, rho_w=1000, h_peak_fn=None, **kwargs):
    # M_p R_p in SI units
    V_p = vol(R_p)
    h = h_peak_fn(M=M_p, R=R_p, **kwargs)  # in m
    # m_cap = rho_w * (V_p - 4/3*np.pi*(R_p - h)**3)
    m_cap = rho_w * rho_m / (rho_m - rho_w) * 4 * np.pi * R_p ** 2 * h  # accounts for water loading
    return m_cap / M_p  # as mass fraction


# give an M, what is the max water content that would make it a waterworld?


masses = np.linspace(0.1, 5, num=10)
gravs = np.ones_like(masses)
wmfs_rock = np.ones_like(masses)
wmfs_dt_cold = np.ones_like(masses)
wmfs_dt_hot = np.ones_like(masses)
for ii, m in enumerate(masses):
    M_p = m * M_E
    R_p = radius_zeng(M_p)  # in km
    gravs[ii] = grav(M_p, R_p) / g_E

    wmf_rock = wmf_max(M_p, R_p, h_peak_fn=h_peak_rock)
    # print(M_p/M_E, 'M_E:', wmf_rock, 'wt %')
    wmfs_rock[ii] = wmf_rock

    wmf_dt_cold = wmf_max(M_p, R_p, h_peak_fn=h_peak_dt_cold)
    # print(M_p/M_E, 'M_E:', wmf_rock, 'wt %')
    wmfs_dt_cold[ii] = wmf_dt_cold

    wmf_dt_hot = wmf_max(M_p, R_p, h_peak_fn=h_peak_dt_hot)
    # print(M_p/M_E, 'M_E:', wmf_rock, 'wt %')
    wmfs_dt_hot[ii] = wmf_dt_hot

pr = np.polyfit(np.log10(masses), np.log10(wmfs_rock), deg=1)
print('rock:', pr)
pd = np.polyfit(np.log10(masses), np.log10(wmfs_dt_cold), deg=1)
print('dt cold:', pd)
pd = np.polyfit(np.log10(masses), np.log10(wmfs_dt_hot), deg=1)
print('dt hot:', pd)


# print(wmfs)

c = '#b13c02ff'
plt.plot(masses, wmfs_rock, c='xkcd:sunshine yellow', lw=lw, label='Crustal strength capacity', zorder=0)
plt.plot(masses, wmfs_dt_cold, c='xkcd:grey green', lw=lw, label='Dynamic topography capacity, cold', zorder=0)
plt.plot(masses, wmfs_dt_hot, c='xkcd:grey green', lw=1, label='Dynamic topography capacity, hot', zorder=0)
plt.xlabel(r'Planet mass ($M_{\oplus}$)', fontsize=labelsize, labelpad=20)
plt.ylabel(r'Surface water mass fraction', fontsize=labelsize, labelpad=20)
plt.loglog()

# fig, ax = dark_background(fig, ax)

handles, labels = ax.get_legend_handles_labels()
# sort both labels and handles by labels
labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
ax.legend(handles, labels, ncol=2, frameon=False, fontsize=legsize,
          bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower center", borderaxespad=0)

# yticks = [1e-4, 1e-1, 1, 10]
xticks = [0.1, 1, 5]
# ax.set_yticks(yticks)
# ax.set_ylim(yticks[0], yticks[-1])
ax.set_xticks(xticks)
ax.set_xlim(xticks[0], xticks[-1])
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
ax.xaxis.set_minor_formatter(ticker.NullFormatter())
ax.tick_params(axis='x', labelsize=ticksize)
ax.tick_params(axis='y', labelsize=ticksize)

plt.savefig(fig_path + 'ww_lim.pdf', bbox_inches='tight',
            # transparent=True
            )
# plt.show()
