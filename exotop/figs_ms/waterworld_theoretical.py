from __future__ import unicode_literals
import numpy as np
import matplotlib.pyplot as plt
from useful_and_bespoke import dark_background
import matplotlib.ticker as ticker
from datetime import date
from matplotlib import rc
from matplotlib.pyplot import rcParams
from model_1D import parameters as p
from matplotlib import use as mpluse
import matplotlib


# matplotlibuse('Agg')  # turn on for running over ssh
# mpluse('pgf')
rc('text', usetex=True)  # turn off for running over ssh
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'CMU Serif'
# plt.rcParams.update({
#     "pgf.texsystem": "xelatex",
#     "pgf.preamble": "\n".join([
#          r"\usepackage[utf8x]{inputenc}",
#          r"\usepackage[T1]{fontenc}",
#          r"\usepackage{cmbright}",
#     ]),
# })



today = date.today().strftime("%b-%d-%Y")
fig_path = '/home/claire/Works/exo-top/exotop/figs_ms/'
rc('text', usetex=True)  # turn off for running over ssh
# rcParams['font.family'] = 'serif'
# rcParams['font.serif'] = 'CMU Serif'
rho_c = 2700
rho_w = 1000
M_E = 5.972e24  # earth mass in kg
R_E = 6371e3
g_E = 9.807  # earth grav
Y = 100e6

labelsize = 25  ##32
ticksize = 25
legsize = 20
lw = 5
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
rcParams['legend.title_fontsize'] = legsize

def plot_trappist_wmf():
    w_Tre = [1e-5, 0.3e-2, 2.9e-2]  # CMF 18, 25, 33%
    e_Tre = [3e-6, [[0], [1.8e-2]], [[1.5e-2], [1.7e-2]]]  # wmf error
    M_Tre = 0.692

    w_Trf = [0.7e-2, 1.9e-2, 4.5e-2]
    e_Trf = [w_Trf[0] - 1e-5 + 3e-6, [[1.3e-2], [1.9e-2]], [[1.2e-2], [1.8e-2]]]
    M_Trf = 1.039

    w_Trg = [0.72e-2, 3.5e-2, 6.4e-2]
    e_Trg = [[[0], [1.3e-2]], [[1.3e-2], [1.6e-2]], [[1.6e-2], [2e-2]]]
    M_Trg = 1.321

    err_c = ['0.1', '0.35', '0.6']  # ms
    # err_c = ['0.4', '0.65', '0.9']  # slides
    err_m = ['o', '^', 's']
    err_kwargs = {'elinewidth': 1, 'capsize': 5, 'ms': 7, 'lw': 0}
    labels = ['CMF = 18\%', 'CMF = 25\%', 'CMF = 33\%']  # ms
    # labels = ['18\% iron', '25\% iron', '33\% iron']  # slides
    annosize = legsize  # - 5
    limkwargs = err_kwargs.copy()
    limkwargs.update({'marker': '_'})
    off = [0, 0.02, 0.04]  # offset

    for i_cmf in range(3):
        c = err_c[i_cmf]
        anno_c = 'k'  # ms
        # anno_c = 'w'  # slides
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


def h_peak_rock(M=None, R=None, Y=100e6, rho_c=2700, C=1 / 2, **kwargs):
    # C is 1/3 to 1/2 - min stress difference supported elastically underneath load (Jeffreys' theorem)
    g = grav(M, R)
    return (C ** -1 * Y) / (rho_c * g)


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
wmfs_200 = np.ones_like(masses)
wmfs_100 = np.ones_like(masses)
wmfs_50 = np.ones_like(masses)
wmfs_dt_cold = np.ones_like(masses)
wmfs_dt_hot = np.ones_like(masses)
for ii, m in enumerate(masses):
    M_p = m * M_E
    R_p = radius_zeng(M_p)  # in km
    gravs[ii] = grav(M_p, R_p) / g_E

    # wmf_granite = wmf_max(M_p, R_p, h_peak_fn=h_peak_rock, Y=131e6, rho_c=2700)
    # print(M_p/M_E, 'M_E:', wmf_rock, 'wt %')
    wmfs_100[ii] = wmf_max(M_p, R_p, h_peak_fn=h_peak_rock, Y=100e6, rho_c=2700)

    # wmf_basalt = wmf_max(M_p, R_p, h_peak_fn=h_peak_rock, Y=58e6, rho_c=3000)
    wmfs_50[ii] = wmf_max(M_p, R_p, h_peak_fn=h_peak_rock, Y=50e6, rho_c=2700)
    wmfs_200[ii] = wmf_max(M_p, R_p, h_peak_fn=h_peak_rock, Y=200e6, rho_c=2700)

    wmf_dt_cold = wmf_max(M_p, R_p, h_peak_fn=h_peak_dt_cold)
    # print(M_p/M_E, 'M_E:', wmf_rock, 'wt %')
    wmfs_dt_cold[ii] = wmf_dt_cold

    wmf_dt_hot = wmf_max(M_p, R_p, h_peak_fn=h_peak_dt_hot)
    # print(M_p/M_E, 'M_E:', wmf_rock, 'wt %')
    wmfs_dt_hot[ii] = wmf_dt_hot

pr = np.polyfit(np.log10(masses), np.log10(wmfs_100), deg=1)
print('granite:', pr)
# pr = np.polyfit(np.log10(masses), np.log10(wmfs_basalt), deg=1)
# print('basalt:', pr)
pd = np.polyfit(np.log10(masses), np.log10(wmfs_dt_cold), deg=1)
print('dt cold:', pd)
pd = np.polyfit(np.log10(masses), np.log10(wmfs_dt_hot), deg=1)
print('dt hot:', pd)

# Mars
M_Mars = 6.39e23
R_Mars = 3389.5e3
h_peak_Mars_gr = h_peak_rock(M=M_Mars, R=R_Mars, Y=100e6, rho_c=2700)
h_peak_Mars_ba = h_peak_rock(M=M_Mars, R=R_Mars, Y=50e6, rho_c=2700)
h_peak_Mars = 21.9e3  # Olympus Mons
print('\nMars\n---')
print('granite max:', h_peak_Mars_gr * 1e-3, 'km')
print('basalt max:', h_peak_Mars_ba * 1e-3, 'km')
print('observed:', h_peak_Mars * 1e-3, 'km', 'error', abs(h_peak_Mars_gr - h_peak_Mars) / h_peak_Mars)

# Earth
M_E = p.M_E
R_E = p.R_E
h_peak_E_gr = h_peak_rock(M=M_E, R=R_E, Y=131e6, rho_c=2700)
h_peak_E_ba = h_peak_rock(M=M_E, R=R_E, Y=58e6, rho_c=3000)
h_peak_E = 8840  # Mt Everest
print('\nEarth\n---')
print('granite max:', h_peak_E_gr * 1e-3, 'km')
print('basalt max:', h_peak_E_ba * 1e-3, 'km')
print('observed:', h_peak_E * 1e-3, 'km', 'error', abs(h_peak_E_gr - h_peak_E) / h_peak_E)

# Venus
M_V = 4.867e24
R_V = 6051.8e3
h_peak_V_gr = h_peak_rock(M=M_V, R=R_V, Y=131e6, rho_c=2700)  # 131e6
h_peak_V_ba = h_peak_rock(M=M_V, R=R_V, Y=58e6, rho_c=3000)
h_peak_V = 11000  # Maxwell Montes
print('\nVenus\n---')
print('granite max:', h_peak_V_gr * 1e-3, 'km')
print('basalt max:', h_peak_V_ba * 1e-3, 'km')
print('observed:', h_peak_V * 1e-3, 'km', 'error', abs(h_peak_V_gr - h_peak_V) / h_peak_V)

# print(wmfs)


""" DO PLOTTING """
plot_trappist_wmf()
handles, labels = ax.get_legend_handles_labels()
# sort both labels and handles by labels
labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
leg_trapp = ax.legend(handles, labels, ncol=1, frameon=False, fontsize=legsize,
                      title=r'\textbf{TRAPPIST-1 system'+'\n'+r'\textbf{water fraction estimates}',
                      bbox_to_anchor=(1.04, 0.5), #(0, 1.02, 1, 0.2),
                      # loc="lower center",
                      borderaxespad=0)

c = '#b13c02ff'
l1, = plt.plot(masses, wmfs_200, c='xkcd:sand', lw=lw, label='Crust strength, 200 MPa', zorder=0)
l2, = plt.plot(masses, wmfs_100, c='xkcd:sand', lw=3, label='Crust strength, 100 MPa', zorder=0)
l3, = plt.plot(masses, wmfs_50, c='xkcd:sand', lw=1, label='Crust strength, 50 MPa', zorder=0)
# plt.plot(masses, wmfs_basalt, c='xkcd:sand', lw=1, label='Crust capacity, basalt', zorder=0)
l4, = plt.plot(masses, wmfs_dt_cold, c='xkcd:grey green', lw=lw, label='Dynamic topography, cold', zorder=0)
l5, = plt.plot(masses, wmfs_dt_hot, c='xkcd:grey green', lw=1, label='Dynamic topography, hot', zorder=0)

leg_scale = ax.legend([l1, l2, l3, l4, l5],
                      ['Crust strength, 200 MPa', 'Crust strength, 100 MPa', 'Crust strength, 50 MPa',
                       'Dynamic topography, cold', 'Dynamic topography, hot'], ncol=1, frameon=False, fontsize=legsize,
                      title=r'\textbf{Topographic flooding limits}',
                      bbox_to_anchor=(1.04, 1), borderaxespad=0)


c_ss = 'xkcd:cobalt'
p1 = plt.scatter(M_V/M_E, 4.727695769579759e+18*1000/M_V, c=c_ss, marker='*', s=60, zorder=1000)
p2 = plt.scatter(1, 4.641828875736904e+18*1000/M_E, c=c_ss, marker='*', s=60, zorder=1000)
p3 = plt.scatter(M_Mars/M_E, 4.013894774478276e+18*1000/M_Mars, c=c_ss, marker='*', s=60, zorder=1000)
plt.annotate(r'\textbf{Venus}', (M_V/M_E+0.07, 4.727695769579759e+18*1000/M_V+0.0003), c='w', fontsize=legsize,
             bbox=dict(boxstyle="square,pad=0.1", fc=c_ss, ec=None, lw=0.5)
             )  # Venus  u'\u2640', chr(0x263f+1)
plt.annotate(r'\textbf{Earth} (limit)', (1+0.1, 4.641828875736904e+18*1000/M_E-0.00008), c='w', fontsize=legsize, va='top',
             bbox=dict(boxstyle="square,pad=0.1", fc=c_ss, ec=None, lw=0.5)
             )  # Earth chr(0x263f+2)
plt.annotate(r'\textbf{Mars}', (M_Mars/M_E+0.01, 4.013894774478276e+18*1000/M_Mars), c='w', fontsize=legsize,
             bbox=dict(boxstyle="square,pad=0.1", fc=c_ss, ec=None, lw=0.5)
             )  # Mars chr(0x263f+3)
p4 = plt.scatter(1, 1.4e21/M_E, c=c_ss, marker='*', s=60, zorder=1000)
plt.annotate(r'\textbf{Earth}'+'\n(modern ocean)', (1-0.07, 1.4e21/M_E-0.00003), c='w', fontsize=legsize, va='top', ha='right',
             bbox=dict(boxstyle="square,pad=0.1", fc=c_ss, ec=None, lw=0.5))
# leg_ss =

leg_scale._legend_box.align = "left"
leg_trapp._legend_box.align = "left"
ax.add_artist(leg_scale)
# ax.add_artist(leg_ss)
ax.add_artist(leg_trapp)


plt.xlabel(r'Planet mass ($M_{\oplus}$)', fontsize=labelsize, labelpad=20)
plt.ylabel(r'Surface water mass fraction', fontsize=labelsize, labelpad=20)
plt.loglog()

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

# fig, ax = dark_background(fig, ax)

plt.savefig(fig_path + 'ww_lim-' + today + '.pdf', bbox_inches='tight',
            # transparent=True,
            facecolor=fig.get_facecolor()
            )
# plt.show()
