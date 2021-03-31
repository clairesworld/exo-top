import pyshtools
# import cartopy.crs as ccrs
import numpy as np
import pandas as pd
from postaspect import aspect_post as ap
from postaspect import aspectdata as post
from postaspect.plt_aspect import plot_save
import matplotlib.pyplot as plt

def hpeak_from_spectrum(power, norm='4pi', lmax=40, n=10, **kwargs):
    ii = 0
    h_peak = []
    while ii < n:
        coeffs_global = pyshtools.SHCoeffs.from_random(power, normalization=norm)
        topo = coeffs_global.expand(lmax=lmax)
        h_peak.append(np.max(topo.data))  # in m
        ii += 1
    return np.mean(h_peak)  # average of a bunch of random spectra consistent with given


def load_spectrum(fpath='', fname='', **kwargs):
    df = pd.read_csv(fpath + fname, header=None, names=['l', 'S_l'], index_col=False)
    degrees = df['l'][1:].to_numpy()
    power = df['S_l'][1:].to_numpy()
    # todo : find nearest l for imperfect digitization
    return power, degrees


def load_spectrum_wavenumber(fpath='', fname='', has_header=True, wl=False, two_d=False, **kwargs):
    names = ['k', 'S']
    if has_header:
        header = 0
    df = pd.read_csv(fpath + fname, header=header, names=names, index_col=False, comment='#')
    k = df['k'][1:].to_numpy()
    power = df['S'][1:].to_numpy()
    if wl:
        k = 1/k
    if two_d:
        power = power*k
    return power, k


def scale_spectrum(h_rms, h_rms0=None, phi0=None, degree=None, pl=None, pl0=None, h_func=None, age=4.5, **kwargs):
    # scale power spectrum phi0 such that it has new rms equal to h_rms

    # if h_func is None:
        # here you are trying to scale dynamic topography spectrum by ratio of h_rms derived from known 1D params
        # get time index nearest to desired snap given in Gyr
        # it = min(enumerate(pl_baseline.t), key=lambda x: abs(age - x[1] * parameters.sec2Gyr))[0]
        # ratio = (pl.d_m * pl.T_m ** 2 * pl.Ra_i_eff ** (-1 / 3)) / (
        #         pl0.d_m * pl0.T_m ** 2 * pl0.Ra_i_eff ** (-1 / 3))  # = phi / phi0
        # ratio = ratio[it]
    if h_rms0 is None:
        h_rms0 = powerspectrum_RMS(power_lm=phi0,degree=degree)  # I hope this is a power spectrum, but this probably needs to be multipled by R or something
    ratio = h_rms/h_rms0
    sqrtphi = np.sqrt(phi0) * ratio
    return sqrtphi ** 2


def vol_from_peak(r0, h_peak, vol_ref=1, **kwargs):
    # getting volume in between trough and peak maxima - as troughs fill in peaks, net volume is 0 and total volume is
    # the spherical annulus between r0 and h_peak
    # R_p = 6051.88e3
    EO = 1.35e9 * 1000**3  # e.g. reference volume based on R_E
    vol = 4 / 3 * np.pi * ((r0 + h_peak) ** 3 - r0 ** 3)
    # print('R_p:', r0*1e-3, 'km, h_peak', h_peak*1e-3, 'km, vol:', vol/EO, 'TO')
    return vol / vol_ref


def vol_from_spectrum(phi0=None, fname='', fpath='', r0=1, n_stats=10, **kwargs):
    # volume between min trough and max peak
    d_ocean_Earth = 3.7e3
    if phi0 is None:
        phi0, l = load_spectrum(fpath=fpath, fname=fname)
    h_peak0 = hpeak_from_spectrum(phi0, l=l, n=n_stats)
    V0 = vol_from_peak(r0, h_peak0)
    return V0


# def powerspectrum_RMS(path=None, power_lm=None, degree=None, amplitude=False, lmax=None): # try to calcuate RMS from digitized power spectrum
#     if path is not None:
#         df = pd.read_csv(path, header=None, names=['degree', 'value'], index_col=False)
#         degree= np.array(df['degree'])
#         S = np.array(df['value'])
#     elif power_lm is not None:
#         degree = np.arange(1, len(power_lm)-1)
#         S = power_lm
#
#     if n is None:
#         n = len(ls)
#
#     RMS_l = []
#     ii = 0
#     while l <= lmax:
#         l = degree[ii]
#         Slm = S[ii]
#         if amplitude:
#             Slm = Slm**2
#         RMS_l.append(np.sqrt(Slm/(2*l + 1)))
#         ii = ii + 1
#     return sum(RMS_l)


def powerspectrum_RMS(path=None, power_lm=None, degree=None, amplitude=False,
                      lmax=None):  # try to calcuate RMS from digitized power spectrum
    if path is not None:
        df = pd.read_csv(path, header=None, names=['degree', 'value'], index_col=False)
        degree = np.array(df['degree'])
        S = np.array(df['value'])
    elif power_lm is not None:
        degree = np.arange(1, len(power_lm) - 1)
        S = power_lm

    if lmax is None:
        lmax = degree[-1]

    RMS_l = []
    ii = 0
    l = degree[ii]

    while (l <= lmax) and ii < len(degree):
        Slm = S[ii]
        if amplitude:
            Slm = Slm ** 2
        RMS_l.append(np.sqrt(Slm / (2 * l + 1)))
        ii = ii + 1
        if ii < len(degree):
            l = degree[ii]

    RMS = sum(RMS_l)
    # print('this RMS', RMS)
    return RMS


def parseval_rms(psd, k):
    # RMS from power spectral density using parseval's theorem
    f = psd
    I = np.trapz(f, k)
    rms = np.sqrt(I/(2*np.pi))
    return rms


def mean_square(L=None, delta_x=None, y=None):
    return (1 / L) * np.sum(y * y * delta_x)


def total_power_parseval(psd=None, delta_k=None):
    return (1.0 / (2.0 * np.pi)) * np.sum(psd * delta_k)


def get_dct_test():
    def f(x):
        # function used for testing routine
        return np.cos(np.pi * x)

    L = 8.0  # length of interval
    N = 128  # number of sample points
    delta_x = L / N  # grid spacing
    x = (np.arange(N) + 0.5) * delta_x  # sample points
    y = f(x)

    return x, y, delta_x, L, N


def dct_spectrum_jfr(case=None, test=False, plot_test=False, L_x=8, x_res=1, data_path='', ts0=None, t0=0.5, dim=False,
                     d=600, dT=442, alpha=4e-5, check_norm=False, **kwargs):
    # 1D power spectral density estimate for Claire
    from scipy.fftpack import dct

    if test:
        x, y, delta_x, L, N = get_dct_test()
    else:
        if ts0 is None:
            ts0 = ap.find_ts(case, t0, data_path=data_path, verbose=False)
        x, y = ap.read_topo_stats(case, ts0, data_path=data_path)
        x = x[::x_res]  # subsample?
        y = y[::x_res]
        N = len(x)
        if dim:
            x = x * d
            y = y * alpha * dT * d
            L = d * L_x
        else:
            L = L_x
        delta_x = L/N

    t = dct(y, type=2, norm='ortho')
    k = (np.pi / L) * np.arange(N)
    delta_k = np.pi / L
    psd = 2 * delta_x * t * t  # 1D power spectral density

    if check_norm:
        ms = mean_square(L, delta_x, y)  # mean square
        print("Mean square of signal", ms)

        # _, rms_prof = ap.peak_and_rms(y)
        # print('Mean square old verj', rms_prof**2)

        total_power = total_power_parseval(psd, delta_k)
        print("Total power (check of Parseval's theorem)", total_power)

        # rms_pars = parseval_rms(psd, k)
        # print('Total power old verj', rms_pars ** 2)

    if plot_test:
        plt.figure()
        plt.plot(x, y, 'bx')
        plt.title("Space domain")
        plt.xlabel("x")
        plt.ylabel("y")

        plt.figure()
        plt.plot(k, psd)
        plt.title("1D Power spectral density")
        plt.xlabel("k")
        plt.ylabel("Power spectral density")

        plt.show()

    return psd, k


# def dct_spectrum_old(case, ts0=None, t0=0.5, x_res=1, norm='ortho', data_path='', plot=False,
#                  L_x=8, dim=False, d=2700, dT=3000, alpha=2e-5, R=6050, test=False, **kwargs):
#     from scipy.fftpack import dct
#
#     if ts0 is None:
#         ts0 = ap.find_ts(case, t0, data_path=data_path, verbose=False)
#
#     x_mids, h = ap.read_topo_stats(case, ts0, data_path=data_path)
#     x_red = x_mids[::x_res]
#     h_red = h[::x_res]
#     if dim:
#         h_red = h_red*alpha*dT*d
#         D_x = d*L_x
#     else:
#         D_x = L_x
#
#     f = dct(h_red, type=2, norm=norm) / 2
#     p = np.arange(len(f))
#     k = np.pi/D_x * p
#     wl = 1/k
#     psd = 2* D_x * f**2
#     # psd_scale = 4*np.pi*R**2 * psd / (2*sh.to_wn(k, R=R) + 1)
#
#     if test:
#         rms_parseval = parseval_rms(psd, k)
#         _, rms_prof = ap.peak_and_rms(h_red)
#         print('frequency rms =', rms_parseval, '| spatial rms =', rms_prof)
#
#     if plot:
#         fig, ax = plot_fit_psd(psd, k, dim=dim, case=case, **kwargs)
#
#     return psd, k


def dct_spectrum_avg(case, ts0=None, tsf=None, t0=None, x_res=1, t_res=100, data_path='', fend='.pkl',
                     plot=False, L_x=8, dim=False, d=2700, dT=3000, alpha=2e-5, load=False, dump=True, **kwargs):
    import os
    import pickle as pkl

    file = data_path+'output-'+case+'/pickle/'+case+'_sph'+fend
    if load and os.path.exists(file):
        psd_mean, k = pkl.load(open(file, "rb"))
    else:
        if ts0 is None:
            ts0 = ap.find_ts(case, t0, data_path=data_path, verbose=False)
        if tsf is None:
            tsf = ap.find_ts(case, 1e9, data_path=data_path, verbose=False)  # use last

        psd_grid = []
        print('Averaging', len(np.arange(ts0, tsf+1, t_res)), 'timesteps')
        for ts in np.arange(ts0, tsf+1, t_res):
            if np.mod(ts, 1000) == 0:
                print('    ts =', ts, '/', tsf)
            psd_i, k = dct_spectrum_jfr(case, ts0=ts, x_res=x_res, data_path=data_path, plot_test=False,
                                        L_x=L_x, dim=dim, d=d, dT=dT, alpha=alpha, **kwargs)
            psd_grid.append(psd_i)

        # take mean
        psd_grid = np.array(psd_grid)
        psd_mean = np.mean(psd_grid, axis=0)

    if dump:
        pkl.dump((psd_mean, k), open(file, "wb"))

    # if dim is False:
    #     psd_mean = psd_mean * d**3 * dT**2 * alpha**2
    #     k = k * d**-1

    if plot:
        fig, ax = plot_fit_psd(psd_mean, k, dim=dim, case=case, alpha=alpha, d=d, data_path=data_path, **kwargs)
        return fig, ax
    else:
        return None, None


def plot_fit_psd(psd, k, dim=True, case='', show_nat_scales=True, save=True, fig_path='', xlim=None, ylim=None,
                 d=2700, dT=3000, alpha=2e-5, l_max=None, l_min=None, labelsize=18, R_p=6050, c_spec='xkcd:slate',
                 x0_guide=7e-4, y0_guide=1e5, x1_guide=2e-3, show_deg=False, show=True, **kwargs):
    import matplotlib.ticker as ticker
    global R

    def to_deg(k):
        return k * 2 * np.pi * R - 0.5

    def to_wn(l):
        return (l + 0.5) / (2 * np.pi * R)
        # return (l + 0.5) / (np.pi * R)  # denominator is planet radius

    if k[0] == 0:  # only wavenumbers greater than 0
        k = k[1:]
        psd = psd[1:]

    plt.figure()
    plt.plot(k, psd, '.-', lw=0.5, alpha=0.5, c=c_spec, label='Power spectral density from DCT-II')
    ax = plt.gca()
    if dim:
        plt.xlabel('$k$, km$^{-1}$', fontsize=labelsize)
        plt.ylabel('$P_k$, km$^3$', fontsize=labelsize)
        R = R_p
    else:
        plt.xlabel('$k$ (nondimensional distance)$^{-1}$', fontsize=labelsize)
        plt.ylabel('$P_k$ (nondimensional distance)$^3$', fontsize=labelsize)
        R = 1
    if show_deg:
        print('k', k)
        print('l', to_deg(k))
        secax = ax.secondary_xaxis('top', functions=(to_deg, to_wn))
        secax.set_xlabel('spherical harmonic degree', fontsize=labelsize)
    # plt.title(case[:12])
    ax.text(0.95, 0.95, case[:12], va='top', ha='right', transform=ax.transAxes, fontsize=labelsize)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.set_yscale('log')
    ax.set_xscale('log')

    if show_nat_scales:
        ax, wl_min, wl_max = nat_scales(case, ax=ax, alpha=alpha, d=d, dim=dim, **kwargs)
    else:
        wl_min, wl_max = nat_scales(case, ax=None, alpha=alpha, d=d, dim=dim, **kwargs)

    if l_min is not None and l_max is not None:
        k_min = to_wn(l_min)
        k_max = to_wn(l_max)
    else:
        k_min = 1/wl_max
        k_max = 1/wl_min
    beta, intercept = fit_slope(psd, k, k_min=k_min, k_max=k_max, ax=ax, fmt='g--', **kwargs)
    show_beta_guide(ax, x0=x0_guide, y0=y0_guide, x1=x1_guide, m=-2, c='k', lw=3, log=True, **kwargs)

    fig = plt.gcf()
    plt.tight_layout()
    if save:
        plot_save(fig, fname='DCT_'+case, fig_path=fig_path, **kwargs)
    elif show:
        plt.show()
    return fig, ax


def fit_slope(S, k, k_min=None, k_max=None, ax=None, fmt='g-', plot=True, **kwargs):
    # find k range
    if k_min is not None and (k_min > np.min(k)):
        i_min = np.argmax(k >= k_min)
    else:
        i_min = 0
    if k_max is not None and (k_max < np.max(k)):
        i_max = np.argmax(k >= k_max) + 1
    else:
        i_max = -1
    kv = k[i_min:i_max]
    Sv = S[i_min:i_max]

    try:
        intercept, slope = np.polynomial.polynomial.polyfit(np.log10(kv), np.log10(Sv), deg=1)
    except TypeError as e:
        print('k', k)
        print('  min, max', k_min, k_max)
        raise e

    intercept = 10**intercept
    print('slope, intercept:', slope, intercept)
    beta = -slope

    if plot:
        if ax is None:
            ax = plt.gca()
        ax.plot(kv, intercept * kv ** -beta, fmt, label=r'naive fit, $\beta$ = ' + '{:.2f}'.format(beta))
        ax.legend()
    return beta, intercept


def show_beta_guide(ax, x0, y0, x1, m=-2, c='xkcd:slate', lw=1, legsize=12, log=True, **kwargs):
    if log:
        b = np.log10(y0) - m*np.log10(x0)
        y1 = 10**b * x1**m
    else:
        b = y0 - m*x0
        y1 = m*x1 + b
    ax.plot((x0, x1), (y0, y1), c=c, lw=lw)
    ax.text((x0+x1)/2, (y0+y1)/2, r'$k^{-2}$', fontsize=legsize, c=c)


def nat_scales(case, ax=None, t1=0, d=2700, alpha=2e-3, c='xkcd:grey', lw=0.5, data_path='', dim=True, **kwargs):

    df = ap.pickleio(case, suffix='_T', t1=t1, load=True, data_path=data_path, **kwargs)
    try:
        T_av, y = ap.time_averaged_profile_from_df(df, 'T_av')
        uv_mag_av, y = ap.time_averaged_profile_from_df(df, 'uv_mag_av')
        dic_av = ap.T_parameters_at_sol(case, n=None, T_av=T_av, uv_mag_av=uv_mag_av, y=y, alpha_m=alpha,
                                        data_path=data_path, **kwargs)
        min_scale = dic_av['delta_rh']
    except KeyError:
        min_scale = np.mean(df.delta_rh.to_numpy())
    if dim:
        min_scale = min_scale * d
        max_scale = 2*d
    else:
        max_scale = 2
    # print('min wl =', min_scale, ', k =', 1/min_scale)
    # print('max wl =', max_scale, ', k =', 1/max_scale)

    if ax is not None:
        ax.axvline(x=1/min_scale, lw=lw, c=c)
        ax.axvline(x=1/max_scale, lw=lw, c=c)
        ylim = ax.get_ylim()
        y_percent = 0.1
        yt = 10**(y_percent*(np.log10(ylim[1]) - np.log10(ylim[0])) + np.log10(ylim[0]))
        # ax.text(1 / max_scale, yt, r'$2d$', va='top', ha='left', fontsize=11, c=c)
        # ax.text(1 / min_scale, yt, r'$\delta_{\rm rh}$', va='top', ha='left', fontsize=11, c=c)

        return ax, min_scale, max_scale
    else:
        return min_scale, max_scale


def k_to_l(k, R):
    return k * 2 * np.pi * R - 0.5


def l_to_k(l, R):
    return (l + 0.5) / (2 * np.pi * R)
