from postaspect.setup_postprocessing import data_path, fig_path
from postaspect import aspect_post as ap
from postaspect import aspectdata as post
from postaspect.plt_aspect import plot_save
import sh_things as sh
from useful_and_bespoke import minmaxnorm
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

def dct_spectrum(case, ts0=None, t0=0.5, x_res=1, norm='ortho', data_path=data_path, plot=False,
                 L_x=8, dim=False, d=2700, dT=3000, alpha=2e-5, R=6050, test=False, **kwargs):
    from scipy import fftpack

    if ts0 is None:
        ts0 = ap.find_ts(case, t0, data_path=data_path, verbose=False)

    x_mids, h = ap.read_topo_stats(case, ts0, data_path=data_path)
    x_red = x_mids.to_numpy()[::x_res]
    h_red = h.to_numpy()[::x_res]
    if dim:
        h_red = h_red*alpha*dT*d
        D_x = d*L_x
    else:
        D_x = L_x

    f = fftpack.dct(h_red, type=2, norm=norm)
    p = np.arange(len(f))
    k = np.pi*p/D_x
    wl = 1/k
    psd = 4*D_x*abs(f)**2
    # psd_scale = 4*np.pi*R**2 * psd / (2*sh.to_wn(k, R=R) + 1)

    if test:
        rms_parseval = sh.parseval_rms(psd, k)
        _, rms_prof = ap.peak_and_rms(h_red)
        print('frequency rms =', rms_parseval, '| spatial rms =', rms_prof)

    if plot:
        fig, ax = plot_fit_psd(psd, k, dim=dim, case=case, **kwargs)

    return psd, k


def plot_fit_psd(psd, k, dim=True, case='', save=True, fig_path=fig_path, **kwargs):
    plt.plot(k, psd, c='xkcd:slate', label='Power spectral density from DCT-II')
    ax = plt.gca()
    if dim:
        plt.xlabel('$k$, km$^{-1}$')
        plt.ylabel('$P_k$, km$^3$')
        secax = ax.secondary_xaxis('top', functions=(sh.to_deg, sh.to_wn))
        secax.set_xlabel('spherical harmonic degree')
    else:
        plt.xlabel('$k$ (nondimensional distance)$^{-1}$')
        plt.ylabel('$P_k$ (nondimensional distance)$^2$')
    plt.title(case[:12])
    plt.xlim([3e-4, 2e-2])
    plt.ylim([1e-6, 1e12])
    ax.set_yscale('log')
    ax.set_xscale('log')
    beta = fit_slope(psd, k, k_min=3e-4, k_max=5e-3, ax=ax, fmt='g--', **kwargs)

    show_beta(ax, x0=4e-4, y0=1e6, x1=7e-4, m=-2, lw=1, c='k', log=True, fontsize=10)

    fig = plt.gcf()
    plt.tight_layout()
    # plt.show()
    if save:
        plot_save(fig, fname='DCT_'+case, fig_path=fig_path, **kwargs)
    return fig, ax


def dct_spectrum_avg(case, ts0=None, tsf=None, t0=None, x_res=1, t_res=100, norm='ortho', data_path=data_path,
                     plot=False, L_x=8, dim=False, d=2700, dT=3000, alpha=2e-5, **kwargs):

    if ts0 is None:
        ts0 = ap.find_ts(case, t0, data_path=data_path, verbose=False)
    if tsf is None:
        tsf = ap.find_ts(case, 1e9, data_path=data_path, verbose=False)  # use last

    psd_grid = []
    for ts in np.arange(ts0, tsf+1, t_res):
        psd_i, k = dct_spectrum(case, ts0=ts, x_res=x_res, norm=norm, data_path=data_path, plot=False,
                                  L_x=L_x, dim=dim, d=d, dT=dT, alpha=alpha)
        psd_grid.append(psd_i)

    # take mean
    psd_grid = np.array(psd_grid)
    psd_mean = np.mean(psd_grid, axis=0)

    if plot:
        fig, ax = plot_fit_psd(psd_mean, k, dim=dim, case=case, **kwargs)
        return fig, ax
    else:
        return None, None


def fit_slope(S, k, k_min=None, k_max=None, ax=None, i_min=0, i_max=-1, fmt='g-', **kwargs):
    # find k range
    if k_min is not None and (k_min > np.min(k)):
        i_min = np.argmax(k >= k_min)
    if k_max is not None and (k_max < np.max(k)):
        i_max = np.argmax(k > k_max)
    kv = k[i_min:i_max]
    Sv = S[i_min:i_max]

    intercept, slope = np.polynomial.polynomial.polyfit(np.log10(kv), np.log10(Sv), deg=1)
    intercept = 10**intercept
    print('slope, intercept:', slope, intercept)
    beta = -slope

    if ax is None:
        ax = plt.gca()
    ax.plot(kv, intercept * kv ** -beta, fmt, label=r'naive fit, $\beta$ = ' + '{:.2f}'.format(beta))
    ax.legend()
    return beta


def show_beta(ax, x0, y0, x1, m=-2, c='xkcd:slate', lw=1, fontsize=12, log=True):
    if log:
        b = np.log10(y0) - m*np.log10(x0)
        y1 = 10**b * x1**m
    else:
        b = y0 - m*x0
        y1 = m*x1 + b
    ax.plot((x0, x1), (y0, y1), c=c, lw=lw)
    ax.text((x0+x1)/2, y0, r'$k^{-2}$', fontsize=fontsize, c=c)


def haarFWT(signal, level=1):
    # https: // stackoverflow.com / questions / 57439509 / implementing - haar - wavelet - in -python - without - packages

    s = .5  # scaling -- try 1 or ( .5 ** .5 )

    h = [1, 1]  # lowpass filter
    g = [1, -1]  # highpass filter
    f = len(h)  # length of the filter

    t = signal;  # 'workspace' array
    l = len(t);  # length of the current signal
    y = [0] * l;  # initialise output

    t = t + [0, 0];  # padding for the workspace

    for i in range(level):

        y[0:l] = [0] * l;  # initialise the next level
        l2 = l // 2;  # half approximation, half detail

        for j in range(l2):
            for k in range(f):
                y[j] += t[2 * j + k] * h[k] * s;
                y[j + l2] += t[2 * j + k] * g[k] * s;

        l = l2;  # continue with the approximation
        t[0:l] = y[0:l];

    return y


def MHF(h, d):
    # mean haar fluctuations along a transect as function of length scale 2n
    # vector d is the distance in m from d0

    L_max = len(d)
    if np.mod(L_max, 2) > 0:
        L_max = L_max - 1  # must be even
    #     print('L_max', L_max) # length of greatest L scale 2n

    MHF_L = []
    L = []
    for twon in range(L_max, 0, -2):
        #         print('\n2n', twon)
        HF = []  # initialise Haar fluctuations at this L
        n = int(twon / 2)
        r = n  # starting right edge of x1 such that r - n >= 0
        while r + 2 + n <= len(d):  # move wavelet along until right edge reaches end
            #             print('r', r)
            #             print('x1 [l:r]', r - n, ':', r + 1)
            #             print('x2 [l:r]', r + 1, ':', r + 2 + n)

            x1 = h[r - n:r + 1]  # first n points
            M1 = np.mean(x1)

            x2 = h[r + 1:r + 2 + n]  # last n points
            M2 = np.mean(x2)

            #             print('x1', len(x1), 'x2', len(x2))
            HF.append(abs(M2 - M1))
            r = r + 1
        mean_HF = np.mean(HF)
        if not np.isnan(mean_HF):
            MHF_L.append(np.mean(HF))
            L.append(twon)
    return MHF_L, L  # mean haar fluctuation with each corresponding length scale


def MHF_from_latlon(lat, lon, h, R, lon_res=20, lat_res=1):
    # planetary mean Haar fluctuations

    n_samples = len(lat)

    # convert latitude to metres from lat[0]
    arc_length = 2 * np.pi * R / 2
    d = minmaxnorm(lat, a=0, b=arc_length)

    # subsample distance vector
    d = d[::lat_res]
    h = h[::lat_res, :]

    haars = []
    Ls = []
    for i in range(0, len(lon), lon_res):
        # take vertical N-S transects after Landais 2019
        hi = h[:, i]

        # calculate mean Haar fluctuations at each L along transect
        MHF_L, L = MHF(hi, d)
        haars.extend(MHF_L)
        Ls.extend(L)
    #         if i == 0:
    #             plt.figure()
    #             plt.plot(L, MHF_L)
    #             plt.ylabel('mean haar fluctuation')
    #             plt.xlabel('distance scale')

    #             plt.figure()
    #             plt.plot(d, hi)
    #             plt.ylabel('altitude')
    #             plt.xlabel('distance from S pole')

    # sort and group
    df = pd.DataFrame({'L': Ls, 'MHF_L': haars})
    df = df.sort_values(by=['L'])
    df = df.groupby(['L']).mean()
    dL = d[1] - d[0]
    return df.index.to_numpy() * dL, df.MHF_L.to_numpy()


def MHF_profiles(case, n_start=None, n_end=None, t_res=20, x_res=1, data_path=data_path, **kwargs):
    # mean haar fluctuation from dyn top profiles
    if os.path.exists(data_path + 'output-' + case):
        dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', read_statistics=False, **kwargs)
        if n_end is None:
            n_end = dat.final_step()
        if n_start is None:
            n_start = dat.final_step() - 10
        ts0 = dat.find_time_at_sol(n_start, return_indices=True)
        ts1 = dat.find_time_at_sol(n_end, return_indices=True)
        times = np.arange(ts0, ts1 + 1, t_res)

        # dat.read_mesh(n_end)
        # print('original x mesh', np.shape(dat.x))
        # x_mesh = dat.x[::x_res]

        x_mids, _ = ap.read_topo_stats(case, ts0, data_path=data_path)
        x_mids = np.array(x_mids[::x_res])

        # load profiles into grid with shape (n_times, n_meshx)
        grid = np.zeros((len(x_mids), len(times)))
        print('grid', np.shape(grid))

        for ii, ts in enumerate(times):
            _, h = ap.read_topo_stats(case, ts, data_path=data_path)
            # print('original h', np.shape(h))
            h = np.array(h[::x_res])
            # print('h', np.shape(h))
            grid[:, ii] = h

        haars = []
        Ls = []
        for ii, ts in enumerate(times):
            # take time slice
            hi = grid[:, ii]

            # calculate mean Haar fluctuations at each L along transect
            MHF_L, L = MHF(hi, x_mids)
            haars.extend(MHF_L)
            Ls.extend(L)

        # sort and group
        df = pd.DataFrame({'L': Ls, 'MHF_L': haars})
        df = df.sort_values(by=['L'])
        df = df.groupby(['L']).mean()
        dL = x_mids[1] - x_mids[0]
        df.head(10)
        return df.index.to_numpy() * dL, df.MHF_L.to_numpy()
    else:
        print('No Aspect Data found:', case)


def plot_MHF(case, fname='MHF_', L_max=None, save=True, **kwargs):
    L, MHF_L = MHF_profiles(case, **kwargs)

    fig, ax = plt.subplots()
    ax.plot(L, MHF_L)
    ax.set_xlabel('$\Delta x$')
    ax.set_ylabel('Mean Haar Fluctuation')
    ax.loglog()

    beta = fit_slope(MHF_L, L, k_max=L_max, fig=fig, ax=ax, i_min=0)
    # slope should = H (because mean)
    print('H =', -beta)
    if save:
        plot_save(fig, fname+case, **kwargs)


# def plot_h_fractal_scaling(case, n=None, rho=3500, kappa=None, k=4, c_p=1200, alpha=3e-5, data_path=data_path, fig_path=fig_path,
#                            figsize=(7, 7), labelsize=16, c='k', lw=3, ni=10, **kwargs):
#
#     if kappa is None and k is not None:
#         kappa = k / (rho * c_p)
#
#     dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', read_statistics=False, **kwargs)
#     if n is None:
#         n = dat.final_step()
#     ts = dat.find_time_at_sol(n)
#
#     x_mids, h = ap.read_topo_stats(case, ts, data_path=data_path)
#     x_mids = np.array(x_mids)
#     h = np.array(h)
#
#     x, y, _, F = dat.read_vertical_heatflux(n)
#     F = post.reduce_dims(F)
#     F_surf = F[:, -1]
#     print('F_surf', np.shape(F_surf))
#     print('h', np.shape(h))
#     print('x', np.shape(x))
#     # todo: F and h are at different x points
#
#     h_t = haar(h, ni)  # apply haar wavelet
#     x_l = x
#     dh_l = np.zeros_like(x_l)
#     phi_l = np.zeros_like(x_l)
#     for i, xi in enumerate(x_l):
#         dh_l[i] = h_t[i] - h_t[0]  # change in topography at that distance scale
#         F_l = np.sum(F_surf[0:i])
#         phi_l[i] = kappa ** 3 * c_p * alpha * rho ** 2 * F_l ** -1
#
#     fig, ax = plt.subplots(figsize=figsize)
#     ax.plot(x_l, phi_l * x_l ** 0.5, c=c, lw=lw, label='predicted scaling')
#     ax.plot(x_l, dh_l, c=c, lw=lw, ls='--', label='ASPECT scaling')
#     ax.set_ylabel(r'$\Delta h', fontsize=labelsize)
#     ax.set_xlabel('L', fontsize=labelsize)
#     fig.savefig(fig_path + 'fractal-scaling.png', bbox_inches='tight')


