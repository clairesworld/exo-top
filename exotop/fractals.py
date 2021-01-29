import sys
sys.path.insert(0, '/home/cmg76/Works/exo-top/')
from exotop.postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    fig_fmt, regime_grid_td, postprocess_kwargs, regime_names_td, load_grid, p_Earth    # noqa: E402
from exotop.postaspect import aspect_scalings as sc  # noqa: E402
from exotop.postaspect import aspect_postprocessing2 as post  # noqa: E402
from exotop.useful_and_bespoke import dark_background, cmap_from_list
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.ticker as ticker

ticksize = 22
axissize = 40
c_fit = 'xkcd:off white'
c_rms = ['xkcd:lime green', 'xkcd:lilac', 'xkcd:orange', 'xkcd:yellow']
lw = 5
ms = 25
elw = 2
ecapsize = 8

def haarFWT ( signal, level=1 ):
    # https: // stackoverflow.com / questions / 57439509 / implementing - haar - wavelet - in -python - without - packages

    s = .5;                  # scaling -- try 1 or ( .5 ** .5 )

    h = [ 1,  1 ];           # lowpass filter
    g = [ 1, -1 ];           # highpass filter
    f = len ( h );           # length of the filter

    t = signal;              # 'workspace' array
    l = len ( t );           # length of the current signal
    y = [0] * l;             # initialise output

    t = t + [ 0, 0 ];        # padding for the workspace

    for i in range ( level ):

        y [ 0:l ] = [0] * l; # initialise the next level
        l2 = l // 2;         # half approximation, half detail

        for j in range ( l2 ):
            for k in range ( f ):
                y [j]    += t [ 2*j + k ] * h [ k ] * s;
                y [j+l2] += t [ 2*j + k ] * g [ k ] * s;

        l = l2;              # continue with the approximation
        t [ 0:l ] = y [ 0:l ] ;

    return y


def haar(y, ni):
    # ni is number of indices in mean
    n = len(y)
    yp = np.zeros_like(y)
    for i in range(ni, n - ni):
        M1 = np.mean(y[i - ni:i])
        M2 = np.mean(y[i:i + ni])
        yp[i - ni:i + ni] = [abs(M2 - M1)]*2*ni
    return yp


def plot_h_fractal_scaling(case, n=None, rho=1, kappa=1, c_p=1, alpha=1, data_path=data_path, fig_path=fig_path,
                           figsize=(7,7), labelsize=16, c='k', lw=3, ni=10, **kwargs):
    ts = -1
    x_mids, h = sc.read_topo_stats(case, ts, data_path=data_path)
    x_mids = np.array(x_mids)
    h = np.array(h)

    dat = post.Aspect_Data(directory=data_path + 'output-' + case + '/', read_statistics=False, **kwargs)
    n = dat.final_step()
    x, y, _, F = dat.read_vertical_heatflux(n)
    F = post.reduce_dims(F)
    F_surf = F[:,-1]
    print('F_surf', np.shape(F_surf))
    print('h', np.shape(h))
    print('x', np.shape(x))
    # todo: F and h are at different x points

    h_t = haar(h, ni)  # apply haar wavelet
    x_l = x
    dh_l = np.zeros_like(x_l)
    phi_l = np.zeros_like(x_l)
    for i, xi in enumerate(x_l):
        dh_l[i] = h_t[i] - h_t[0]  # change in topography at that distance scale
        F_l = np.sum(F_surf[0:i])
        phi_l[i] = kappa**3 * c_p * alpha * rho**2 * F_l**-1

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(x_l, phi_l*x_l**0.5, c=c, lw=lw, label='predicted scaling')
    ax.plot(x_l, dh_l, c=c, lw=lw , ls='--', label='ASPECT scaling')
    ax.set_ylabel(r'$\Delta h', fontsize=labelsize)
    ax.set_xlabel('L', fontsize=labelsize)
    fig.savefig(fig_path+'fractal-scaling.png', bbox_inches='tight')

rho = 3500
alpha = 3e-5
c_p = 1200
k = 4
kappa = k / (rho * c_p)
plot_h_fractal_scaling(case='Ra3e8-eta1e7-wide-ascii', ni=10, rho=rho, alpha=alpha, c_p=c_p, kappa=kappa)