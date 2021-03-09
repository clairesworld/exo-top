from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, data_path, fig_path, c_rms, c_peak, \
    c_regimes_td, fig_fmt, \
    regime_grid_td, regime_names_td, load_grid, p_Earth, postprocess_kwargs
import fractals as fract

ticksize = 22
axissize = 40
c_fit = 'xkcd:off white'
c_rms = ['xkcd:lime green', 'xkcd:lilac', 'xkcd:orange', 'xkcd:yellow']
lw = 5
ms = 25
elw = 2
ecapsize = 8
# data_path = '/home/claire/Works/aspect/runs/model-output/'

# fract.plot_h_fractal_scaling(case='Ra3e8-eta1e7-wide-ascii', ni=10, rho=rho, alpha=alpha, c_p=c_p, kappa=kappa)

# cases = ['Ra3e8-eta1e6-wide', 'Ra3e8-eta1e7-wide-ascii', 'Ra3e8-eta1e8-wide-ascii',
#              'Ra1e8-eta1e6-wide', 'Ra1e8-eta1e7-wide', 'Ra1e8-eta1e8-wide-ascii']
cases = ['Ra1e8-eta1e8-wide-ascii', 'Ra3e8-eta1e8-wide-ascii']
ts = [[133000, 133900],[137000, 137900]]
for ii, case in enumerate(cases):
    # fract.plot_MHF(case=case, x_res=1, t_res=10)
    fract.dct_spectrum_avg(case, ts0=ts[ii][0], tsf=ts[ii][1], t_res=100, x_res=1, norm='ortho', data_path=data_path,
                           fig_path=fig_path, plot=True, dim=True)
