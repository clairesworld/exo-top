from postaspect.setup_postprocessing import Ra_ls, eta_ls, t1_grid, end_grid, load_grid, fig_fmt
from postaspect import plt_aspect as plat
from postaspect.aspect_post import T_parameters_at_sol, get_cases_list, pickleio, h_at_ts, Nu_at_ts  # noqa: E402
from postaspect import ani_aspect as anims
from useful_and_bespoke import cmap_from_ascii
import os
import matplotlib.ticker as ticker
from matplotlib import rc
from matplotlib.pyplot import rcParams
from datetime import date
# from colormath.color_objects import sRGBColor

""" manuscript gridspec """
data_path = '/home/claire/Works/aspect/runs/model-output/'
today = date.today().strftime("%b-%d-%Y")
fig_path = '/home/claire/Works/exo-top/exotop/figs_ms/'
rc('text', usetex=True)  # turn off for running over ssh
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'CMU Serif'
cmap_path = '/home/claire/Works/exo-top/exotop/plot/cmaps/'
cmap_name = 'pm3d20'
# cmap = cmap_from_ascii(cmap_name, path=cmap_path, end='.txt').reversed()
# cmap = cmap_from_ascii('smooth-cool-warm', path=cmap_path, end='.csv', delimiter=',')
cmap = 'coolwarm'

labelsize = 30
ticksize = 20
legsize = 30  #20
wspace = 0.3
hspace = 0.3
c_h = 'xkcd:silver' # 'xkcd:light purple'  #'xkcd:dark pastel green'  # 'xkcd:tealish green'
hlim = (-3e-2, 3e-2)

fig, axes = anims.T_h_gridspec(case='Ra1e8-eta1e8-wide-ascii', data_path=data_path, fig_path=fig_path,
                               labelsize=labelsize, ticksize=ticksize, legsize=legsize, wspace=wspace, hspace=hspace,
                               legtext=r'Ra = $1 \times 10^8$' + '\n' + r'$\Delta \eta = 1 \times 10^8$',
                               cmap=cmap, save=True, c='k', c_h=c_h, hlim=hlim, avg_prof=False,
                               rasterized=True,
                               w=30, h=5.4, nr=6, nc=30, leg_x=0.12, c_bl='xkcd:silver', vmin=0, vmax=1.1,
                               i_n=209, i_ts=183795, t1=0.09, fig_fmt='.pdf')
