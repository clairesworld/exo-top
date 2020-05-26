""" mostly plotting results for rocky planet evolution + topographies """

import numpy as np
import parameters as p
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from collections.abc import Iterable
import collections
import six
from scipy import interpolate
import pandas as pd
from IPython.display import display, Math
import random as rand
import terrestrialplanet as tp
import thermal as therm
import topography as top
import rheology as rh
import inputs as ins
from mpl_tools import colorize
import matplotlib
from matplotlib.font_manager import FontProperties
from matplotlib import rc

import pyshtools
import cartopy.crs as ccrs
# matplotlib.font_manager._rebuild()
rc('font',**{'family':'serif','serif':['CMU Serif Roman']}) # latex font for matplotlib
rc('text', usetex=True)
# np.seterr('raise')







###### BUILD PLANETS IN BULK #####

def bulk_planets(N=1, name=None, mini=None, maxi=None, like=None, visc_type='Thi', t_eval=None, random=False,
                 T_m0=1750, T_c0=2250, D_l0=600e3, tf=4.5, 
                 update_kwargs=None, **kwargs): 
    """varying single parameter 'name' between mini and maxi, use default values otherwise."""
    
    if like is not None:
        pl_kwargs = eval('ins.'+like+'_in') 
        model_kwargs =  eval('ins.'+like+'_run')
    else:
        pl_kwargs = {}
        model_kwargs = dict(T_m0=T_m0, T_c0=T_c0, D_l0=D_l0, tf=tf, visc_type=visc_type) # model params
        
    if update_kwargs is not None:
        pl_kwargs.update(update_kwargs)
    
    planets = []
    ii=0
    arr = np.linspace(mini, maxi, num=N)
    while ii<N:
        if random:
            val = rand.uniform(mini, maxi)
        else:
            val = arr[ii]
        new_kwargs = pl_kwargs.copy()
        new_kwargs.update({name:val})
        pl = tp.TerrestrialPlanet(**new_kwargs)
        
        pl = therm.solve(pl, t_eval=t_eval, **model_kwargs) # T_m, T_c, D_l
        pl = top.topography(pl, C=1)
        planets.append(pl)
        ii+=1
    return planets

def build_planet(ident='Earthbaseline', run_args=None, update_args=None):
    planet_kwargs = eval('ins.'+ident+'_in') 
    model_kwargs = eval('ins.'+ident+'_run')
    if run_args is not None:
        model_kwargs.update(run_args)
    if update_args is not None:
        planet_kwargs.update(update_args)
    pl = tp.TerrestrialPlanet(**planet_kwargs)
    pl = therm.solve(pl, **model_kwargs) # T_m, T_c, D_l
    pl = top.topography(pl, C=2)
    return pl

def build_solarsystem(run_args=None, ident_list=['Moon1', 'Mercury1', 'Mars1', 'Venus', 'Earth'], dicts=False):
    planets = []
    for ii, ident in enumerate(ident_list):
        pl = build_planet(ident, run_args)
        planets.append(pl)
    if dicts:
        # Create a zip object from two lists
        z = zip(ident_list, planets)
        # Create a dictionary from zip object
        return dict(z)
    return planets


###### PLOTTING ######

def plot_output(pl, names, ncols=6, tspan=None, title=None, plots_save=False, verbose=False,
                compare_dir=None, fig_path='figs/', labelpad=None,labelsize=15, legsize=10, fname=None,
                line_args=None, cmp_line_args=None, annotate_colour='xkcd:bright purple',
                print_tf=False, colorbar=False, legend=True, hidex=False, fformat='.pdf',
                ident=None, fig=None, axes=None, label=None, cmp_label=None,  ticksize=12,
                fontname=None, suptitlepad=1.04, legax=0, **kwargs):
     # names: y param
    if ident is None:
        ident = pl.ident
    if fname is None:
        fname = ident
    if line_args is None:
        line_args = {'lw':2, 'ls':'-', 'c':'k', 'marker':None, 'ms':5}
    if cmp_line_args is None:
        cmp_line_args = {'lw':1, 'ls':'-', 'c':'r', 'marker':None, 'ms':5}
    
    t = pl.t # time points of ode solutions in s
    
    nrows = int(np.ceil(len(names)/ncols))
    if (fig is None) and (axes is None):
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*3))
    if tspan is None:
        tspan = (0, t[-1]*1e-9/p.years2sec)
    out_vars = list(names.keys())
    ylabels = list(names.values()) # tuple (ylabel, yscale)
    if label is None:
        label = 'this work'
    for n, par in enumerate(out_vars):
        #loop across axes: each y axis variable
#         print('n', n, 'par', par)
        ax = axes.flatten()[n]
        y = eval('pl.'+par)
        if np.size(y)==1:
            y = [y]*len(t)
        try:
#             print('y', y)
            # if the name of the parameter you want to plot exists
            yl = str(ylabels[n][0])
            if par=='eta_m': # always log scale for viscosity
                y = np.log10(y)
            plot_one(ax, t*1e-9/p.years2sec, y*ylabels[n][1], xlabel='', ylabel=yl, ticksize=ticksize, labelpad=labelpad,
                     label=label, fontname=fontname, labelsize=labelsize, legsize=legsize, line_args=line_args)
            if compare_dir is not None: 
                # if data exists to benchmark this param
                try:
                    if (isinstance(compare_dir, collections.Iterable)) and (not isinstance(compare_dir, six.string_types)):
                        for cc, cdir in enumerate(compare_dir):
                            df = pd.read_csv(cdir+'/'+par+'.csv', header=None, names=['time', 'value'],
                                            index_col=False)
                            if cmp_label is None:
                                cmp_label=cdir
                            plot_one(ax, df['time'], df['value'], 
                                     '', yl, labelsize=labelsize, legsize=legsize,  ticksize=ticksize,
                                     label=cmp_label[cc], fontname=fontname, line_args=cmp_line_args[cc])
                    else:
                        # not iterable
                        df = pd.read_csv(compare_dir+'/'+par+'.csv', header=None, names=['time', 'value'],
                                            index_col=False)
                        if cmp_label is None:
                            cmp_label=compare_dir
                        plot_one(ax, df['time'], df['value'], 
                                 '', yl, labelsize=labelsize, legsize=legsize,  ticksize=ticksize,
                                 label=cmp_label, fontname=fontname, line_args=cmp_line_args)
                except IOError:
                    print('file', str(par+'.csv'), 'not found')
                    pass
            if par=='urey' and print_tf: # print final value of urey ratio
                ii = np.where(t*1e-9/p.years2sec<=tspan[-1])
                ax.annotate('%.2f'%(y[ii][-1]), xy=(tspan[-1], y[ii][-1]), fontsize=legsize, 
                            color=annotate_colour,
                            textcoords="axes fraction", xytext=(0.95, 0.2),
                            ha='right', va='bottom',
                            arrowprops=dict(arrowstyle='->', connectionstyle="arc3,rad=-0.1",
                                           ec=annotate_colour))
            ax.set_xlim(tspan)
            if legend and (n==legax):
                ax.legend(frameon=False, fontsize=legsize)
        except ValueError as e:
            print('could\'t plot', par)
            print(e)
    
    while n+1 < ncols*nrows :
        fig.delaxes(axes.flatten()[n+1])
        n += 1 # hide unused axes
    

    plot_setxlabel(axes, 'Age (Gyr)', 'bottom', fontname=fontname, labelpad=labelpad, labelsize=labelsize)
    if title is None:
        title = pl.ident
    fig.suptitle(title, fontsize=labelsize, y=suptitlepad, fontname=fontname)
    
    if colorbar:
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        fig.colorbar(sc, cax=cbar_ax)
        sc.set_visible=False
    
    
    plt.tight_layout()
    if plots_save:
        plt.savefig(fig_path+fname+fformat, bbox_inches='tight') 
    if verbose:
        print('\n        n timesteps =', len(t))
        print('$t_f$ =', t[-1]*p.sec2Gyr, 'Gyr')
        print(r'$R_p$ =', '%.2f'%(pl.R_p/p.R_E), 'R_E =', '%.2f'%(pl.R_p*1e-3), 'km')
#         print(r'$R_c$ =', '%.2f'%(kwargs['CRF']*kwargs['R_p']*1e-3), 'km')
        print('M_c', '%.2e'%pl.M_c, 'kg')   
        print(r'$M_{m+lid}$ =', '%.2e'%(pl.M_m), 'kg')
        print(r'$g_{sfc}$ =', '%.2f'%pl.g_sfc, 'm s^-2')
        print(r'$g_{cmb}$ =', '%.2f'%pl.g_cmb, 'm s^-2')
        print(r'$\kappa_m$', '%.6f'%pl.kappa_m, 'm^2 s^-1')
        print(r'CRF =', '%.2f'%pl.CRF)
        print(r'$h_0$ =', '%.2f'%(pl.h_rad_m[0]*1e12), 'pW kg^-1')
        print(r'$h_{4.5}$ =', '%.2f'%(pl.h_rad_m[-1]*1e12), 'pW kg^-1')
#         print(r'$H_0$ =', '%.2f'%(H_rad_m[0] + H_rad_lid[0]), 'TW')
#         print(r'$H_{4.5}$ =', '%.2f'%(H_rad_m[-1] + H_rad_lid[-1]), 'TW')
        print(r'Urey ratio @ $t_f$ =', '%.2f'%pl.urey[-1])
        print('q_sfc(t=0)', '%.2f'%(pl.q_sfc[0]*1e3), 'mW m^-3')
    return fig, axes

def snaps(pl, plot_snapshots=None, fig_path=None, plots_save=False, ident=None, **kwargs):
    if ident is None:
        ident = pl.ident
    
    t = pl.t # time points of ode solutions in s
    try:
        n_col = len(plot_snapshots)
    except:
        n_col = 1
    fig2, axes2 = plt.subplots(1, n_col,figsize=(3*n_col,5))
    for iax, tx in enumerate(plot_snapshots): # tx is the time value u want nearest
        ii = min(enumerate(t), key=lambda x: abs(tx - x[1]*p.sec2Gyr))[0]
        plot_structure(ax=axes2[iax], t=t[ii], T_m=pl.T_m[ii], T_c=pl.T_c[ii], T_s=pl.T_s,
                       T_l=pl.T_l[ii], R_l=pl.R_l[ii], R_p=pl.R_p, R_c=pl.R_c, h_rad_m=pl.h_rad_m[ii],
                       d_lbl = pl.TBL_c[ii], d_ubl = pl.TBL_u[ii], q_ubl = pl.q_ubl[ii], a0=pl.a0[ii],
                       k_m=pl.k_m, legsize=10, **kwargs)
    plt.tight_layout()
    if plots_save:
        fig2.savefig(fig_path+pl.ident+'_profiles.pdf', bbox_inches='tight') 
    return fig2, axes2

def plot_structure(ax=None, t=None, T_m=None, T_c=None, R_p=None, R_l=None, R_c=None, T_l=None, 
                   T_s=None, h_rad_m=None, d_lbl=None, d_ubl=None, q_ubl=None, a0=None, k_m=None,
                   labelsize=16, legsize=14, Tlid_ini=None, **kwargs):
    """ plot temp structure (for a given time) """
    r_c = np.linspace(0, R_c*1e-3)
    r_lbl = np.linspace(R_c*1e-3, (R_c+d_lbl)*(1e-3))
    r_m = np.linspace((R_c+d_lbl)*1e-3, (R_l-d_ubl)*1e-3) # radius for mantle in km
    r_ubl = np.linspace((R_l-d_ubl)*1e-3, (R_l)*1e-3)
    r_l = np.linspace(R_l*1e-3, R_p*1e-3) # radius for lid
    T_cond = therm.sph_conduction(r_l*1e3, a0=a0, T_l=T_l, R_p=R_p, R_l=R_l, T_s=T_s, k_m=k_m, **kwargs)
    q = therm.sph_flux(r_l*1e3, a0=a0, T_l=T_l, T_s=T_s, R_p=R_p, R_l=R_l, k_m=k_m, **kwargs)
    if Tlid_ini=='linear':
        T_cond = therm.sph_conduction(r_l*1e3, a0=0, T_l=T_l, R_p=R_p, R_l=R_l, T_s=T_s, k_m=k_m,**kwargs)
        q = therm.sph_flux(r_l*1e3, a0=0, T_l=T_l, T_s=T_s, R_p=R_p, R_l=R_l, k_m=k_m,**kwargs)
    
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(3,5))
    ax.axhline(y=R_l*1e-3, ls='--', lw=1, c='xkcd:bluish purple')
    ax.axhline(y=R_c*1e-3, ls='--', lw=1, c='xkcd:clay')

    ax.plot(T_cond, r_l, c='xkcd:bluish purple')
    ax.plot(T_l + q_ubl/k_m*(R_l - r_ubl*1e3), r_ubl, c='xkcd:greenish')
    ax.plot([T_m]*len(r_m), r_m, c='xkcd:greenish')
    ax.plot([T_m, T_c], [(R_c+d_lbl)*1e-3, R_c*1e-3], c='xkcd:greenish', marker = None)
    ax.plot([T_c]*len(r_c), r_c, c='xkcd:clay')
        
    ax.set_ylabel('Radius (km)', fontsize=labelsize)
    ax.set_xlabel('Temperature (K)', fontsize=labelsize)
    ax.set_ylim([0, R_p*1e-3])
#     #ax.fill_between(x=x, y0=[0]*len(x), y1=[R_cmb*1e-3]*len(x), color='xkcd:gold') # core
    ax.text(T_cond[-1], 0, 'Core', ha='left', va='bottom', fontsize=legsize, c='xkcd:clay')
#     #ax.fill_between(x=x, y0=[R_cmb*1e-3]*len(x), y1=[R_p*1e-3]*len(x), color='xkcd:tomato') # mantle
    ax.text(T_cond[-1], R_c*1e-3, 'Convecting region', ha='left', va='bottom', fontsize=legsize, c='xkcd:greenish')
    ax.text(T_cond[-1], R_l*1e-3, 'Lid', ha='left', va='bottom', fontsize=legsize, c='xkcd:bluish purple')
    
    ax2 = ax.twiny()
    ax2.set_xlabel('Flux, steady-state (mW m$^{-2}$)', color='xkcd:grey')  
    ax2.plot(q*1e3, r_l, color='xkcd:grey')
    ax2.plot(q_ubl*1e3, r_ubl[0], marker='o', color='xkcd:grey')
    ax2.annotate('$q_{ubl}$', (q_ubl*1e3, r_ubl[-1]), color='xkcd:grey', fontsize=12, ha="left", va="top")
    ax2.tick_params(axis='x', labelcolor='xkcd:grey')
    
    ax.set_title(('%.1f'%(t*1e-9/p.years2sec))+' Gyr', fontsize=labelsize)

    return ax

def interp_benchmark(path, yscale=1):
    df = pd.read_csv(path, header=None, names=['time', 'value'], index_col=False) 
    f = interpolate.interp1d(np.array(df['time']), np.array(df['value'])*yscale, kind='linear')
    times = df['time'] # in Gyr 
    return np.array(times), f

def plot_qsfc_error(pl, ax3=None, compare_dir=None, fig_path=None, plots_save=False, ident=None, **kwargs):
    """ sanity check on q_sfc """
    if ident is None:
        ident = pl.ident
    
    t = pl.t # time points of ode solutions in s
    if ax3 is None:
        fig3, ax3 = plt.subplots(1, 1, figsize=(5,5))
        
    t_D_l, f_D_l_interp = interp_benchmark(path=compare_dir+'/D_l.csv', yscale=1e3) # in Gyr, m
    temp = t*1e-9/p.years2sec # in Gyr
    
    try:
        t_T_l, f_T_l_interp = interp_benchmark(path=compare_dir+'/T_l.csv')# in Gyr, K
        iii = np.where((temp>=t_T_l.min()) & (temp<=t_T_l.max()))
        times0 = temp[iii] # time points of ODE solver subset to interpolation range
        T_l_interp = f_T_l_interp(times0) # published plot interpolated to model times, in K
    except FileNotFoundError as e:
        t_T_l, f_T_avg_interp = interp_benchmark(path=compare_dir+'/T_avg.csv') # in Gyr, K
        iii = np.where((temp>=t_T_l.min()) & (temp<=t_T_l.max()))
        times0 = temp[iii] # time points of ODE solver subset to interpolation range
        T_avg_interp = f_T_avg_interp(times0) # published plot interpolated to model times, in K
        T_l_interp = Tl_from_Tmean(R_l=pl.R_l[iii], T_avg=T_avg_interp, a0=pl.a0[iii], **kwargs)
    

    D_l_interp = f_D_l_interp(times0) # published plot interpolated to model times, in m
    R_l_interp = pl.R_p - D_l_interp
    q_sfc_interp = therm.sph_flux(pl.R_p, a0=pl.a0[iii], T_l=T_l_interp, T_s=pl.T_s, R_l=R_l_interp, 
                                  R_p=pl.R_p, k_m=pl.k_m, **kwargs) # sfc flux in W m^-2

    ax3.plot(times0, pl.q_sfc[iii]*1e3, c='xkcd:black', label='this work')
    ax3.plot(times0, q_sfc_interp*1e3, c='xkcd:blue', label='Thiriet interp')
    df = pd.read_csv(compare_dir+'/q_sfc.csv', header=None, names=['time', 'value'],
                     index_col=False) # in Gyr, mW m^-3
    ax3.plot(df['time'], df['value'], c='xkcd:red', label='Thiriet digitised')
    ax3.legend(frameon=False, fontsize=14)
    ax3.set_xlabel('Time (Gyr)', fontsize=16)
    ax3.set_ylabel('$q_{sfc}$ (mW m$^{-2}$)', fontsize=16)

#     plt.tight_layout()

    fig0, ax0 = plt.subplots(1,1,  figsize=(4,4))
    ax0.plot(times0, pl.q_sfc[iii]*1e3 - q_sfc_interp*1e3, c='xkcd:grey')
    ax0.set_xlabel('Time (Gyr)', fontsize=14)
    ax0.set_ylabel('$\Delta q_{sfc}$ (mW m$^{-2}$)', fontsize=14)
    ax0.set_title('Mean error: $\pm$'+'%.2f'%np.mean(np.absolute(pl.q_sfc[iii]*1e3 - q_sfc_interp*1e3))+' mW m$^{-2}$', 
                  fontsize=14)
    plt.tight_layout()
    if plots_save:
        fig3.savefig(fig_path+pl.ident+'_test_qsfc.pdf')
        fig0.savefig(fig_path+pl.ident+'_q_error.pdf')

def plot_Tavg(pl, ax3=None, compare_dir=None, fig_path=None, plots_save=False, ident=None, **kwargs):
    """ sanity check on T_avg """
    if ident is None:
        ident = pl.ident
    
    t = pl.t # time points of ode solutions in s

    if ax3 is None:
        fig3, ax3 = plt.subplots(1, 1, figsize=(5,5))
        
    # plot your T_avg calculation but using interpolated published D_l    
    t_D_l, f_D_l_interp = interp_benchmark(path=compare_dir+'/D_l.csv', yscale=1e3) # in Gyr, m
    
    # select model time points in interpolation range
    temp = t*1e-9/p.years2sec # in Gyr
    iii = np.where((temp>=t_D_l.min()) & (temp<=t_D_l.max()))
    times0 = temp[iii] # time points of ODE solver subset to interpolation range
    D_l_interp = f_D_l_interp(times0) # published D_l at model time points in m
    R_l_interp = pl.R_p - D_l_interp # m
    
    T_avg_interp = therm.T_mean(T_m=pl.T_m[iii], T_l=pl.T_l[iii], R_p=pl.R_p, R_l=R_l_interp, 
                                R_c=pl.R_c, a0=pl.a0[iii], T_s=pl.T_s, k_m=pl.k_m, **kwargs)
    ax3.plot(times0, pl.T_avg[iii], c='xkcd:black', label='this work')
    ax3.plot(times0, T_avg_interp, c='xkcd:blue', label='this work with Thiriet D_l')
    df = pd.read_csv(compare_dir+'/T_avg.csv', header=None, names=['time', 'value'],
                     index_col=False) # in Gyr, mW m^-3
    ax3.plot(df['time'], df['value'], c='xkcd:red', label='Thiriet digitised')
    ax3.legend(frameon=False, fontsize=14)
    ax3.set_xlabel('Time (Gyr)', fontsize=16)
    ax3.set_ylabel('$T_{avg}$ (K)', fontsize=16)  
    ax3.set_title('blue should match red')
    
    plt.tight_layout()
    if plots_save:
        fig3.savefig(fig_path+pl.ident+'_test_Tavg.pdf')
        
def plot_one(ax, x, y, xlabel, ylabel, labelsize=15, legsize=16, ticksize=12, line_args=None,
             text=None, xticks=True, ylim=None, label=None,labelpad=None, fontname=None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(3,3))
    ax.set_xlabel(xlabel, fontsize=labelsize, fontname=fontname, labelpad=labelpad)
    ax.set_ylabel(ylabel, fontsize=labelsize, fontname=fontname, labelpad=labelpad)
    if not xticks:
        ax.set_xticks([])
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.plot(x, y, label=label, **line_args)
    ax.set_xlim(x.min(), x.max())
    ax.tick_params(axis='x', labelsize=ticksize) 
    ax.tick_params(axis='y', labelsize=ticksize) 
    if text is not None:
        ax.text(0.05, 0.95, text, ha='left', va='top', transform=ax.transAxes, fontsize=legsize)
    return ax

def plot_setxlabel(axes, label, style, labelsize=15, fontname=None, labelpad=None, **kwargs):
    try:
        ny, nx = np.shape(axes) # assume 2D
        for ii, ax in enumerate(axes):
            for jj, a in enumerate(ax):
                if (style=='all') or ((style=='bottom') and (ii==ny-1)):
                    a.set_xlabel(label, fontsize=labelsize, fontname=fontname)
                if (style=='bottom') and (ii<ny-1):
                    a.set_xticks([])
    except ValueError: # 1D
        ny = 1
        try:
            for jj, a in enumerate(axes):
                a.set_xlabel(label, fontsize=labelsize, fontname=fontname, labelpad=labelpad)
        except TypeError: # single
            axes.set_xlabel(label, fontsize=labelsize, fontname=fontname, labelpad=labelpad)
                
def Tl_from_Tmean(R_c=None, R_l=None, R_p=None, T_avg=None, T_s=None, a0=None, k_m=None, **kwargs):
    """solved for T_l using sympy"""
    return (k_m*(4*R_l**4 - 4*R_l**3*R_p - 3*R_l**2*R_p**2 + 2*R_l*R_p**3 + R_p**4)*(-60*Ea*R_c**6*k_m + 60*Ea*R_c**3*R_l**3*k_m + 30*Ea*R_c**3*R_l**2*R_p*k_m + 30*Ea*R_c**3*R_l*R_p**2*k_m - 30*Ea*R_l**5*R_p*k_m - 30*Ea*R_l**4*R_p**2*k_m + 120*R_b*R_c**3*R_l**3*T_avg*a_rh*k_m - 60*R_b*R_c**3*R_l**2*R_p*T_avg*a_rh*k_m - 60*R_b*R_c**3*R_l*R_p**2*T_avg*a_rh*k_m - 8*R_b*R_l**8*a0*a_rh + 14*R_b*R_l**7*R_p*a0*a_rh + 9*R_b*R_l**6*R_p**2*a0*a_rh - 20*R_b*R_l**5*R_p**3*a0*a_rh - 60*R_b*R_l**5*R_p*T_s*a_rh*k_m - 10*R_b*R_l**4*R_p**4*a0*a_rh - 30*R_b*R_l**4*R_p**2*T_s*a_rh*k_m + 18*R_b*R_l**3*R_p**5*a0*a_rh - 120*R_b*R_l**3*R_p**3*T_avg*a_rh*k_m + 180*R_b*R_l**3*R_p**3*T_s*a_rh*k_m + R_b*R_l**2*R_p**6*a0*a_rh + 60*R_b*R_l**2*R_p**4*T_avg*a_rh*k_m - 30*R_b*R_l**2*R_p**4*T_s*a_rh*k_m - 4*R_b*R_l*R_p**7*a0*a_rh + 60*R_b*R_l*R_p**5*T_avg*a_rh*k_m - 60*R_b*R_l*R_p**5*T_s*a_rh*k_m) - 2*sqrt(15)*sqrt(Ea*k_m**3*(60*Ea*R_c**6*k_m - 60*Ea*R_c**3*R_l**2*R_p*k_m - 60*Ea*R_c**3*R_l*R_p**2*k_m + 15*Ea*R_l**4*R_p**2*k_m + 30*Ea*R_l**3*R_p**3*k_m + 15*Ea*R_l**2*R_p**4*k_m - 240*R_b*R_c**3*R_l**3*T_avg*a_rh*k_m + 120*R_b*R_c**3*R_l**2*R_p*T_avg*a_rh*k_m + 120*R_b*R_c**3*R_l*R_p**2*T_avg*a_rh*k_m + 16*R_b*R_l**8*a0*a_rh - 28*R_b*R_l**7*R_p*a0*a_rh - 18*R_b*R_l**6*R_p**2*a0*a_rh + 40*R_b*R_l**5*R_p**3*a0*a_rh + 120*R_b*R_l**5*R_p*T_s*a_rh*k_m + 20*R_b*R_l**4*R_p**4*a0*a_rh + 60*R_b*R_l**4*R_p**2*T_s*a_rh*k_m - 36*R_b*R_l**3*R_p**5*a0*a_rh + 240*R_b*R_l**3*R_p**3*T_avg*a_rh*k_m - 360*R_b*R_l**3*R_p**3*T_s*a_rh*k_m - 2*R_b*R_l**2*R_p**6*a0*a_rh - 120*R_b*R_l**2*R_p**4*T_avg*a_rh*k_m + 60*R_b*R_l**2*R_p**4*T_s*a_rh*k_m + 8*R_b*R_l*R_p**7*a0*a_rh - 120*R_b*R_l*R_p**5*T_avg*a_rh*k_m + 120*R_b*R_l*R_p**5*T_s*a_rh*k_m))*(R_c - R_l)*(R_l - R_p)**2*(2*R_l + R_p)**2*(R_c**2 + R_c*R_l + R_l**2))/(30*R_b*R_l**2*a_rh*k_m**2*(R_l - R_p)**2*(2*R_l + R_p)**2*(4*R_l**4 - 4*R_l**3*R_p - 3*R_l**2*R_p**2 + 2*R_l*R_p**3 + R_p**4))

def powerspectrum_RMS(path, amplitude=False): # try to calcuate RMS from digitized power spectrum
    df = pd.read_csv(path, header=None, names=['degree', 'value'], index_col=False) 
    ls = np.array(df['degree'])
    S = np.array(df['value'])
    RMS_l = []
    for ii, l in enumerate(ls):
        val = np.sqrt(S[ii]/(2*l + 1))
        if amplitude:
            val = val**2
        RMS_l.append(val)
    return sum(RMS_l)

def eta_from_Ra(rho=None, g=None, alpha=None, dT=None, d=None, kappa=None, Ra=None):
    return rho*g*alpha*dT*d**3/(kappa*Ra)

def Ra_from_RaB(F=None, dT=1000, l=None, k=None, Ra_B=None):
    # convert between basal heating Ra and internal Ra, given F flux into base
    return Ra_B * k*dT/(F*l) # basal heating Ra_B

def Ra_from_RaF_2(F=74e-3, dT=1000, l=750e3, kappa=8e-7, rho=3300, c_p=1200, Ra_F=2.4e6): # moresi and parsons
    return Ra_F*(rho*c_p*kappa*dT)/(l*F)

def Ra_from_RaF(F=None, dT_m=None, k=None, l=None, Ra_F=None, **kwargs): # F is surface flux
    return Ra_F/(l*F/(k*dT_m))

def benchmark_thermal_plots(ident, show_qsfc_error=False, show_Tavg=False, names=None, pl_update_args=None,
                            model_update_args=None, **kwargs):
    if names is None:
            names = {'T_avg':('$T_{avg}$ (K)',1), 
                     'q_sfc':('$q_{sfc}$ (mW m$^{-2}$)',1e3),
                     'D_l':('$D_l$ (km)',1e-3),
                     'q_core':('$q_{B}$ (mW m$^{-2}$)', 1e3), 
                     'urey':('Ur',1),
                     'T_l':('$T_l$ (K)',1), 
                    }
    planet_kwargs = eval('ins.'+ident+'_in') 
    model_kwargs = eval('ins.'+ident+'_run')
    if pl_update_args is not None:
        planet_kwargs.update(pl_update_args)
    if model_update_args is not None:
        model_kwargs.update(model_update_args)
    pl = tp.TerrestrialPlanet(**planet_kwargs)
    pl = therm.solve(pl, **model_kwargs) # T_m, T_c, D_l

    print('T_mf', pl.T_m[-1])

    fig, axes = plot_output(pl, names, verbose=False, **kwargs)
    if show_qsfc_error:
        plot_qsfc_error(pl, ident=ident, **kwargs)
    if show_Tavg:
        plot_Tavg(pl, **kwargs)
    return fig, axes

def plot_vs_x(scplanets=None, lplanets=None, xname=None, ynames=None, planets2=None, fig=None, axes=None,
              labels=False, labelsize=15, legsize=12, alpha=1, legend=False, snap=4.5,labelpad=None,
              plots_save=False, s=30, ls='-', lw=1, cmap='rainbow', marker='o', legtitle=None, legendtop=False,
              colorbar=False, c='k', ylabel=True, ymin=None, ymax=None, set_ylim=True, set_xlim=False,
              zorder_l=None, zorder_sc=None, label_l=None, fname=None, ticksize=12, xmin=None, xmax=None, **kwargs):
        # for a list of planets, plot some parameter on the y axis vs. parameter x
    if (c is None) and (scplanets is not None):
        c = np.arange(len(scplanets))
#         colour = cm.get_cmap(cmap)
#         norm = colors.Normalize(vmin=0, vmax=len(planets))
    nax = len(ynames)
    if xmin is not None:
        set_xlim=True
    if axes is None:
        fig, axes = plt.subplots(1, nax, figsize=(5*nax, 4))
    
    xparam = list(xname.keys())[0]
    xlabels = list(xname.values())[0]
    yparam = list(ynames.keys())
    ylabels = list(ynames.values()) # tuple (ylabel, yscale)
    
    
    ii=0
    while ii < nax:
        try:
            ax = axes[ii]
        except TypeError: # single ax
            ax = axes
        if scplanets is not None:
            x = []
            y = []
            for ip, pl in enumerate(scplanets): # planets to plot as scatter
                data_x = eval('pl.'+xparam)*xlabels[1]
                if isinstance(data_x, Iterable):
                    data_x = data_x[-1] # get final value
                x.append(data_x)
                data_y = eval('pl.'+yparam[ii])*ylabels[ii][1]
                if isinstance(data_y, Iterable):
                    data_y = data_y[-1] # get final value
                y.append(data_y) 
                if labels:
                    ax.annotate(xy=(data_x,data_y), s=pl.ident[0:2], fontsize=legsize)
            sc = ax.scatter(x, y, c=c, s=s, marker=marker, cmap=cmap, zorder=zorder_sc)
            if colorbar and (ii==0):
                plt.colorbar(sc)
        if lplanets is not None:
            x = []
            y = []
            try:
                for ip, pl in enumerate(lplanets): # planets to plot as line
                    t = pl.t
                    it = min(enumerate(t), key=lambda x: abs(snap - x[1]*p.sec2Gyr))[0] # get time index nearest to desired snap given in Gyr
                    data_x = eval('pl.'+xparam)*xlabels[1]
                    if isinstance(data_x, Iterable):
                        data_x = data_x[it] # if an evolution model then take certain snap
                    x.append(data_x)
                    data_y = eval('pl.'+yparam[ii])*ylabels[ii][1]
                    if isinstance(data_y, Iterable):
                        data_y = data_y[it] # if an evolution model then take certain snap
                    y.append(data_y) 
            except TypeError: # if given a single planet (not iterable) - get values across evol
                x = eval('lplanets.'+xparam)*xlabels[1]
                y = eval('lplanets.'+yparam[ii])*ylabels[ii][1]
            # sort
            if isinstance(y, Iterable) and isinstance(x, Iterable):
                x, y = zip(*sorted(zip(x, y)))
            elif isinstance(y, Iterable) and (not isinstance(x, Iterable)):
                x, y = zip(*sorted(zip([x]*np.ones_like(y), y)))
            elif (not isinstance(y, Iterable)) and isinstance(x, Iterable):
                x, y = zip(*sorted(zip(x, [y]*np.ones_like(x))))    

            ax.plot(x, y, ls=ls, c=c, lw=lw, alpha=alpha, zorder=zorder_l, label=label_l)
        if ylabel:
            ax.set_ylabel(ylabels[ii][0], fontsize=labelsize, labelpad=labelpad)
            ax.tick_params(axis='y', labelsize=ticksize) 
        else:
            ax.yaxis.set_ticklabels([])
        ax.tick_params(axis='x', labelsize=ticksize) 
        
        # log scale for viscosity
        if (yparam[ii] is 'eta_m') or (yparam[ii] is 'nu_m') or (yparam is 'Ra_i'):
            ax.set_yscale('log')
        if (xparam is 'eta_m') or (xparam is 'nu_m') or (xparam is 'eta_0') or (xparam is 'Ra_i'):
            ax.set_xscale('log')
        
        if set_ylim:
            ax.set_ylim(ymin, ymax)
        if set_xlim:
            ax.set_xlim(xmin, xmax)

        if legend:

            if legendtop:
                legend=ax.legend(frameon=False, fontsize=legsize,  
                                 borderaxespad=0, title=legtitle,  #mode="expand", 
                                 loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=2,)
           
            else:
                legend=ax.legend(frameon=False, fontsize=legsize, loc='upper left', 
                             bbox_to_anchor= (1.05, 0.9),
                             borderaxespad=0, ncol=1, title=legtitle)
            if legtitle is not None:
                    plt.setp(legend.get_title(),fontsize=legsize)
                    legend._legend_box.align = "left"
        ii+=1
    plot_setxlabel(axes, xlabels[0], 'all', labelsize=labelsize, labelpad=labelpad, **kwargs)
    
    plt.tight_layout()
    if plots_save:
        if fname is None:
            fname = 'scatter_'+xparam
        plt.savefig(fig_path+fname+'.pdf', bbox_inches='tight') 
    return fig, axes