import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import random as rand
import pickle as pkl
import collections
import six
from scipy.optimize import curve_fit
from exotop import aspect_postprocessing2 as post # from exotop
cwd = os.getcwd()

def read_topo_stats(case, snap, path='model-output/'):
    df = pd.read_csv(path+'output-'+case+'/dynamic_topography_surface.'+'{:05}'.format(snap), header=None,
                     names=['x', 'y', 'h'], 
                     skiprows=1,
                     index_col=False, delimiter=r"\s+", engine='python')
    return df['x'], df['h']

def trapznorm(A):
    mean = np.trapz(A)/(len(A)-1)
#     print('original mean:',mean)
    return A - mean

def trapzmean(A):
    return np.trapz(A)/(len(A)-1)

def peak_and_rms(h):
    return np.max(h), np.sqrt(trapzmean(h**2))

def read_evol(case, i, path='model-output/', skiprows=None):
    # e.g. return time, column i, nsteps (from statistics file)
    df = pd.read_csv(path+'output-'+case+'/statistics', header=None,
                 skiprows=skiprows,
                 index_col=False, delimiter=r"\s+", engine='python', comment='#')
#     print(df)
    return np.array(df.iloc[:, 1]), np.array(df.iloc[:, i-1]), len(df.index)

def top_profile(case, savefig=True, fig_path='/raid1/cmg76/aspect/figs/', path='model-output/'):
    time, y, nsteps = read_evol(case, i=2, path=path)
    snap = nsteps - 2 # honestly not sure why it's -2 but seems to work haha
    x, h = read_topo_stats(case, snap)
    # normalize to 0 mean
    h_norm = trapznorm(h)
    fig = plt.figure()
    plt.plot(x, h_norm)
    plt.xlabel('x')
    plt.ylabel('dynamic topography')
    plt.title(case)
    print('mean:',trapzmean(h_norm))
    print('max:', np.max(h_norm))
    print('min:', np.min(h_norm))
    if savefig:
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)    
        fig.savefig(fig_path+case+'_h_'+'{:05}'.format(snap)+'.png')
    
def pd_quasiss(case, i, t1, path='model-output/'):
    time, y, nsteps = read_evol(case, i, path=path)
    x, h = read_topo_stats(case, nsteps-2)
    h_norm = trapznorm(h)
    peak, rms = peak_and_rms(h_norm)
    # what is the probability distribution of i from t1 to end?
    i_time = np.nonzero(time > t1)
    print('pretending transience ends at timestep',i_time[0][0])
    fig = plt.figure()
    plt.gca().hist(y[i_time])

def get_T_params(case, t1=0, path='model-output/', pickleto=None, picklefrom=None, plotTz=False,
                 setylabel=True, setxlabel=True, savefig=True,
                 fig_path='/raid1/cmg76/aspect/figs/', fig=None, ax=None, 
                 legend=True, labelsize=16, **kwargs):
    if (os.path.exists(path+'output-'+case)):
        flag=False
        n_quasi = []
        if (picklefrom is not None) and (os.path.exists(fig_path+'data/'+picklefrom)):
            try:
                dT_rh, dT_m, delta, d_m, D_l, T_l = pkl.load(open( fig_path+'data/'+picklefrom, "rb" ))
                print('loading T profiles from', case)
            except ValueError:
                dT_rh, dT_m, delta, d_m, D_l, T_l = pkl.load(open( fig_path+'data/'+picklefrom, "rb" ), 
                                                             protocol=2)
                print('loading T profiles from', case)
            if not dT_rh: # if empty
                flag=True
        else:
            flag=True

        if flag:  
            print('calculating T profiles from', case)
            dT_rh = []
            dT_m = []
            delta = []
            d_m = []
            D_l = []
            T_l = []
            dat = post.Aspect_Data(directory=path+'output-'+case+'/')
            dat.read_statistics(verbose=False)
            time = dat.stats_time
            snaps = dat.read_stats_sol_files()
            i_time = np.nonzero(time > t1)[0]
            n_quasi = np.unique(snaps[i_time]) # find graphical snapshots within time range

            for n in n_quasi:
                n = int(n)
                x, y, z, u, v, _ = dat.read_velocity(n, verbose=False)
                x, y, z, T = dat.read_temperature(n, verbose=False)
                dT_rh_n, dT_m_n, delta_n, d_m_n, D_l_n, T_l_n = dat.h_components(n, T=T, u=u, v=v, cut=True)
                dT_rh.append(dT_rh_n)
                dT_m.append(dT_m_n)
                delta.append(delta_n)
                d_m.append(d_m_n)
                D_l.append(D_l_n)
                T_l.append(T_l_n) 

        if pickleto is not None:
            pkl.dump((dT_rh, dT_m, delta, d_m, D_l, T_l), open( fig_path+'data/'+pickleto, "wb" ))  
            
        if plotTz and ((not flag) or (len(n_quasi)>2)):
            if not flag: # load last time step if haven't yet
                dat = post.Aspect_Data(directory=path+'output-'+case+'/')
                n = dat.final_step()
                x, y, z, u, v, _ = dat.read_velocity(n, verbose=False)
                x, y, z, T = dat.read_temperature(n, verbose=False)
                dT_rh_n = dT_rh[-1]
                delta_n = delta[-1]
                D_l_n = D_l[-1]
                T_l_n = T_l[-1]
            fig, ax = dat.plot_profile(T, fig=fig, ax=ax, xlabel='', ylabel='', c='k', lw=1)
            ax.axhline(D_l_n, label='$z_{lid}$', c='xkcd:red orange', lw=0.5)
            ax.axhline(D_l_n-delta_n, label=r'$\delta$', c='xkcd:tangerine', lw=0.5)
            ax.text(0, D_l_n-delta_n, r'$\delta = $'+'{:04.2f}'.format(delta_n), ha='left', va='top',
                   color='xkcd:tangerine')
            ax.plot([T_l_n, T_l_n], [0, D_l_n], ls='--', alpha=0.5, lw=0.5, c='xkcd:red orange')
            ax.plot([T_l_n+dT_rh_n, T_l_n+dT_rh_n], [0, D_l_n-delta_n], ls='--', alpha=0.5, lw=0.5, 
                    c='xkcd:tangerine')
            if legend:
                ax.legend(frameon=False)
            if setxlabel:
                ax.set_xlabel('temperature', fontsize=labelsize)
            if setylabel:
                ax.set_ylabel('depth', fontsize=labelsize)
            if savefig:
                fig.savefig(fig_path+case+'-T_z.png', bbox_inches='tight')  

        return (dT_rh, dT_m, delta, d_m, D_l, T_l), fig, ax   
    else:
        print('case', case, 'not found')
        return ([np.nan], [np.nan], [np.nan], [np.nan], [np.nan], [np.nan]), fig, ax
    
def get_h(case, t1=0, path='model-output/', pickleto=None, picklefrom=None, 
          fig_path='/raid1/cmg76/aspect/figs/', **kwargs):
    flag=False
    if (picklefrom is not None) and (os.path.exists(fig_path+'data/'+picklefrom)):
        try:
            peak_list, rms_list = pkl.load(open( fig_path+'data/'+picklefrom, "rb" ))
            print('loaded h for case', case)
        except ValueError:
            peak_list, rms_list = pkl.load(open( fig_path+'data/'+picklefrom, "rb" ), protocol=2)
            print('loaded h for case', case)
        if (not peak_list) or (not rms_list): # if stored stuff is empty
            flag=True
            pickleto = picklefrom
    else:    
        flag=True
        pickleto = picklefrom

    if flag: # load
        print('building distribution of h for case', case)
        time, v_rms, nsteps = read_evol(case, i=11, path=path)
        # what is the probability distribution of i from t1 to end?
        i_time = np.nonzero(time > t1)[0]
        rms_list = []
        peak_list = []
        t_used = []
        for ii in i_time:
            try:
                x, h = read_topo_stats(case, ii)
                h_norm = trapznorm(h)
                peak, rms = peak_and_rms(h_norm)
                rms_list.append(rms)
                peak_list.append(peak)
                t_used.append(ii)
            except FileNotFoundError:
                pass

    if pickleto is not None:
        pkl.dump((peak_list, rms_list), open( fig_path+'data/'+pickleto, "wb" ))
    return peak_list, rms_list
    
def pd_top(case, t1=0, path='model-output/', fig_path='/raid1/cmg76/aspect/figs/', sigma=2, 
           plot=True, fig=None, ax=None, savefig=True, settitle=True, setxlabel=True, legend=True, 
           labelsize=16, pickleto=None, picklefrom=None,
           c_peak='xkcd:forest green', c_rms='xkcd:periwinkle', peak_list=None, rms_list=None):
    if sigma==2:
        qs = [2.5, 50, 97.5]
    elif sigma==1:
        qs = [16, 50, 84]    

    if (peak_list is None) or (rms_list is None):
        peak_list, rms_list = get_h(case, t1, path, pickleto, picklefrom)
    
    if plot:
        if ax is None:
            fig = plt.figure()
            ax = plt.gca()
        ax.hist(rms_list, color=c_rms, histtype='step', label='rms')
        ax.hist(peak_list, color=c_peak, histtype='step', label='peak')
        ax.axvline(x=np.median(rms_list), color='k', ls='--', label='median')
        ax.axvline(x=np.median(peak_list), color='k', ls='--')
        ax.axvline(x=np.mean(rms_list), color='k', ls='-', lw=1, label='mean')
        ax.axvline(x=np.mean(peak_list), color='k', ls='-', lw=1)
        ax.yaxis.set_ticks([])
        ax.text(0.05, 0.05, 'n = {:d}'.format(len(rms_list)), ha='left', va='bottom', transform = ax.transAxes)
        if legend:
            ax.legend(frameon=False)
        if setxlabel:
            ax.set_xlabel('dynamic topography', fontsize=labelsize)
        if settitle:
            ax.set_title(case, fontsize=labelsize)
        if savefig:
            if not os.path.exists(fig_path):
                os.makedirs(fig_path)    
            fig.savefig(fig_path+case+'_h_hist.png')

    if (not peak_list) or (not rms_list): # empty
        print(case, '- h list is empty')
        return np.array([np.nan, np.nan, np.nan]), np.array([np.nan, np.nan, np.nan]), fig, ax
    return np.percentile(peak_list, qs), np.percentile(rms_list, qs), fig, ax

def pd_h_components(case, t1=0, data_path='model-output/', fig_path='/raid1/cmg76/aspect/figs/', sigma=2,
                    pickleto=None, picklefrom=None, plotTz=False, savefig=False, 
                    settitle=True, setxlabel=True, c='xkcd:pale purple',
                    legend=True, plotpd=False, labelsize=16, fig=None, ax=None):
    # probability distribution of h' = f(x) for single case
    # it's a bit annoying because you're switching methods...
    if sigma==2:
        qs = [2.5, 50, 97.5]
    elif sigma==1:
        qs = [16, 50, 84]    
   
    (dT_rh, dT_m, delta, d_m, D_l, T_l), fig, ax = get_T_params(case, t1=t1, path=path, setxlabel=setxlabel,
                                                                pickleto=pickleto, picklefrom=picklefrom, 
                                                                savefig=savefig,fig=fig, ax=ax, 
                                                                plotTz=plotTz, fig_path=fig_path)
    
    x_list = (np.array(dT_rh)/np.array(dT_m))*(np.array(delta)/np.array(d_m)) 
    if plotpd and (not plotTz):
        if ax is None:
            fig = plt.figure()
            ax = plt.gca()
        ax.hist(x_list, color=c, histtype='step')
        ax.axvline(x=np.median(x_list), color='k', ls='--', label='median')
        ax.axvline(x=np.mean(x_list), color='k', ls='-', lw=1, label='mean')
        ax.yaxis.set_ticks([])
        ax.text(0.95, 0.95, 'n = {:d}'.format(len(x_list)), ha='right', va='top', transform = ax.transAxes)
        if legend:
            ax.legend(frameon=False)
        if setxlabel:
            ax.set_xlabel(r'$\Delta T_{rh}/\Delta T_m \delta/d_m$', fontsize=labelsize)
        if settitle:
            ax.set_title(case, fontsize=labelsize)
        if savefig:
            if not os.path.exists(fig_path):
                os.makedirs(fig_path)    
            fig.savefig(fig_path+case+'_x_hist.png')
            
    if (not x_list): # empty
        print(case, '- x list is empty')
        return np.array([np.nan, np.nan, np.nan]), fig, ax
    return np.percentile(x_list, qs), fig, ax  

def plot_evol(case, i, fig=None, ax=None, savefig=True, fend='_f.png',
              ylabel='rms velocity', xlabel='time', yscale=1, c='k', settitle=True, setxlabel=True,
              setylabel=True, legend=False, labelsize=16, labelpad=5, label=None, fig_path='/raid1/cmg76/aspect/figs/'):
    if not setxlabel:
        xlabel=''
    if not setylabel:
        ylabel=''
    if ax is None:
        fig = plt.figure()
        ax = plt.gca()
    time, y, nsteps = read_evol(case, i=i)
    ax.plot(time, y*yscale, c=c, lw=0.5, label=label)
    ax.set_xlim(0,ax.get_xlim()[1])
    ax.set_xlabel(xlabel, fontsize=labelsize, labelpad=labelpad)
    ax.set_ylabel(ylabel, fontsize=labelsize, labelpad=labelpad)
    if settitle:
        ax.set_title(case, fontsize=labelsize)
    if legend:
        ax.legend(frameon=False)
    if savefig:
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path+case+fend)
    return fig, ax

def case_subplots(cases, labels=None, labelsize=16, labelpad=5, t1=None, save=True, dt_xlim=(0.0,0.065), 
                  fname='cases.png', data_path='model-output/', fig_path='/raid1/cmg76/aspect/figs/', 
                  loadpickle=False, dumppickle=False, dumppicklex=False,
                  suptitle='', includepd=True, includeTz=True, loadpicklex=False):
    # rows are cases, columns are v_rms, q, hist
    ncases = len(cases)
    ncols = 2
    if includepd:
        ncols = ncols+1
    if includeTz:
        ncols = ncols+1
    fig, axes = plt.subplots(ncases, ncols, figsize=(17, ncases*2.5))
    if t1 is None:
        t1 = [0]*ncases
    for ii, case in enumerate(cases):
        icol = 0
        if (os.path.exists(data_path+'output-'+case)):
            print('plotting summary for',case)
            if ii==ncases-1:
                setxlabel=True
            else:
                setxlabel=False
            if ii==int(np.median(range(ncases))):
                setylabel=True
            else:
                setylabel=False
            legend=True
            if ii>0:
                legend=False

            ax = axes[ii, icol]
            fig, ax = plot_evol(case, i=11, fig=fig, ax=ax, savefig=False, ylabel='rms velocity', 
                                c='k', settitle=False, setxlabel=setxlabel, setylabel=setylabel, 
                                labelsize=labelsize, labelpad=labelpad, legend=False)

            # Create a Rectangle patch to mark "transient" times
            rect = patches.Rectangle((ax.get_xlim()[0],ax.get_ylim()[0]), t1[ii], 
                                     ax.get_ylim()[1] - ax.get_ylim()[0],
                                     edgecolor='None',facecolor='k', alpha=0.2, zorder=0)
            ax.add_patch(rect)

            try:
                ax.text(0.01, 0.95, labels[ii],horizontalalignment='left', verticalalignment='top',
                        transform = ax.transAxes, fontsize=labelsize)
            except:
                pass

            icol = icol+1
            ax = axes[ii, icol]
            fig, ax = plot_evol(case, i=20, fig=fig, ax=ax, savefig=False, ylabel='heat flux', 
                                c='xkcd:light red', settitle=False, setxlabel=setxlabel, setylabel=setylabel,
                                labelsize=labelsize, labelpad=labelpad, label='top')
            fig, ax = plot_evol(case, i=19, fig=fig, ax=ax, savefig=False, ylabel='heat flux', yscale=-1,
                                c='xkcd:purple blue', settitle=False, setxlabel=setxlabel, setylabel=setylabel,
                                labelsize=labelsize, labelpad=labelpad, label='bottom', legend=legend)

            # Create a Rectangle patch to mark "transient" times
            rect = patches.Rectangle((ax.get_xlim()[0],ax.get_ylim()[0]), t1[ii], ax.get_ylim()[1] - ax.get_ylim()[0],
                                     edgecolor='None',facecolor='k', alpha=0.2, zorder=0)
            ax.add_patch(rect)

            if includeTz: # final timestep only
                icol = icol+1
                ax = axes[ii, icol]
                picklefile = case+'_pdx.pkl'
                if loadpicklex:
                    picklefrom = picklefile
                else:
                    picklefrom = None
                if dumppicklex:
                    pickleto = picklefile
                else:
                    pickleto = None
                _, fig, ax = get_T_params(case, t1=t1[ii], path=data_path, savefig=False,
                                          pickleto=pickleto, picklefrom=picklefrom, plotTz=True,
                                          setxlabel=False, setylabel=False,
                                          legend=False, 
                                          fig_path=fig_path, fig=fig, ax=ax)   
                if legend:
                    ax.legend(frameon=False)
                if setxlabel:
                    ax.set_xlabel('temperature', fontsize=labelsize)
                if setylabel:
                    ax.set_ylabel('depth', fontsize=labelsize)
            
            if includepd:
                icol = icol+1
                ax = axes[ii, icol]
                picklefile = case+'_pdtop.pkl'
                if loadpickle:
                    picklefrom = picklefile
                else:
                    picklefrom = None
                if dumppickle:
                    pickleto = picklefile
                else:
                    pickleto = None
                _, _, fig, ax = pd_top(case, t1[ii], path=data_path, sigma=2, plot=True, 
                                       fig=fig, ax=ax, savefig=False, settitle=False, setxlabel=setxlabel, 
                                       legend=legend, labelsize=labelsize, 
                                       pickleto=pickleto, picklefrom=picklefrom)
                ax.set_xlim(dt_xlim[0], dt_xlim[1]) # for fair comparison
                    
        else:
            print(case, 'not found')
    plt.suptitle(suptitle, fontsize=labelsize, y=1.02)
    fig.tight_layout()
    if save:
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path+fname, bbox_inches='tight')
    return fig, axes

def get_cases_list(Ra, eta):
    # either Ra or eta is iterable
    if (isinstance(Ra, collections.Iterable) and not isinstance(Ra, six.string_types)):
        x_var = Ra
        cases = ['Ra'+r+'-eta'+eta+'-wide' for r in Ra]
    elif (isinstance(eta, collections.Iterable) and not isinstance(eta, six.string_types)):
        x_var = eta
        cases = ['Ra'+Ra+'-eta'+e+'-wide' for e in eta]
    else:
        raise Exception('Ra or eta must be iterable')
    return cases, x_var  

def plot_h_vs(Ra=None, eta=None, t1=None, path='model-output/', fig_path='/raid1/cmg76/aspect/figs/', 
              loadpickle=False, dumppickle=False, loadpicklex=False, plotTz=False,
              save=True, fname='h.png', plotpd=False, sigma=2, labelsize=16, xlabel='', title='',
              c_peak='xkcd:forest green', c_rms='xkcd:periwinkle', x_components=False, legend=True,
              fit=False, fitRa=None, fitfn='line', cases=None, x_var=None, logx=True, logy=True,
              fig=None, ax=None):
    # Ra or eta is list of strings, t1 is a list of numbers the same length
    if cases is None:
        cases, x_var = get_cases_list(Ra, eta)  
    if t1 is None:
        t1 = [0]*len(x_var)
    if fitRa is None:
        fitRa = Ra

    h_peak = np.zeros((len(x_var), 3))
    h_rms = np.zeros((len(x_var), 3))
    x = np.zeros((len(x_var), 3)) 
    peak_all = []
    rms_all = []
                 
    for ii, case in enumerate(cases):
        picklefile = case+'_pdtop.pkl'
        if loadpickle:
            picklefrom = picklefile
        else:
            picklefrom = None
        if dumppickle:
            pickleto = picklefile
        else:
            pickleto = None
        # assume time-dependent convection (if steady state then shouldn't matter)
        peak_list, rms_list = get_h(case, t1=t1[ii], path=path, pickleto=pickleto, picklefrom=picklefrom)
        peak_all.append((peak_list, float(x_var[ii])))
        rms_all.append((rms_list, float(x_var[ii])))
                 
        try:
            h_peak[ii,:], h_rms[ii,:], _, _ = pd_top(case, plot=False, 
                                                     peak_list=peak_list, rms_list=rms_list) 
        except Exception as e:
            print('aspect_scalings.py:', e, '\n setting h all nan for case', case)
            h_peak[ii,:], h_rms[ii,:] = ([np.nan,np.nan,np.nan],[np.nan,np.nan,np.nan])        
                 
        if x_components:
            # instead of plotting vs Ra or eta, plot vs theoretical components of scaling relationship
            picklefile = case+'_pdx.pkl'
            if loadpicklex:
                picklefrom = picklefile
            else:
                picklefrom = None
            if dumppickle:
                pickleto = picklefile
            else:
                pickleto = None
            # assume time-dependent convection (if steady state then shouldn't matter)
            try:
                x[ii,:], _, _ = pd_h_components(case, t1=t1[ii], data_path=path, sigma=sigma, pickleto=pickleto,
                                                picklefrom=picklefrom, fig_path=fig_path, plotTz=plotTz)
            except Exception as e:
                print('aspect_scalings.py line 350:', e, '\n setting x all nan for case', case)
                x[ii,:] = (np.nan,np.nan,np.nan)
        else:
            x[ii,:] = float(x_var[ii])
   
    yerr_peak = [h_peak[:,1]-h_peak[:,0], h_peak[:,2]-h_peak[:,1]]
    yerr_rms = [h_rms[:,1]-h_rms[:,0], h_rms[:,2]-h_rms[:,1]]
    if x_components:
        try:
            xerr = [x[:,1]-x[:,0], x[:,2]-x[:,1]]
        except Exception as e1:
            print('aspect_scalings.py line 364', e1)
    else:
        xerr = None
    
    if ax is None:
        fig = plt.figure()
        ax = plt.gca()
    ax.errorbar(x[:,1], h_peak[:,1], yerr=yerr_peak, xerr=xerr,
                label='peak', fmt='-^', c=c_peak, alpha=0.9, capsize=5)
    ax.errorbar(x[:,1], h_rms[:,1], yerr=yerr_rms, xerr=xerr,
                label='rms', fmt='-o', c=c_rms, alpha=0.9, capsize=5)
    
    if fit:
#         fitx = x_var # todo: only fit subset?
#         fitidx = np.where(np.intersect1d(x[:,1], fitx))[0]
        if (len(x_var)>1):
            fitx = [[a[1]]*len(a[0]) for a in rms_all]
            fith = [a[0] for a in rms_all]
            flatfitx = [item for sublist in fitx for item in sublist]
            flatfith = [item for sublist in fith for item in sublist]
            expon, const = fit_h(flatfitx, flatfith)
            xprime = [a[1] for a in rms_all]
            hprime = const*xprime**expon
            ax.plot(xprime, hprime, c=c_rms, ls='--', lw=0.5, 
                    label='{:.2e} Ra^{:.3f}'.format(const, expon), zorder=100)
        else:
            print('not enough points to fit -- Ra', Ra, 'eta', eta)
    
    if legend:
        ax.legend(frameon=False)
    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')
    ax.set_ylabel('dynamic topography', fontsize=labelsize)
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_title(title, fontsize=labelsize)
    
    if save:
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path+fname, bbox_inches='tight')
    return fig, ax

from scipy import stats
def fit_h(x, h, plot=True):
    x1 = np.log10(np.array(x))
    h1 = np.log10(np.array(h))
    slope, intercept, r_value, p_value, std_err = stats.linregress(x1,h1)
    return slope, 10**intercept

def fit_h_sigma(x, h, h_err=None, fn='line'):
    def line(x, a, b):
        return a * x + b
#     def expon(x, C, n):
#         return C * x**n
    
    idx = np.nonzero(np.isnan(x)==False)[0]
    x_fit = np.log10(x[idx])
    h_fit = np.log10(h[idx])
    if h_err is not None:
        h_err = np.log10(h_err[idx])
    print('fitting x =',x_fit,'h =',  h_fit)
    if fn=='line':
        popt, pcov = curve_fit(line, x_fit, h_fit, sigma=h_err)
        print('slope:', popt[0], 'intercept:', popt[1])

    return 10**(popt[1] + popt[0]*x) # h evaluated at x

def subplots_h_vs(Ra_ls, eta_ls, regime_grid, c_regimes, loadpickle=True, dumppickle=False, save=True, 
                  sigma=2, t1=None, fit=False, loadpicklex=False, 
                  path='model-output/', fig_path='/raid1/cmg76/aspect/figs/', fname='h_Ra_all.png'):
    fig, axes = plt.subplots(2,2, figsize=(10,10))
    flaxes = axes.flatten()
 
    for ii, eta in enumerate(eta_ls): 
        ax = flaxes[ii]
        Ra_steady = [Ra_ls[j] for j in np.nonzero(regime_grid[ii]=='steady')[0]] 
        Ra_trans = [Ra_ls[j] for j in np.nonzero(regime_grid[ii]=='transitional')[0]] 
        Ra_chaos = [Ra_ls[j] for j in np.nonzero(regime_grid[ii]=='chaotic')[0]] 

        # steady
        fig, ax = plot_h_vs(Ra=Ra_steady, eta=eta, t1=t1, sigma=sigma, fig=fig, ax=ax, 
                               c_rms=c_regimes[0], c_peak=c_regimes[0],
                               loadpickle=loadpickle, dumppickle=dumppickle, plotpd=False, 
                               save=False, fit=True)
        # trans
        fig, ax = plot_h_vs(Ra=Ra_trans, eta=eta, t1=t1, sigma=sigma, fig=fig, ax=ax, 
                               c_rms=c_regimes[1], c_peak=c_regimes[1],
                               loadpickle=loadpickle, dumppickle=dumppickle, plotpd=False, 
                               save=False, fit=True)       
        # chaotic
        fig, ax = plot_h_vs(Ra=Ra_chaos, eta=eta, t1=t1, sigma=sigma, fig=fig, ax=ax, 
                               c_rms=c_regimes[2], c_peak=c_regimes[2],
                               loadpickle=loadpickle, dumppickle=dumppickle, plotpd=False, 
                               save=False, fit=True,
                               title='$\Delta \eta$='+eta, xlabel='Ra')
        
        # show regime boundaries
        try:
            ax.axvline(float(Ra_steady[-1])*2, c='k', lw=0.5)
            ax.text(float(Ra_steady[-1]), ax.get_ylim()[0], 'steady state', 
                    fontsize=8, va='bottom', ha='center')
        except:
            print('no steady regime for eta', eta)
        try:
            ax.axvline(float(Ra_trans[-1])*2, c='k', lw=0.5)
            ax.text(float(Ra_trans[-1]), ax.get_ylim()[0], 'transitional', 
                    fontsize=8, va='bottom', ha='center')
        except:
            print('no transitional regime for eta', eta)
        try:
            ax.text(float(Ra_chaos[-1]), ax.get_ylim()[0], 'chaotic', 
                    fontsize=8, va='bottom', ha='center')
        except:
            print('no chaotic regime for eta', eta)
    if save:
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path+fname, bbox_inches='tight')  
    return fig, axes

def plot_convection_regimes(Ra, eta, regime_grid, path='model-output/', fig_path='/raid1/cmg76/aspect/figs/', loadpickle=False, 
                            dumppickle=False, save=True, fname='regimes.png', labelsize=16, sigma=2,
                            overploth=False, nlevels=10, clist=None, cmap_contours='spring', **kwargs):
    # Ra and eta are lists of strings
    if clist is None:
        cmap = plt.cm.get_cmap('jet', 3)
    else:
        cmap = cmap_from_list(clist, cmap_name='regimes')
    
    plot_grid = np.zeros((len(eta), len(Ra)))
    for x, xval in enumerate(Ra):
        for y, yval in enumerate(eta):
            desc = regime_grid[y,x]
            if desc=='chaotic':
                plot_grid[y,x] = 3
            elif desc=='transitional':
                plot_grid[y,x] = 2
            elif desc=='steady':
                plot_grid[y,x] = 1
            else: # i.e. no convection
                plot_grid[y,x] = np.nan # label grid with numbers
    m = np.ma.masked_where(np.isnan(plot_grid),plot_grid)
    im = plt.imshow(m, origin='bottom', aspect='equal', interpolation='None', 
                    vmin=1, vmax=3)
    ax = plt.gca()
    
    for x, xval in enumerate(Ra):
        ax.axvline(x=x-0.5, c='k', lw=0.5)
    for y, yval in enumerate(eta):
        ax.axhline(y=y-0.5, c='k', lw=0.5)
           
    ax.set_xlabel('Ra', fontsize=labelsize)
    ax.set_ylabel(r'$\Delta \eta$', fontsize=labelsize)
    ax.set_xticks(np.arange(len(Ra)))
    ax.set_yticks(np.arange(len(eta)))
    ax.set_xticklabels(Ra)
    ax.set_yticklabels(eta)
    
    cbar = plt.colorbar(im, ticks=[1.5, 2, 2.5])
    cbar.ax.set_yticklabels(['steady', 'transitional', 'chaotic'])
    
    if overploth: # do h rms contours, only if already stored rms
        h_grid = np.zeros((len(eta), len(Ra)))
        h_grid[:] = np.nan
        for iir in range(1,4): # separate contours for each regime (don't link)
            for x, Raval in enumerate(Ra):
                for y, etaval in enumerate(eta):
                    if plot_grid[y,x] == iir: # if not np.isnan( plot_grid[y,x] ):
                        case = 'Ra'+Raval+'-eta'+etaval+'-wide'
                        peak, rms, _, _ = pd_top(case, path=path, fig_path=fig_path, 
                                                 sigma=sigma, plot=False, picklefrom=case+'_pdtop.pkl')
                        h_grid[y,x] = rms[1]
            CS = ax.contour(h_grid, nlevels, cmap=cmap_contours) 
            ax.clabel(CS, inline=1, fontsize=10)
              
    if save:
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        plt.savefig(fig_path+fname, bbox_inches='tight')

from matplotlib.colors import LinearSegmentedColormap
def cmap_from_list(clist, n_bin=None, cmap_name=''):
    if n_bin is None:
        n_bin = len(clist)
    cm = LinearSegmentedColormap.from_list(cmap_name, clist, N=n_bin)
    return cm
