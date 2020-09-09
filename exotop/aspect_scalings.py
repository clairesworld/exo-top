import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import random as rand
import pickle as pkl
import collections
import six
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
    
def pd_top(case, t1, path='model-output/', fig_path='/raid1/cmg76/aspect/figs/', sigma=2, plot=True, fig=None, 
           ax=None, savefig=True, settitle=True, setxlabel=True, legend=True, labelsize=16, pickleto=None, picklefrom=None,
           c_peak='xkcd:forest green', c_rms='xkcd:periwinkle'):
    if sigma==2:
        qs = [2.5, 50, 97.5]
    elif sigma==1:
        qs = [16, 50, 84]    


    if (picklefrom is not None) and (os.path.exists(fig_path+'data/'+picklefrom)):
        try:
            peak_list, rms_list = pkl.load(open( fig_path+'data/'+picklefrom, "rb" ))
        except ValueError:
            peak_list, rms_list = pkl.load(open( fig_path+'data/'+picklefrom, "rb" ), protocol=2)
    
    else:    
        time, v_rms, nsteps = read_evol(case, i=11, path=path)
        # what is the probability distribution of i from t1 to end?
        i_time = np.nonzero(time > t1)[0]
    #     print('i_time', i_time)
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
    #             print('D.T. file', ii, 'not found')
                pass
    
    if pickleto is not None:
        pkl.dump((peak_list, rms_list), open( fig_path+'data/'+pickleto, "wb" ))
    
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

#     print('timesteps using:', np.array(t_used).min(), 'to', np.array(t_used).max())
    return np.percentile(peak_list, qs), np.percentile(rms_list, qs), fig, ax


def pd_h_components(case, t1, data_path='model-output/', fig_path='/raid1/cmg76/aspect/figs/', sigma=2,
                    pickleto=None, picklefrom=None, plot=False, savefig=False, settitle=True, setxlabel=True, legend=True, 
                    labelsize=16,  fig=None, ax=None):
    # probability distribution of h' = f(x) for single case
    # it's a bit annoying because you're switching methods...
    if sigma==2:
        qs = [2.5, 50, 97.5]
    elif sigma==1:
        qs = [16, 50, 84]    
   

    if (picklefrom is not None) and (os.path.exists(fig_path+'data/'+picklefrom)):
        try:
            x_list = pkl.load(open( fig_path+'data/'+picklefrom, "rb" ))
        except ValueError:
            x_list = pkl.load(open( fig_path+'data/'+picklefrom, "rb" ), protocol=2)
    
    else:    
        dat = post.Aspect_Data(directory=data_path+'output-'+case+'/')
        dat.read_statistics(verbose=False)
        time = dat.stats_time
        snaps = dat.read_stats_sol_files()
        i_time = np.nonzero(time > t1)[0]
        n_quasi = np.unique(snaps[i_time]) # find graphical snapshots within time range

        x_list = []
        for n in n_quasi:
            n = int(n)
            x, y, z, u, v, _ = dat.read_velocity(n, verbose=False)
            x, y, z, T = dat.read_temperature(n, verbose=False)
            x = dat.h_components(n, T=T, u=u, v=v)
            x_list.append(x)
        
    if pickleto is not None:
        pkl.dump(x_list, open( fig_path+'data/'+pickleto, "wb" ))
    
    if plot:
        if ax is None:
            fig = plt.figure()
            ax = plt.gca()
        ax.hist(x_list, color='k', histtype='step')
        ax.axvline(x=np.median(x_list), color='k', ls='--', label='median')
        ax.axvline(x=np.mean(x_list), color='k', ls='-', lw=1, label='mean')
        ax.yaxis.set_ticks([])
        if legend:
            ax.legend(frameon=False)
        if setxlabel:
            ax.set_xlabel('$\Delta T_{rh} \delta_u$', fontsize=labelsize)
        if settitle:
            ax.set_title(case, fontsize=labelsize)
        if savefig:
            if not os.path.exists(fig_path):
                os.makedirs(fig_path)    
            fig.savefig(fig_path+case+'_x_hist.png')
    
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

def case_subplots(cases, labels=None, labelsize=16, labelpad=5, t1=None, save=True, dt_xlim=(0.0,0.065), fname='cases.png', data_path='model-output/',
                  fig_path='/raid1/cmg76/aspect/figs/', loadpickle=False, dumppickle=False, suptitle='', includepd=True):
    # rows are cases, columns are v_rms, q, hist
    ncases = len(cases)
    fig, axes = plt.subplots(ncases, 3, figsize=(15, ncases*2.5))
    if t1 is None:
        t1 = [0]*ncases
    for ii, case in enumerate(cases):
        if (os.path.exists(data_path+'output-'+case)):
            print('loading',case)
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

            ax = axes[ii, 0]
            fig, ax = plot_evol(case, i=11, fig=fig, ax=ax, savefig=False, ylabel='rms velocity', 
                                c='k', settitle=False, setxlabel=setxlabel, setylabel=setylabel, 
                                labelsize=labelsize, labelpad=labelpad, legend=False)

            # Create a Rectangle patch to mark "transient" times
            rect = patches.Rectangle((ax.get_xlim()[0],ax.get_ylim()[0]), t1[ii], ax.get_ylim()[1] - ax.get_ylim()[0],
                                     edgecolor='None',facecolor='k', alpha=0.2, zorder=0)
            ax.add_patch(rect)

            try:
                ax.text(0.01, 0.95, labels[ii],horizontalalignment='left', verticalalignment='top',
                        transform = ax.transAxes, fontsize=labelsize)
            except:
                pass

            ax = axes[ii, 1]
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

            if includepd:
                ax = axes[ii, 2]
                picklefile = case+'_pdtop.pkl'
                if loadpickle:
                    picklefrom = picklefile
                else:
                    picklefrom = None
                if dumppickle:
                    pickleto = picklefile
                else:
                    pickleto = None
                _, _, fig, ax = pd_top(case, t1[ii], path='model-output/', sigma=2, plot=True, fig=fig, ax=ax, 
                                      savefig=False, settitle=False, setxlabel=setxlabel, legend=legend,
                                      labelsize=labelsize, pickleto=pickleto, picklefrom=picklefrom)
                ax.set_xlim(dt_xlim[0], dt_xlim[1]) # for fair comparison
        
        else:
            print(case, 'not found')
    plt.suptitle(suptitle, fontsize=labelsize, y=1.02)
    fig.tight_layout()
    if save:
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path+fname, bbox_inches='tight')
    return fig, ax

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
    
    

def plot_h_vs(Ra, eta, t1=None, path='model-output/', fig_path='/raid1/cmg76/aspect/figs/', loadpickle=False, dumppickle=False,
              save=True, fname='h.png', plotpd=True, sigma=2, labelsize=16, xlabel='', title='',
              c_peak='xkcd:forest green', c_rms='xkcd:periwinkle', x_components=False, fit=False):
    # Ra or eta is list of strings, t1 is a list of numbers the same length
    cases, x_var = get_cases_list(Ra, eta)  
    if t1 is None:
        t1 = [0]*len(x_var)

    h_peak = np.zeros((len(x_var), 3))
    h_rms = np.zeros((len(x_var), 3))
    x = np.zeros((len(x_var), 3)) 
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

        try:
            h_peak[ii,:], h_rms[ii,:], _, _ = pd_top(case, t1=t1[ii], path=path, plot=plotpd, sigma=sigma,
                                                 pickleto=pickleto, picklefrom=picklefrom) 
        except Exception as e:
            print('aspect_scalings.py line 332:', e)
            h_peak[ii,:], h_rms[ii,:] = ([np.nan,np.nan,np.nan],[np.nan,np.nan,np.nan])
        if x_components:
            # instead of plotting vs Ra or eta, plot vs theoretical components of scaling relationship
            picklefile = case+'_pdx.pkl'
            oldload = loadpickle
            olddump = dumppickle
            loadpickle = True
            dumppickle = True
            if loadpickle:
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
                                            picklefrom=picklefrom, fig_path=fig_path)
            except Exception as e:
                print('aspect_scalings.py line 350:', e)
                x[ii,:] = (np.nan,np.nan,np.nan)
            loadpickle = oldload
            dumppickle = olddump

        else:
            x[ii,:] = float(x_var[ii])

    
    yerr_peak = [h_peak[:,1]-h_peak[:,0], h_peak[:,2]-h_peak[:,1]]
    yerr_rms = [h_rms[:,1]-h_rms[:,0], h_rms[:,2]-h_rms[:,1]]
    if x_components:
        print('x', x)
        try:
            xerr = [x[:,1]-x[:,0], x[:,2]-x[:,1]]
        except Exception as e1:
            print('aspect_scalings.py line 364', e1)

    else:
        xerr = None
    
    fig = plt.figure()
    ax = plt.gca()
#     ax.plot([float(s) for s in Ra], h_peak[:,1], label='peak')
#     ax.plot([float(s) for s in Ra], h_rms[:,1], label='RMS')
    ax.errorbar(x[:,1], h_peak[:,1], yerr=yerr_peak, xerr=xerr,
                label='peak', fmt='-o', c=c_peak, alpha=0.9, capsize=5)
    ax.errorbar(x[:,1], h_rms[:,1], yerr=yerr_rms, xerr=xerr,
                label='rms', fmt='-o', c=c_rms, alpha=0.9, capsize=5)
    ax.legend(frameon=False)
    ax.set_xscale('log')
    ax.set_ylabel('dynamic topography', fontsize=labelsize)
    ax.set_xlabel(xlabel, fontsize=labelsize)
    ax.set_title(title, fontsize=labelsize)
    
    if fit:
        print('fitting x =',x[:,1],'h =',  h_rms[:,1])
        try:
            m, b = fit_h(x[:,1], h_rms[:,1], h_err=None)
            ax.plot(x[:,1], m*x[:,1] + b, c='xkcd:grey', alpha=0.5, ls='--')
            print('C =', b)
        except Exception as e:
            print('aspect_scalings.py line 386:', e)
    
    
    if save:
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path+fname, bbox_inches='tight')
    return fig, ax


def fit_h(x, h, h_err):
    idx = np.argwhere(np.isnan(x)==False)
    p = np.polyfit(x[idx], h[idx], 1)
    return p[1], p[0]
    

# def plot_h_eta(Ra, eta, t1=None, path='model-output/', fig_path='/raid1/cmg76/aspect/figs/', loadpickle=False, dumppickle=False,
#                save=True, fname='h_eta.png', plotpd=True, sigma=2, labelsize=16,
#                c_peak='xkcd:forest green', c_rms='xkcd:periwinkle'):
#     # Ra is list of strings, t1 is a list of numbers the same length
#     if t1 is None:
#         t1 = [0.005]*len(Ra)
#     else:
#         h_peak = np.zeros((len(Ra), 3))
#         h_rms = np.zeros((len(Ra), 3))
#         for jj, e in enumerate(eta):
#             case = 'Ra'+Ra+'-eta'+e+'-wide'
#             picklefile = case+'_pdtop.pkl'
#             if loadpickle:
#                 picklefrom = picklefile
#             else:
#                 picklefrom = None
#             if dumppickle:
#                 pickleto = picklefile
#             else:
#                 pickleto = None
#             # assume time-dependent convection (if steady state then shouldn't matter)
#             h_peak[jj, :], h_rms[jj, :], _, _ = pd_top(case, t1=t1[jj], path=path, plot=plotpd, sigma=sigma,
#                                                        pickleto=pickleto, picklefrom=picklefrom) 

#     fig = plt.figure()
#     ax = plt.gca()
#     ax.errorbar([float(s) for s in eta], h_peak[:,1], yerr=[h_peak[:,1]-h_peak[:,0], h_peak[:,2]-h_peak[:,1]], 
#                 label='peak', fmt='-o', c=c_peak, alpha=0.9, capsize=5)
#     ax.errorbar([float(s) for s in eta], h_rms[:,1], yerr=[h_rms[:,1]-h_rms[:,0], h_rms[:,2]-h_rms[:,1]], 
#                 label='rms', fmt='-o', c=c_rms, alpha=0.9, capsize=5)
#     ax.legend(frameon=False)
#     ax.set_xscale('log')
#     ax.set_ylabel('dynamic topography', fontsize=labelsize)
#     ax.set_xlabel('$\Delta \eta = $', fontsize=labelsize)
#     ax.set_title('Ra ='+Ra, fontsize=labelsize)
#     if save:
#         if not os.path.exists(fig_path):
#             os.makedirs(fig_path)
#         fig.savefig(fig_path+fname, bbox_inches='tight')
#     return fig, ax

