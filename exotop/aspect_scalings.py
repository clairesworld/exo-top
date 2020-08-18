import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import random as rand
import pickle as pkl
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
           ax=None,savefig=True,settitle=True, setxlabel=True, legend=True, labelsize=16):
    if sigma==2:
        qs = [2.5, 50, 97.5]
    elif sigma==1:
        qs = [16, 50, 84]
    time, v_rms, nsteps = read_evol(case, i=11, path=path)
    # what is the probability distribution of i from t1 to end?
    i_time = np.nonzero(time > t1)
#     print('i_time', i_time)
    rms_list = []
    peak_list = []
    t_used = []
    for ii in i_time[0]:
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
    if plot:
        if ax is None:
            fig = plt.figure()
            ax = plt.gca()
        ax.hist(rms_list, color='r', histtype='step', label='rms')
        ax.hist(peak_list, color='b', histtype='step', label='peak')
        ax.axvline(x=np.median(rms_list), color='k', ls='--', label='median')
        ax.axvline(x=np.median(peak_list), color='k', ls='--')
        ax.axvline(x=np.mean(rms_list), color='k', ls='-', lw=0.5, label='mean')
        ax.axvline(x=np.mean(peak_list), color='k', ls='-', lw=0.5)
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

def plot_evol(case, i, fig=None, ax=None, savefig=True, fig_path='figs/', fend='_f.png',
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
    ax.plot(time, y*yscale, c=c, label=label)
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

def case_subplots(cases, fig_path='figs/', labels=None, labelsize=16, labelpad=5, t1=None, save=True,
                  fig_path='/raid1/cmg76/aspect/figs/'):
    # rows are cases, columns are v_rms, q, hist
    ncases = len(cases)
    fig, axes = plt.subplots(ncases, 3, figsize=(15, ncases*2.5))
    for ii, case in enumerate(cases):
        print('reading',case)
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
        rect = patches.Rectangle((0,0), t1[ii], ax.get_ylim()[1] - ax.get_ylim()[0],
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
        rect = patches.Rectangle((0,0), t1[ii], ax.get_ylim()[1] - ax.get_ylim()[0],
                                 edgecolor='None',facecolor='k', alpha=0.2, zorder=0)
        ax.add_patch(rect)
        
        ax = axes[ii, 2]
        _, _, fig, ax = pd_top(case, t1[ii], path='model-output/', sigma=2, plot=True, fig=fig, ax=ax, 
                              savefig=False, settitle=False, setxlabel=setxlabel, legend=legend,
                              labelsize=labelsize)
        
    fig.tight_layout()
    if save:
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path+'cases-Ra.png', bbox_inches='tight')
    return fig, ax

def plot_h_Ra(Ra, eta, t1=None, path='model-output/', fig_path='/raid1/cmg76/aspect/figs/', picklefrom=None, pickleto=None,
              save=True, fname='h_Ra.png', plotpd=True, sigma=2, labelsize=16, 
              c_peak='xkcd:forest green', c_rms='xkcd:periwinkle'):
    # Ra is list of strings, t1 is a list of numbers the same length
    if t1 is None:
        t1 = [0.005]*len(Ra)
    if picklefrom is not None:
        h_peak, h_rms = pkl.load(open( fig_path+'data/'+picklefrom, "rb" ))
    else:
        h_peak = np.zeros((len(Ra), 3))
        h_rms = np.zeros((len(Ra), 3))
        for ii, r in enumerate(Ra):
            case = 'Ra'+r+'-eta'+eta+'-wide'
            print(case)
            # assume time-dependent convection (if steady state then shouldn't matter)
            # assuming quasi steady state from 0.005 but defo need to check
            h_peak[ii,:], h_rms[ii,:], _, _ = pd_top(case, t1=t1[ii], path=path, plot=plotpd, sigma=sigma) 
            # also plot probability distribution and velocity evolution to check steady state 
    if pickleto is not None:
        pkl.dump((h_peak, h_rms), open( fig_path+'data/'+pickleto, "wb" ))
           
    fig = plt.figure()
    ax = plt.gca()
#     ax.plot([float(s) for s in Ra], h_peak[:,1], label='peak')
#     ax.plot([float(s) for s in Ra], h_rms[:,1], label='RMS')
    ax.errorbar([float(s) for s in Ra], h_peak[:,1], yerr=[h_peak[:,1]-h_peak[:,0], h_peak[:,2]-h_peak[:,1]], 
                label='peak', fmt='-o', c=c_peak, alpha=0.9)
    ax.errorbar([float(s) for s in Ra], h_rms[:,1], yerr=[h_rms[:,1]-h_rms[:,0], h_rms[:,2]-h_rms[:,1]], 
                label='rms', c=c_rms, alpha=0.9)
    ax.legend(frameon=False)
    ax.set_xscale('log')
    ax.set_ylabel('dynamic topography', labelsize=labelsize)
    ax.set_xlabel('Ra', labelsize=labelsize)
    ax.set_title('$\Delta \eta = $'+eta, labelsize=labelsize)
    if save:
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path+fname, bbox_inches='tight')
    return fig, ax

def plot_h_eta(Ra, eta, t1=None, path='model-output/', fig_path='/raid1/cmg76/aspect/figs/', picklefrom=None, pickleto=None,
               save=True, fname='h_eta.png', plotpd=True, sigma=2, labelsize=16,
               c_peak='xkcd:forest green', c_rms='xkcd:periwinkle'):
    # Ra is list of strings, t1 is a list of numbers the same length
    if t1 is None:
        t1 = [0.005]*len(Ra)
    if picklefrom is not None:
        h_peak, h_rms = pkl.load(open( fig_path+'data/'+picklefrom, "rb" ))
    else:
        h_peak = np.zeros((len(Ra), 3))
        h_rms = np.zeros((len(Ra), 3))
        for jj, e in enumerate(eta):
            print(Ra)
            case = 'Ra'+Ra+'-eta'+e+'-wide'
            print(case)
            # assume time-dependent convection (if steady state then shouldn't matter)
            # assuming quasi steady state from 0.005 but defo need to check
            h_peak[jj, :], h_rms[jj, :], _, _ = pd_top(case, t1=t1[jj], path=path, plot=plotpd, sigma=sigma) 
            # also plot probability distribution and velocity evolution to check steady state 
    if pickleto is not None:
        pkl.dump((h_peak, h_rms), open( fig_path+'data/'+pickleto, "wb" ))
           
    fig = plt.figure()
    ax = plt.gca()
    ax.errorbar([float(s) for s in Ra], h_peak[:,1], yerr=[h_peak[:,1]-h_peak[:,0], h_peak[:,2]-h_peak[:,1]], 
                label='peak', fmt='-o', c=c_peak, alpha=0.9)
    ax.errorbar([float(s) for s in Ra], h_rms[:,1], yerr=[h_rms[:,1]-h_rms[:,0], h_rms[:,2]-h_rms[:,1]], 
                label='rms', c=c_rms, alpha=0.9)
    ax.legend(frameon=False)
    ax.set_xscale('log')
    ax.set_ylabel('dynamic topography', labelsize=labelsize)
    ax.set_xlabel('$\Delta \eta = $', labelsize=labelsize)
    ax.set_title('Ra ='+Ra, labelsize=labelsize)
    if save:
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig.savefig(fig_path+fname, bbox_inches='tight')
    return fig, ax
