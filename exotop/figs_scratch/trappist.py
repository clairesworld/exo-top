import numpy as np
from model_1D import the_results as pltt
from model_1D import evolve as evol
from model_1D import parameters as p
from model_1D import oceans, topography
import sh_things as sh
from model_1D import inputs


def single(default, dist_res=100, yvars=['dyn_top_rms', 'max_ocean'], yscales=[1, p.TO], t_eval=None, age=6,
           maxi=None, mini=None, names=None, run_kwargs=None, update_kwargs=None, n_sigma=1, it=-1, **kwargs):
    if maxi is None:
        maxi = [300e3, 2.5e12, 2000, 200e3]
    if mini is None:
        mini = [240e3, 1.5e10, 1000, 50e3]
    if names is None:
        names = ['Ea', 'eta_pre', 'T_m0', 'D_l0']

    if default is not None:
        pl_kwargs = eval('inputs.' + default + '_in').copy()
        model_kwargs = eval('inputs.' + default + '_run').copy()
    else:
        pl_kwargs= {}  # use defaults given in terrestrialplanet.py
        model_kwargs = {}  # initial conditions defaults given in thermal.py
    if update_kwargs is not None:
        pl_kwargs.update(update_kwargs)
    if run_kwargs is not None:
        model_kwargs.update(run_kwargs)
    pl_ensemble = evol.bulk_planets_mc(n=dist_res, names=names, mini=mini, maxi=maxi, pl_kwargs=pl_kwargs,
                                       model_kwargs=model_kwargs, t_eval=t_eval, log=False, **kwargs)

    if it is None:
        it = pltt.age_index(pl_ensemble[0].t, age, age_scale=p.sec2Gyr)  # should all be evaluated at same time

    print('Ra', np.array([np.log10(vars(pl)['Ra_i_eff'][it]) for pl in pl_ensemble]))

    for yvar, yscale in zip(yvars, yscales):
        col = np.array([vars(pl)[yvar][it] for pl in pl_ensemble]) * yscale
        y_av = np.mean(col)
        y_std = np.std(col)
        y_upper = y_av + y_std * n_sigma  # todo for log scape
        y_lower = y_av - y_std * n_sigma
        print(yvar, y_av, '+', y_upper, '-', y_lower)
    return None


# degree, phi0 = sh.load_model_spectrum_pkl(fname='base_spectrum_l1.pkl', path='/home/claire/Works/exo-top/exotop/figs_scratch/')
# postprocessors = ['topography', 'ocean_capacity']
# ids = ['TRAPPIST1e', 'TRAPPIST1f', 'TRAPPIST1h', 'TRAPPIST1h']
# for plid in ids:
#     print('\n', plid)
#     single(plid, dist_res=100, yvars=['dyn_top_rms', 'max_ocean'], yscales=[1, p.TO**-1],
#            run_kwargs=None, update_kwargs=None, postprocessors=postprocessors, phi0=phi0)


# planet e
print('e')
R_p = 0.788
M_p = 0.692
g = 8.01
x_h2o = [0.3*1e-2]  # CMF=0.25
for x in x_h2o:
    oceans.min_topo(x, R_p*p.R_E, M_p, rms_1=1000, tol=0.01, spectrum_fname='spectrum_-2.pkl', verbose=False)
    print('h_max=', topography.strength_max(g))

# planet f
print('f')
R_p = 1.045
M_p = 1.039
g = 9.32
x_h2o = [1.9*1e-2]  # CMF=0.25
for x in x_h2o:
    oceans.min_topo(x, R_p*p.R_E, M_p, rms_1=1000, tol=0.01, spectrum_fname='spectrum_-2.pkl', verbose=False)
    print('h_max=', topography.strength_max(g))

# planet g
print('g')
R_p = 1.129
M_p = 1.321
g = 10.15
x_h2o = [3.5*1e-2]  # CMF=0.25
for x in x_h2o:
    oceans.min_topo(x, R_p*p.R_E, M_p, rms_1=1000, tol=0.01, spectrum_fname='spectrum_-2.pkl', verbose=False)
    print('h_max=', topography.strength_max(g))
