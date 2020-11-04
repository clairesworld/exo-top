import sys
sys.path.insert(0,'/home/cmg76/Works/exo-top/')
import numpy as np
from exotop import aspect_scalings as sc
from exotop import aspect_postprocessing2 as asp
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt

Ra_ls = ['1e6', '3e6', '1e7', '3e7', '1e8', '3e8'] # 3e5 no convection
eta_ls = ['1e5', '1e6', '1e7', '1e8']
t1 = np.array([[0.5, 0.3, 0.25, 0.7, 0.2, 1], # eta 1e5
               [0.9, 0.3, 0.3, 0.4, 0.3, 1],  # eta 1e6
               [0, 0.7, 0.35, 0.4, 0.3, 1],  # eta 1e7
               [0, 0, 0.6, 0.55, 1, 1]]) # eta 1e8

# plot just T(z) at final timestep, and extract all temperature params

# for ii, eta in enumerate(eta_ls): # eta_ls
#     for jj, Ra in enumerate(Ra_ls):
#         case = 'Ra'+Ra+'-eta'+eta+'-wide'
#         sc.get_T_params(case, t1=t1[ii,jj], pickleto=case+'_pdx.pkl', 
#                         picklefrom=case+'_pdx.pkl', plotTz=True,
#                         xlabel='temperature', ylabel='depth', savefig=True)


# case = 'Ra3e6-eta1e7-wide'
# sc.get_T_params(case, t1=0.7, 
#                 picklefrom=case+'_pdx.pkl', plotTz=True,
#                 xlabel='temperature', ylabel='depth', savefig=True)

# case = 'Ra3e8-eta1e5-wide'
# sc.get_T_params(case, t1=0.045, 
#                 picklefrom=case+'_pdx.pkl', plotTz=True,
#                 xlabel='temperature', ylabel='depth', savefig=True)

# plot summaries

# for ii, eta in enumerate(eta_ls): # eta_ls
#     cases_ii = ['Ra'+Ra+'-eta'+eta+'-wide' for Ra in Ra_ls]
#     labels_ii = ['Ra='+Ra for Ra in Ra_ls]
#     fig, ax = sc.case_subplots(
#         cases_ii,
#         labels=labels_ii,
#         t1=t1[ii], save=True, loadpickle=True, dumppickle=False,
#         includeTz=True,
#         loadpicklex=True, dumppicklex=True,
#         fname='all-eta'+eta+'.png', suptitle='$\Delta \eta$='+eta,
#         includepd=True, # turn on once you know where steady state starts
#         includegraphic=True,
#        )

# compare 64 and 129 resolution for Ra=3e7
# fig, ax = sc.case_subplots(
#     ['Ra3e7-eta1e5-wide', 'Ra3e7-eta1e5-wide-128'],
#         ['64', '128'], 
#     t1=[0.25, 0.25], loadpickle=True, dumppickle=dump, fname='all-Ra3e7-res.png', suptitle='$\Delta \eta$=1e5, Ra=3e7',
#     includepd=True, # turn on once you know where steady state starts
#    )

eta_ls = ['1e5', '1e6', '1e7', '1e8']
fig, ax = sc.case_subplots(
    ['Ra3e8-eta'+eta+'-wide' for eta in eta_ls],
    labels=[r'$\Delta \eta$='+eta for eta in eta_ls],
    t1=None, save=True, loadpickle=False, dumppickle=True,
    includeTz=False, loadpicklex=False, dumppicklex=True,
    fname='all-Ra3e8.png', suptitle='Ra = 3e8',
    includepd=False,  # turn on once you know where steady state starts
       )

print('summary plots complete')
#
# cases_to_ascii = ['Ra1e8-eta1e7-wide']
# for case in cases_to_ascii:
#     dat = asp.Aspect_Data(directory='/raid1/cmg76/aspect/model-output/output-'+case+'/', verbose=False)
#     dat.write_ascii(A=None, fname='Tf_'+case, ext='.txt', path='/raid1/cmg76/aspect/model-input/ascii/', n=None, default_field='T')