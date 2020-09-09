import sys, os
sys.path.insert(0,'/home/cmg76/Works/exo-top/')
import numpy as np
from exotop import aspect_scalings as sc


Ra_ls = ['3e5', '1e6', '3e6', '1e7', '3e7', '1e8', '3e8']
eta_ls = ['1e5', '1e6', '1e7', '1e8']

load = True
dump = False

t1 = np.array([[0.5, 0.5, 0.3, 0.25, 0.3, 0.1, 0.04],
               [0, 0.4, 0.3, 0.3, 0.17, 0.15, 0.005], 
               [0, 0, 0.45, 0.35, 0.175, 0.15, 0.002], 
               [0, 0, 0, 0.4, 0.15, 0.002, 0.002]]) # eta outside, Ra inside
# t1=None

#for ii, eta in enumerate(eta_ls):

 #   fig, ax = sc.case_subplots(
  #      ['Ra3e5-eta'+eta+'-wide','Ra1e6-eta'+eta+'-wide','Ra3e6-eta'+eta+'-wide','Ra1e7-eta'+eta+'-wide','Ra3e7-eta'+eta+'-wide','Ra1e8-eta'+eta+'-wide','Ra3e8-eta'+eta+'-wide'], 
   #     labels=['Ra=3e5','Ra=1e6', 'Ra=3e6', 'Ra=1e7', 'Ra=3e7', 'Ra=1e8', 'Ra=3e8'],
    #    t1=t1[ii], loadpickle=load, dumppickle=dump, fname='all-eta'+eta+'.png', suptitle='$\Delta \eta$='+eta,
     #   includepd=True, # turn on once you know where steady state starts
      #  )

    #print(eta, 'Ra subplots complete')

    
includeRa = np.array([['1e6', '3e6', '1e7', '3e7'], 
                      ['1e6', '3e6', '1e7', '3e7', '1e8'], 
                      ['3e6', '1e7', '3e7', '1e8'], 
                      ['1e7', '3e7']]) # eta outside, Ra inside
    
for ii, eta in enumerate(eta_ls): 
    try:
        fig, ax = sc.plot_h_vs(Ra=includeRa[ii], eta=eta, t1=t1[ii], sigma=2, title='$\Delta \eta$='+eta, xlabel='Ra',
                           loadpickle=True, dumppickle=dump,  plotpd=False, fname='h-Ra_eta'+eta+'.png', )
        print(eta, 'Ra scaling complete')
    except Exception as e:
        print('ERROR', e)

for ii, eta in enumerate(eta_ls):
    try:
        fig, ax = sc.plot_h_vs(Ra=includeRa[ii], eta=eta, t1=t1[ii], sigma=2, title='$\Delta \eta$='+eta, 
                               xlabel='$\Delta T_{rh} \delta_u d_m$',
                               loadpickle=True, dumppickle=dump,  plotpd=False, fname='h-x_eta'+eta+'.png', x_components=True,
                               fit=False)
        print(eta, 'dT*delta*d scaling complete')
    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fename = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fename, exc_tb.tb_lineno, e)




