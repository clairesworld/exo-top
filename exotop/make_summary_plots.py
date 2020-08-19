import sys
sys.path.insert(0,'/home/cmg76/Works/exo-top/')
import numpy as np
from exotop import aspect_scalings as sc

# set times at which quasi steady-state begins
t1_Ra = [0.6, 0.6, 0.35, 0.25, 0.25, 0.15, 0.035, 0.017]
t1_eta = [0.15, 0.15, 0.15]

fig, ax = sc.case_subplots(
			   ['Ra3e5-eta1e5-wide','Ra1e6-eta1e5-wide','Ra3e6-eta1e5-wide','Ra1e7-eta1e5-wide','Ra3e7-eta1e5-wide','Ra1e8-eta1e5-wide','Ra3e8-eta1e5-wide','Ra1e9-eta1e5-wide'], 
			   labels=['Ra=3e5','Ra=1e6', 'Ra=3e6', 'Ra=1e7', 'Ra=3e7', 'Ra=1e8', 'Ra=3e8', 'Ra=1e9'],
                           t1=t1_Ra, loadpickle=True, dumppickle=True, fname='cases-Ra.png', suptitle='$\Delta \eta$=1e5')
print('Ra subplots complete')

fig, ax = sc.plot_h_Ra(Ra=['3e5', '1e6', '3e6', '1e7', '3e7', '1e8'], eta='1e5', t1=t1_Ra, sigma=2, 
                       loadpickle=True, dumppickle=True,  plotpd=False)
print('Ra scaling complete')

fig, ax = sc.case_subplots(
		          ['Ra1e8-eta1e5-wide', 'Ra1e8-eta1e6-wide', 'Ra1e8-eta1e7-wide'], 
			  labels=['$\Delta \eta$=1e5', '$\Delta \eta$=1e6', '$\Delta \eta$=1e7'],
                          loadpickle=True, dumppickle=True, t1=t1_eta, fname='cases-eta.png', suptitle='Ra = $1\times 10^8$')
print('eta subplots complete')

fig, ax = sc.plot_h_eta(Ra='1e8', eta=['1e5', '1e6', '1e7'], t1=t1_eta, sigma=2, 
                        loadpickle=True, dumppickle=True, plotpd=False)
print('eta scaling complete')
