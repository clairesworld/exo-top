import numpy as np
import scalings as sc

fig, ax = sc.case_subplots(
			   ['Ra3e5-eta1e5-wide','Ra1e6-eta1e5-wide','Ra3e6-eta1e5-wide','Ra1e7-eta1e5-wide','Ra3e7-eta1e5-wide','Ra1e8-eta1e5-wide','Ra3e8-eta1e5-wide','Ra1e9-eta1e5-wide'], 
			   labels=['Ra=3e5','Ra=1e6', 'Ra=3e6', 'Ra=1e7', 'Ra=3e7', 'Ra=1e8', 'Ra=3e8', 'Ra=1e9'],
                           t1=[0.6, 0.6, 0.35, 0.25, 0.25, 0.15, 0.035, 0.017])


fig, ax = sc.case_subplots(
		          ['Ra1e8-eta1e5-wide', 'Ra1e8-eta1e6-wide', 'Ra1e8-eta1e7-wide'], 
			  labels=['$\Delta \eta$=1e5', '$\Delta \eta$=1e6', '$\Delta \eta$=1e7'],
                          t1=[0.1, 0.1, 0.12])
