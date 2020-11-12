import sys
sys.path.insert(0, '/home/cmg76/Works/exo-top/')
import numpy as np
from exotop import aspect_postprocessing2 as asp
from exotop import aspect_scalings as sc
from setup_postprocessing import Ra_ls, eta_ls, t1, end

# plot just T(z) at final timestep, and extract all temperature params

for ii, eta in enumerate(eta_ls):  # eta_ls
    for jj, Ra in enumerate(Ra_ls):
        case = 'Ra' + Ra + '-eta' + eta + '-wide'
        sc.get_T_params(case, t1=t1[ii, jj], pickleto=case + '_pdx.pkl',
                        picklefrom=case + '_pdx.pkl', plotTz=True,
                        xlabel='temperature', ylabel='depth', savefig=True)

# case = 'Ra3e6-eta1e7-wide'
# sc.get_T_params(case, t1=0.7,
#                 picklefrom=case+'_pdx.pkl', plotTz=True,
#                 xlabel='temperature', ylabel='depth', savefig=True)
#
# case = 'Ra3e8-eta1e5-wide'
# sc.get_T_params(case, t1=0.045,
#                 picklefrom=case+'_pdx.pkl', plotTz=True,
#                 xlabel='temperature', ylabel='depth', savefig=True)
