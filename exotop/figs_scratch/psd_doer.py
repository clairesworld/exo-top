import fractals as fract
import numpy as np

data_path = '/home/claire/Works/aspect/runs/model-output/'
case = 'Ra3e8-eta1e8-wide-ascii'
fract.dct_spectrum_avg(case, ts0=137000, tsf=137900, x_res=1, t_res=100, norm='ortho',
                       test=True, data_path=data_path,
                       plot=True, load=False, dump=False, save=False, L_x=8, l_min=30, l_max=300,
                       dim=True, d=600, dT=442, alpha=4e-5  # Lees table 1-2 Ra=1e6 ; Venus: d=2700, dT=3000, alpha=2e-5
                       )