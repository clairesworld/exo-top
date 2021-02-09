import sys
sys.path.insert(0, '/home/cmg76/Works/exo-top/')
import numpy as np
from exotop import aspect_postprocessing2 as asp

cases_to_ascii = ['Ra1e8-eta1e7-wide']
for case in cases_to_ascii:
    dat = asp.Aspect_Data(directory='/raid1/cmg76/aspect/model-output/output-'+case+'/', verbose=False)
    dat.write_ascii(A=None, fname='Tf_'+case, ext='.txt', path='/raid1/cmg76/aspect/model-input/ascii/', n=None, default_field='T')