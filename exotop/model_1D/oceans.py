from . import parameters
import numpy as np

def ocean_vol_from_hpeak(h_peak, EO=1, R_p=parameters.R_E):
        # R_p = 6051.88e3):
        # EO = 1.35e9 * 1000**3
        vol = 4 / 3 * np.pi * ((R_p + h_peak) ** 3 - R_p ** 3)
        return vol / EO
