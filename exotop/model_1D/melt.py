import pyMelt as m
import numpy as np

def decompression_melting(P_GPa, Tp, rho_m=None, alpha_m=None, cp_m=None, X_H2O_wtpt=0.1, P_lithosphere=None, p_max_melt=10):

    Tp_C = Tp - 273.15

    # init lithologies
    lz = m.lithologies.matthews.klb1()
    # px = m.lithologies.matthews.kg1()
    hz = m.lithologies.shorttle.harzburgite()

    hlz = m.hydrousLithology(lz, H2O=X_H2O_wtpt, continuous=True,
                             phi=0.5)  # continuous melting - only matters for hydrous melting because water is extracted
    # hlz_batch = m.hydrousLithology(lz, 0.1)

    # put into pymelt mantle object - mostly lherzolite mantle
    mantle = m.mantle([hlz, hz], [8, 2], ['HLz', 'Hz'])

    # start = time.time()
    column = mantle.adiabaticMelt(Tp_C)
    # end = time.time()
    # print('decompression melt calculation in', (end - start), 'seconds')
    # f, a = column.plot()

    oib = m.geosettings.intraPlate(column, P_lithosphere=P_lithosphere, relative_density=0.2)
    melt_flux_vol = oib.melt_flux  # m3 s-1
    melt_flux_mass = melt_flux_vol * rho_melt
    mdot = # fractional rate of melt generation
    V = # I think this is volume of crust

    melt_cooling = rho_melt * mdot * (cp * (T - Ts) * L) * V

