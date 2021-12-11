import numpy as np
import model_1D.parameters as p

class InteriorStructure():
    # adapted from noack 2021 valid to 2 M_E

    def __init__(self, verbose=False, **kwargs):
        if verbose:
            print('getting interior structure')
        # define default attributes
        # default_attr = {
        #     # constants
        #     'M_p': 1 * p.M_E,
        #     'X_Fe': 0.35815,  # weight iron fraction of planet from 0.15 to 0.8
        #     'num_Fe_m': 0.063632,  # iron number in mantle (Mg/Fe ratio in olivine)
        # }
        #
        # # add input parameters, use default if not given
        # default_attr.update(kwargs)
        # self.__dict__.update((k, v) for k, v in default_attr.items())


    def init_structure(self, M_p=p.M_E, X_Fe=0.35815, num_Fe_m=0.063632, which='hot', CMF=None, T_c0=None,
                       m_Fe=55.845e-3, m_Mg=24.305e-3, m_Si=28.0855e-3, m_O=15.999e-3, **kwargs):
        self.M_p = M_p
        self.X_Fe = X_Fe
        self.num_Fe_m = num_Fe_m
        self.R_p = (7030 - 1840 * self.X_Fe) * 1e3 * (self.M_p/p.M_E)**0.282  # planet radius in km
        self.X_Fe_m = 2 * self.num_Fe_m * m_Fe / ((2 * (1 - self.num_Fe_m) * m_Mg + self.num_Fe_m * m_Fe) + m_Si + 4 * m_O)  # iron mass fraction in mantle

        if CMF is None:
            self.CMF = (self.X_Fe - self.X_Fe_m) / (1 - self.X_Fe_m)  # core mass fraction
        else:
            self.CMF = CMF

        if which == 'hot' or which == 'warm':
            self.R_c = 4850e3 * self.CMF**0.328 * (self.M_p / p.M_E)**0.266
        elif which == 'cold':
            self.R_c = 4790e3 * self.CMF ** 0.328 * (self.M_p / p.M_E) ** 0.266
        else:
            raise Exception('specify hot, warm or cold start')

        self.rho_c_av = self.CMF * self.M_p / (4*np.pi/3 * self.R_c**3)  # core average density
        self.rho_m_av = (1 - self.CMF) * self.M_p / (4*np.pi/3 * (self.R_p**3 - self.R_c**3))  # mantle avg density

        self.g_0 = p.G * self.M_p / (self.R_p**2)  # surface gravity
        self.g_c = p.G * self.CMF * self.M_p / (self.R_c**2)
        self.g_m_av = (self.g_0 + self.g_c) / 2  # mantle avg gravity, off by a few % from radial avg profile

        self.p_cmb = self.g_m_av * self.rho_m_av * (self.R_p - self.R_c)  # CMB pressure in Pa - avg error 3%

        # initial temperatures
        if T_c0 is None:
            if which == 'hot':
                self.T_cmb = 5400 * (self.p_cmb / 140e9) ** 0.48 / (1 - np.log(1 - self.num_Fe_m))
        else:
            self.T_cmb = T_c0
            print('fixing input T_cmb0 =', T_c0)




# test = InteriorStructure()
# test.init_structure(CMF=None, which='hot')
# print(test.__dict__)