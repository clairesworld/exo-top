""" adapted from Nimmo+ 2004 G J Int"""

import numpy as np
from scipy.special import erf
# from mpmath import nsum, exp, inf, sinh
from model_1D import parameters as p
from model_1D import thermal as th
from model_1D import parameterised_interior as intr
import matplotlib.pyplot as plt
from scipy import integrate
from useful_and_bespoke import find_nearest_idx

# testing using Nimmo+ 2004 parameters - based on Table 4
core_params_default = {
    # constants
    'alpha_c': 1.35e-5,  # thermal expansion in 1/K
    'K_0': 634112127376.8201,  # compressibility at 0 pressure, back-calculated from Table 1
    'rho_0': 7019,  # density at 0 pressure
    'c_c': 840,
    'L_h': 750e3,  # latent heat melting iron in J/kg
    'k_c': 50,  # core thermal conductivity in W/m/K
    'gamma_melt': 9e-9,  # melting temp gradient at 350 GPa (~ICB) in K/Pa, Nimmo p 367
    'gamma_ad': 7.5e-9,  # adiabatic gradient at ICB in K/Pa, Nimmo p 368

    # radioisotope defaults using O'Neill+ 2020 SSR equation (1)
    'U_part': 0.0655,  # Faure+ 2020
    'Th_part': 0.0588,  # Faure+ 2020
    'K_part': 0.1558,  # Xiong + 2018
    'x_Eu': 1,  # concentration of r-process elements wrt solar (i.e. Eu, U, Th)

    # test variations?
    'chi_0': 4.2e-2,  # initial concentraton of light elements (mostly O?), Nimmo 2004 p 368
    'delta_rho': 590,  # density contrast at ICB - TODO: depends on chi_0 - see Alfe+ 2002
    'theta': 0.11,
    # +- 0.015, can use melting curve eqn to estimate knowing mag of temperature reduction from contaminants, mostly S & Si
    'rho_cen0': 12500,  # central density in kg/m3 - initial guess
}

test_params_default = {
    # 'p_cmb': 136e9,  # for testing - couple later
    # 'T_cmb': 4155,
    'Q_cmb': 9e12,
    # 'R_c': 3480e3,  # outer core radius
    # 'h_rad_c': 1.5e-12,
    # 'test_cooling_rate': -33 * p.sec2Gyr,
    # 'R_ic': 1220e3,
    # 'T_ic': 5581,
}


def temperature_profile(r, T_cen=None, D=None, **kwargs):
    return T_cen * np.exp(-r ** 2 / D ** 2)


def pressure_profile(r, L=None, p_cmb=None, R_c=None, rho_cen=None, **kwargs):
    """ p_c : pressure at CMB """
    # should this evolve?? no??
    # print('r', r, 'R_c', R_c, 'p_c', p_c, 'rho_cen', rho_cen, 'L', L)
    I = lambda x: (((3 * x ** 2 / 10) - (L ** 2 / 5)) * np.exp(-x ** 2 / L ** 2))
    return p_cmb + 4 * np.pi * p.G * rho_cen ** 2 / 3 * (I(R_c) - I(r))


def mass_profile(r, rmin=0, L=None, rho_cen=None, **kwargs):
    # print('in mass profile: r', r * 1e-3, 'km, ', 'rmin', rmin * 1e-3, 'km')
    I = lambda x: (-(L ** 2 / 2) * x * np.exp(-(x ** 2) / L ** 2) + L ** 3 / 4 * np.sqrt(np.pi) * erf(
        x / L))
    m = 4 * np.pi * rho_cen * (I(r) - I(rmin))
    # print('   m', m / p.M_E, 'M_E')
    return m


def gravity_profile(r, rho_cen=None, L=None, **kwargs):
    g = 4 * np.pi / 3 * p.G * rho_cen * r * (1 - (3 * r ** 2 / (5 * L ** 2)))
    return g


def density_profile(r, rho_cen=None, L=None, **kwargs):
    rho = rho_cen * np.exp(-r ** 2 / L ** 2)
    return rho


def gravpotential_profile(r, rho_cen=None, L=None, R_c=None, **kwargs):
    # gravitational potential energy, relative to 0 potential at the cmb
    I = lambda x: 2 / 3 * np.pi * p.G * rho_cen * x ** 2 * (1 - (3 * x ** 2 / (10 * L ** 2)))
    psi = I(r) - I(R_c)
    return psi


######      CLASS DEF     ######

class TerrestrialCore:

    def __init__(self, test=None, T_c0=None, r_res=200, **kwargs):

        # define default attributes
        default_attr = {
            # constants
            'age': 4.5,
            'alpha_c': 1.35e-5,  # thermal expansion in 1/K - should this avg value change a la noack parameterisation
            'K_0': 634112127376.8201,  # compressibility at 0 pressure, back-calculated from Table 1
            'rho_0': 7019,  # density at 0 pressure
            'c_c': 840,  # should this avg value change a la noack parameterisation
            'L_h': 750e3,  # latent heat melting iron in J/kg
            'k_c': 50,  # core thermal conductivity in W/m/K
            'gamma_melt': 9e-9,  # melting temp gradient at 350 GPa (~ICB) in K/Pa, Nimmo p 367
            'gamma_ad': 7.5e-9,  # adiabatic gradient at ICB in K/Pa, Nimmo p 368
            'delta': 0.5,  # liquid-solid light element partitioning
            'C_r': -13.454545e3,  # ratio of dT_i/dt and dR_i/dt - todo: vary??

            # radioisotope defaults using O'Neill+ 2020 SSR equation (1)
            'U_part': 0.0655,  # Faure+ 2020
            'Th_part': 0.0588,  # Faure+ 2020
            'K_part': 0.1558,  # Xiong + 2018
            'x_Eu': 1,  # concentration of r-process elements wrt solar (i.e. Eu, U, Th)

            # test variations?
            'chi_0': 4.2e-2,  # initial concentraton of light elements (mostly O?), Nimmo 2004 p 368
            'delta_rho': 590,  # density contrast at ICB - TODO: depends on chi_0 - see Alfe+ 2002
            'theta': 0.11,
            # +- 0.015, can use melting curve eqn to estimate knowing mag of temperature reduction from contaminants, mostly S & Si
            # 'R_c': 3480e3,  # outer core radius
            'rho_cen0': 12500,  # central density in kg/m3- initial guess
        }

        if test is not None and 'Nimmo' in test:
            self.CMF = 0.32264755186184685
            self.rho_cen = 12500
            self.p_cmb = 139e9  # 136e9  # for testing - couple later
            self.Q_cmb = 9e12
            self.R_c = 3480e3  # outer core radius
            self.h_rad_c = 1.5e-12
            # test_cooling_rate = -33 * p.sec2Gyr
            self.R_ic = 1220e3
            self.T_ic = 5581
            self.K_part = 400e-6 / 0.00025982905982905983  # 400 ppm K in core
            self.U_part = 0
            self.Th_part = 0
        if test == 'Nimmo_evol':
            self.T_cmb = 4800
        elif test == 'Nimmo_static':
            self.T_cmb = 4155


        elif test is None:
            # get initial conditions at CMB - spcifically p_cmb and R_c
            pre = intr.InteriorStructure()
            pre.init_structure(which='hot', T_c0=T_c0, **kwargs)
            print('M_p in init structure', pre.M_p / p.M_E, 'M_E | p_cmb', pre.p_cmb * 1e-9, 'GPa | T_cmb', pre.T_cmb,
                  'K')
            self.__dict__.update(**pre.__dict__)

        # add other input parameters, override if given
        default_attr.update(kwargs)
        for k, v in default_attr.items():
            if not hasattr(self, k):
                # print('using', k, v, 'from default_attr')
                self.__dict__.update({k: v})
            # else:
            # print('using', k, 'set above')
        self.rv_core = np.linspace(0, self.R_c, num=r_res)

        # add derived parameters - that won't evolve
        self.L = np.sqrt((3 * self.K_0 * np.log(self.rho_cen0 / self.rho_0 + 1)) / (
                2 * np.pi * p.G * self.rho_0 * self.rho_cen0))  # initial guess
        if test is None:
            # find central density to match input CMF
            self.rho_cen = self.get_central_density(m1=self.CMF * p.M_E, **kwargs)
            self.L = np.sqrt((3 * self.K_0 * np.log(self.rho_cen / self.rho_0 + 1)) / (
                    2 * np.pi * p.G * self.rho_0 * self.rho_cen))  # length scale
        self.D = np.sqrt(3 * self.c_c / (2 * np.pi * self.alpha_c * self.rho_cen * p.G))  # another length scale
        self.M_c = mass_profile(r=self.R_c, L=self.L, rho_cen=self.rho_cen)  # mass of core

        print('calculated core mass =', self.M_c / p.M_E, 'M_E')
        if test is None:
            print('input from CMF =', self.CMF * pre.M_p / p.M_E,
                  'M_E')
        self.pv_core = pressure_profile(r=self.rv_core, L=self.L, rho_cen=self.rho_cen, R_c=self.R_c,
                                        p_cmb=self.p_cmb)

        # initial conditions
        self.dT_cdt = 0
        self.dRdt = 0
        self.evolve_profiles(rho_cen=self.rho_cen, **kwargs)

    def melt_curve(self, p=None, T_m0=1695, T_m1=10.9e-12, T_m2=-8.0e-24, which='nimmo', **kwargs):
        """ theta depends especially on concentration of S and Si (reduces T at ICB), but other data are for pure Fe
        uncertainty on Tm1 and Tm2 is about 1.1e-12 Pa^-1 and 1.5e-24 Pa^-2
        stixrude parameterisation depends directly on concentration of light elements (can assume constant if mostly S/Si contributing)"""
        if which == 'stixrude':
            # core melt temperature parameterisation from Stixrude (2014)
            return 6500 * (p * 1e-9 / 340) ** 0.515 / (1 - np.log(1 - self.chi))
        elif which == 'nimmo':
            return T_m0 * (1 - self.theta) * (1 + T_m1 * p + T_m2 * p ** 2)

    def adiabat(self, p=None, T_a1=3.5e-12, T_a2=-1.8e-24, **kwargs):
        """ adiabatic temperature as a function of pressure determined by NImmo from fitting adiabatic T(r) profile and relating z to pressure
        therefore constants T_a1 etc may vary with different central densities rho_cen and T_cen ?? test
        T_c, p_c : temperature and pressure at CMB"""

        # try:
        return self.T_cmb * (1 + T_a1 * p + T_a2 * p ** 2) / (1 + T_a1 * self.p_cmb + T_a2 * self.p_cmb ** 2)
        # except ValueError:
        #     p = np.expand_dims(p, axis=0)
        #     return self.T_cmb * (1 + T_a1 * p + T_a2 * p ** 2) / (1 + T_a1 * self.p_cmb + T_a2 * self.p_cmb ** 2)

    def get_inner_core_wrapper(self, plot='error', **kwargs):
        # todo: fit to adiabat equation for different T, p at outer core radius

        Tm = self.melt_curve(p=self.pv_core, **kwargs)

        # solve adiabatic T(p) given T_c
        if np.size(self.T_cmb) == 1:
            Ta_p = self.adiabat(p=self.pv_core, **kwargs)
            Ta_z = self.Tv_core
            R_ic, T_ic = self.get_inner_core(Ta_p=Ta_p, Ta_z=Ta_z, Tm=Tm, T_cen=self.T_cen,
                                             T_cmb=self.T_cmb, plot=plot, **kwargs)

        else:
            print('\ngetting inner core w full solution, T_cmb', np.shape(self.T_cmb))
            # if self.T_cmb.ndim == 1:
            #     # patch
            #     self.T_cmb = np.expand_dims(self.T_cmb, axis=0)  # patch
            #     self.T_cen = np.expand_dims(self.T_cen, axis=0)  # patch
            #     self.pv_core = np.expand_dims(self.pv_core, axis=1)  # patch
            #     self.rv_core = np.expand_dims(self.rv_core, axis=1)  # patch
            #     Tm = np.expand_dims(Tm, axis=1)
            #     print('T_cmb', np.shape(self.T_cmb))
            #     print('pv_core', np.shape(self.pv_core))
            #     print('T_cen', np.shape(self.T_cen))
            #     print('Tm', np.shape(Tm))

            # using full solution for T_cmb
            T_ic = np.zeros(np.shape(self.T_cmb))
            R_ic = np.zeros(np.shape(self.T_cmb))
            Ta_p = self.adiabat(p=self.pv_core, **kwargs)
            for ii, ti in enumerate(self.T_cmb):
                # print('T cen shape in inner core wrapper', np.shape(ti))
                Ta_z = temperature_profile(self.rv_core, T_cen=self.T_cen[ii], D=self.D)
                # print('Ta_p[ii]', np.shape(Ta_p[ii]))
                R_ic_ii, T_ic_ii = self.get_inner_core(Ta_p=Ta_p[ii], Ta_z=Ta_z, Tm=Tm, T_cen=self.T_cen[ii],
                                                       T_cmb=self.T_cmb[ii], plot=plot,
                                                       **kwargs)
                R_ic[ii] = R_ic_ii
                T_ic[ii] = T_ic_ii

        if plot is True:
            show_inner_core(pv=self.pv_core, Ta_p=Ta_p, Tm=Tm, Ta_z=Ta_z, r=self.rv_core, c=self, **kwargs)
        return R_ic, T_ic

    def get_inner_core(self, Ta_p=None, Ta_z=None, Tm=None, T_cen=None, T_cmb=None, plot='error', **kwargs):
        idx_p = np.argwhere(np.diff(np.sign(Ta_p - Tm))).flatten()
        idx_z = np.argwhere(np.diff(np.sign(Ta_z - Tm))).flatten()

        # select which adiabat
        idx = idx_z
        Ta = Ta_z

        try:
            if np.size(idx) == 1:
                T_ic = Ta[idx].item()
                R_ic = self.rv_core[idx].item()
            elif np.size(idx) > 1:
                # use last (shallowest) occurrence
                # show_inner_core(p=self.pv_core, Ta_p=Ta_p, Tm=Tm, Ta_z=Ta_z, r=self.rv_core,
                #                 title='mulitple intersections', **kwargs)
                T_ic = Ta[idx[-1]].item()
                R_ic = self.rv_core[idx[-1]].item()
            elif np.less(Ta, Tm).all():
                # entirely below solidus
                T_ic = T_cmb
                R_ic = self.R_c
            else:
                # no solid inner core
                T_ic = T_cen
                R_ic = 0
        except ValueError as e:
            if plot == 'error':
                show_inner_core(pv=self.pv_core, Ta_p=Ta_p, Tm=Tm, Ta_z=Ta_z, r=self.rv_core, **kwargs)
            raise e
        return R_ic, T_ic

    def inner_core_age(self, t, **kwargs):
        if np.size(self.R_ic) <= 1:
            raise Exception('error in inner core age: only one radius given')
        idx = np.argmax(self.R_ic > 0)
        # print('R_ic[idx]', self.R_ic[idx]*1e-3, 'km', 'R_ic[idx-1]', self.R_ic[idx-1]*1e-3, 'km', 'idx', idx, )
        return t[idx]

    def grav_heat(self, **kwargs):
        if (np.size(self.R_ic) == 1) and (self.R_ic == 0 or self.R_ic == self.R_c):
            self.Q_g_bar = 0
        else:
            # gravitational heat - change in energy associated with the separation of the light component
            C_c = 4 * np.pi * self.R_ic ** 2 * self.rho_ic * self.chi / self.M_oc
            self.beta_c = 1 / self.chi * self.delta_rho / self.rho_ic  # compositional expansion coefficient - think rho here is rho_ic but maybe rho at cmb? NImmo eq 21

            # analytic nimmo way
            C_g2 = 3 * self.L ** 2 / 16 - self.R_c ** 2 / 2 * (1 - (3 * self.R_c ** 2 / (10 * self.L ** 2)))
            f_g = lambda x: ((3 / 20 * x ** 5 - self.L ** 2 / 8 * x ** 3 - self.L ** 2 * C_g2 * x) * np.exp(
                -x ** 2 / self.L ** 2) + C_g2 / 2 * self.L ** 3 * np.sqrt(np.pi) * erf(x / self.L))
            I_g = 8 * np.pi ** 2 * self.rho_cen ** 2 * p.G / 3 * (f_g(self.R_c) - f_g(self.R_ic))
            self.Q_g_bar = (I_g - self.M_oc * self.psi_ic) * self.beta_c * C_c * self.C_r

            if np.size(self.R_ic) > 1:
                self.Q_g_bar = np.squeeze(self.Q_g_bar)

        if np.size(self.R_ic) > 1:
            idx = np.where(self.R_ic == self.R_c)[0]
            self.Q_g_bar[idx] = 0
        return self.Q_g_bar

    def latent_heat(self, **kwargs):
        # latent heat
        if np.size(self.R_ic) == 1 and (self.R_ic == self.R_c):
            return 0
        Q_L_bar = 4 * np.pi * self.R_ic ** 2 * self.L_h * self.rho_ic * self.C_r
        if np.size(self.R_ic) > 1:
            Q_L_bar = np.squeeze(Q_L_bar)
            idx = np.where(self.R_ic == self.R_c)[0]
            Q_L_bar[idx] = 0
        self.Q_l_bar = Q_L_bar
        return self.Q_l_bar

    def specific_heat(self, **kwargs):
        # specific heat
        A = np.sqrt((1 / self.L ** 2 + 1 / self.D ** 2) ** -1)
        I_s = 4 * np.pi * self.T_cen * self.rho_cen * (
                -(A ** 2) * self.R_c / 2 * np.exp(-(self.R_c ** 2) / A ** 2) + A ** 3 * np.sqrt(np.pi) / 4 * erf(
            self.R_c / A))
        Q_s_bar = -self.c_c / self.T_cmb * I_s  # Q_s_bar = Q_s / dT_c/dtrho_ic
        if np.size(self.T_cmb) > 1:
            Q_s_bar = np.squeeze(Q_s_bar)
        self.Q_s_bar = Q_s_bar
        return self.Q_s_bar

    def radioactive_heat(self, t=None,  # h_rad_c=None,
                         tau_i=np.array([1250, 4468, 703.8, 14050]),  # half life in Myr
                         h_i=np.array([28.761e-6, 94.946e-6, 568.402e-6, 26.368e-6]),  # heat production in W/kg
                         c_i=np.array([30.4e-9, 22.7e-9, 0.16e-9, 85e-9]),  # BSE concentration in kg/kg
                         # part_i=np.array([0.1158, 0.0655, 0.0655, 0.0588]),
                         **kwargs):
        # negligible? Faure 2020

        # radioactive heat
        # if D >= L:
        #     B = np.sqrt((1 / L ** 2 - 1 / D ** 2) ** -1)
        #     I_T = 4 * np.pi * rho_cen / T_cen * (
        #                 -B ** 2 * R_c / 2 * np.exp(-R_c ** 2 / B ** 2) + B ** 3 * np.sqrt(np.pi) / 4 * erf(R_c / B))
        # else:
        #     B = np.sqrt((1 / D ** 2 - 1 / L ** 2) ** -1)
        #     S_n = R_c / (2 * np.sqrt(np.pi)) + B / np.sqrt(np.pi) * nsum(lambda x: exp(-x ** 2 / 4) / x * sinh(x * R_c / B),
        #                                                                  [1, inf])
        #     I_T = 4 * np.pi * rho_cen / T_cen * (
        #             -B ** 2 * R_c / 2 * np.exp(-R_c ** 2 / B ** 2) - B ** 2 * S_n / 2)
        part_i = [self.K_part, self.U_part, self.U_part, self.Th_part]
        h_rad_c = th.h_rad(t, c_i * part_i, h_i, tau_i, self.age, x_Eu=self.x_Eu)
        # part_i is the partitioning coeffiient between BSE and core

        self.Q_r = self.M_c * h_rad_c
        # print('   in rad heat: M_c', M_c, 'h_rad_c', h_rad_c )
        return self.Q_r

    def evolve_profiles(self, **kwargs):
        if np.size(self.T_cmb) > 1:
            print('using full soluton for T_cmb with shape', np.shape(self.T_cmb))
            self.T_cmb = np.expand_dims(self.T_cmb, axis=0).T
            print(' new T_cmb shape', np.shape(self.T_cmb))
        # else:
        self.T_cen = self.T_cmb / np.exp(-(self.R_c ** 2) / (self.D ** 2))
        # print('T_cen shape', np.shape(self.T_cen))

        # update profiles in core
        try:
            self.Tv_core = temperature_profile(r=self.rv_core, T_cen=self.T_cen, D=self.D)
        except ValueError:
            print('catching valueerror - updating profiles with full time evol')
            # print('tile T cen shape', np.shape(np.tile(self.T_cen,(len(self.rv_core), 1))))
            self.Tv_core = temperature_profile(r=self.rv_core, T_cen=self.T_cen, D=self.D)
        # print('  Tv_core solution shape', np.shape(self.Tv_core))

        # update inner core boundary
        self.R_ic, self.T_ic = self.get_inner_core_wrapper()

        self.M_oc = mass_profile(r=self.R_c, rmin=self.R_ic, L=self.L, rho_cen=self.rho_cen)  # mass of outer core
        self.rho_ic = density_profile(r=self.R_ic, rho_cen=self.rho_cen, L=self.L)
        self.g_ic = gravity_profile(r=self.R_ic, rho_cen=self.rho_cen, L=self.L)
        self.psi_ic = gravpotential_profile(r=self.R_ic, R_c=self.R_c, rho_cen=self.rho_cen, L=self.L)

        # update compositions
        if self.C_r is None:
            self.C_r = 1 / (self.gamma_melt - self.gamma_ad) * self.T_ic / self.T_c / (
                    self.rho_ic * self.g_ic)  # Gubbins 2003 eq 36 relating dR_i/dt to dT/dt
        self.chi = self.chi_0 / (1 - (self.R_ic / self.R_c) ** 3 + self.delta * (self.R_ic / self.R_c) ** 3)

    def get_central_density(self, m1=None, tol=1e-9, max_iter=1000, **kwargs):
        """ iterate rho_cen until mass profile produces the desired m1 value from CMF parameterisation. tol in kg"""
        rho_cen0 = self.rho_cen0  # initial guess
        m0 = mass_profile(r=self.R_c, L=self.L, rho_cen=rho_cen0)  # current mass of core
        n = 0
        # print('initial error:', abs(m1 - m0)/m1)
        while (tol < abs(m1 - m0) / m1) and (n < max_iter):
            n = n + 1
            rho_cen1 = rho_cen0 * m1 / m0
            L = np.sqrt((3 * self.K_0 * np.log(rho_cen1 / self.rho_0 + 1)) / (
                    2 * np.pi * p.G * self.rho_0 * rho_cen1))  # length scale
            # D = np.sqrt(3 * c.C_p / (2 * np.pi * c.alpha_c * rho_cen1 * p.G))  # another length scale
            m0 = mass_profile(r=self.R_c, rho_cen=rho_cen1, L=L)  # mass of core
            # print('guess', n, ':', rho_cen1, 'kg/m3', 'rho1/rho0:', rho_cen1/rho_cen0, 'm1/m0', m1/m0)
            rho_cen0 = rho_cen1

        return rho_cen0

    def solve(self, t0=0, tf=4.5, plot=True, **kwargs):
        if np.size(self.T_cmb == 1):
            T_c0 = self.T_cmb  # initialise
        else:
            raise Exception('cannot solve core temperature evolution starting from non-scalar T_cmb')
        if not hasattr(self, 'Q_cmb'):
            self.Q_cmb = test_params_default['Q_cmb']
            print('using fixed input Q_cmb =', self.Q_cmb * 1e-12, 'TW')
        print('T_cmb0', T_c0)
        f = integrate.solve_ivp(fun=lambda t, y: temperature_LHS(t, y, **dict(c=self, **kwargs)),
                                t_span=(t0 * 1e9 * p.years2sec, tf * 1e9 * p.years2sec), y0=T_c0,
                                max_step=100e6 * p.years2sec,
                                method='RK45',
                                )

        print('>>>>finished core integration >>>')

        # update profiles with T(t) solution
        self.T_cmb = f.y[0]
        self.R_ic = f.y[1]
        self.evolve_profiles(t=f.t)

        # patch array shapes
        self.T_cmb = np.squeeze(self.T_cmb)
        self.T_cen = np.squeeze(self.T_cen)

        if plot:
            plot_flux_evol(c=self, t=f.t, **kwargs)

    def core_gradient(self, test_cooling_rate=None, t=None, Q_r=None, **kwargs):
        if Q_r is None:
            self.Q_r = self.radioactive_heat(t=t, )
        else:
            self.Q_r = Q_r
        # print('in core_gradient() C_r', self.C_r)
        self.Q_s_bar = self.specific_heat(**kwargs)
        self.Q_l_bar = self.latent_heat(**kwargs)
        self.Q_g_bar = self.grav_heat(**kwargs)
        if test_cooling_rate is not None:
            print('Q_s =', self.Q_s_bar * test_cooling_rate * 1e-12, 'TW')
            print('Q_l =', self.Q_l_bar * test_cooling_rate * 1e-12, 'TW')
            print('Q_g =', self.Q_g_bar * test_cooling_rate * 1e-12, 'TW')
            print('Q_r =', self.Q_r * 1e-12, 'TW')

        dTdt = (self.Q_r - self.Q_cmb) / -(
                self.Q_l_bar + self.Q_g_bar + self.Q_s_bar)  # note typo in NImmo paper where Q bar should be -ve

        self.Q_s = self.specific_heat() * dTdt
        self.Q_l = self.latent_heat() * dTdt
        self.Q_g = self.grav_heat() * dTdt
        self.dTdt = dTdt
        return dTdt


def K_0(L=7272e3, rho_0=7019, rho_cen=12500):
    # compressibility at 0 pressure -> let this be fixed at Nimmo value (todo: iterate?)
    return L ** 2 * 2 * np.pi * p.G * rho_0 * rho_cen / (3 * np.log(rho_cen / rho_0 + 1))


def temperature_LHS(t, y, c=None, **kwargs):
    """ ODE equation to solve """
    # print('  LHS in core: t =', t*p.sec2Gyr, 'Gyr', 'T_cmb', c.T_cmb, 'K')
    c.T_cmb = y
    # try:
    #     print('in temperature LHS', 'dRdt', c.dRdt, 'dTdt', c.dT_cdt)
    #     c.C_r = c.dRdt / c.dT_cdt  # update
    # except ZeroDivisionError:
    #     c.C_r = 0

    c.evolve_profiles(**kwargs)
    c.dT_cdt = c.core_gradient(t=t, test_cooling_rate=None, **kwargs)

    # if c.R_ic == 0:
    #     c.dRdt = 0
    # else:
    #     c.dRdt = c.Q_l / (4 * np.pi * c.R_ic**2 * c.L_h * c.rho_ic)

    return c.dT_cdt  # c.dRdt


# def core_heat_flow(c,  **kwargs):
#     """
#     L_h: latent heating
#     delta_rho: density difference across inner core boundary due to presence of light element in outer core (measured for earth??)
#     gamma_melt: dT/dp of melting for core
#     gamma_ad: adiabatic dT/dp fore core
#     chi_0: initial concentration of light element in outer core (expelled from inner core during freezing)
#     delta: mass partition coefficient of light element between solid and liquid - results should be insensitive
#     """
#
#     c.evolve_profiles(**kwargs)
#
#     # calculate heat fluxes in governing equation, Nimmo+ 2004 eq 30, where Q_bar = Q / dT_c/dt
#     c.dTdt = core_gradient(c, **kwargs)
#     return c


def print_core(c, ii=-1):
    print('\n----parameters')
    print(' R_c = ', c.R_c * 1e-3, 'km')
    print(' L = ', c.L * 1e-3, 'km')
    print(' D = ', c.D * 1e-3, 'km')
    print(' p_cmb =', c.p_cmb * 1e-9, 'GPa')
    print(' C_r = ', c.C_r)
    # print(' p_cen =', c.p_cen*1e-9, 'GPa')
    print(' rho_cen =', c.rho_cen, 'kg/m3')

    print('\n----timestep')
    if np.size(c.T_cmb) > 1:
        print(' T_cmb = ', np.squeeze(np.squeeze(c.T_cmb))[ii], 'K')
        print(' R_ic = ', np.squeeze(np.squeeze(c.R_ic))[ii] * 1e-3, 'km')
        print(' T_ic =', np.squeeze(np.squeeze(c.T_ic))[ii], 'K')
        print(' chi = ', np.squeeze(np.squeeze(c.chi))[ii])
        print(' beta_c =', np.squeeze(np.squeeze(c.beta_c))[ii])
        print(' dT_c/dt =', np.squeeze(np.squeeze(c.dTdt))[ii] / p.sec2Gyr, 'K/Gyr')
        print(' Q_s', np.squeeze(np.squeeze(c.Q_s))[ii] * 1e-12, 'TW')
        print(' Q_L', np.squeeze(np.squeeze(c.Q_l))[ii] * 1e-12, 'TW')
        print(' Q_g', np.squeeze(np.squeeze(c.Q_g))[ii] * 1e-12, 'TW')
        print(' Q_R', np.squeeze(np.squeeze(c.Q_r))[ii] * 1e-12, 'TW')
    else:
        print(' T_cmb = ', c.T_cmb, 'K')
        print(' R_ic = ', c.R_ic * 1e-3, 'km')
        print(' T_ic =', c.T_ic, 'K')
        print(' chi = ', c.chi)
        print(' beta_c =', c.beta_c)
        print(' dT_c/dt =', c.dTdt / p.sec2Gyr, 'K/Gyr')
        print(' Q_s', c.Q_s * 1e-12, 'TW')
        print(' Q_L', c.Q_l * 1e-12, 'TW')
        print(' Q_g', c.Q_g * 1e-12, 'TW')
        print(' Q_R', c.Q_r * 1e-12, 'TW')
    print('\n')


def show_inner_core(pv=None, Ta_p=None, Tm=None, Ta_z=None, r=None, c=None, title=None, t_plot=None, times=None,
                    **kwargs):
    if pv is None:
        # calculate from object
        if t_plot is None:
            # final timestep
            ii = -1
        else:
            ii = find_nearest_idx(times, t_plot)
        pv = c.pv_core
        r = c.rv_core
        Ta_z = temperature_profile(c.rv_core, T_cen=c.T_cen[ii], D=c.D)
        Tm = c.melt_curve(p=pv, **kwargs)

        T_cmb_holder = c.T_cmb
        if np.size(T_cmb_holder) > 1:
            c.T_cmb = T_cmb_holder[-1]
        Ta_p = c.adiabat(p=np.expand_dims(pv, axis=0), **kwargs).T
        # print('T_cmb', c.T_cmb)
        c.T_cmb = T_cmb_holder

    fig, (ax, ax2) = plt.subplots(1, 2)
    if Ta_z is not None:
        ax2.plot(Ta_z, r * 1e-3, 'r--', label='z-adiabat')
        ax2.plot(Ta_p, r * 1e-3, 'k--', label='p-adiabat')
        ax2.plot(Tm, r * 1e-3, 'k-', label='solidus')

        idx2 = np.argwhere(np.diff(np.sign(Ta_z - Tm))).flatten()
        if np.size(idx2) == 1:
            ax2.axhline(r[idx2] * 1e-3, c='r', lw=0.5)
            ax.axhline(pv[idx2] * 1e-9, c='r', lw=0.5)

    if Ta_p is not None:
        ax.plot(Ta_z, pv * 1e-9, 'r--', label='z-adiabat')
        ax.plot(Ta_p, pv * 1e-9, 'k--', label='p-adiabat')
        ax.plot(Tm, pv * 1e-9, 'k-', label='solidus')

        idx = np.argwhere(np.diff(np.sign(Ta_p - Tm))).flatten()  # find intersection in T-p axes
        if np.size(idx) == 1:
            ax.axhline(pv[idx] * 1e-9, c='k', lw=0.5)
            ax2.axhline(r[idx] * 1e-3, c='k', lw=0.5)
        # else:
        #     print('idx', idx)

    ax.invert_yaxis()
    ax.set_xlabel('T (K)')
    ax.set_ylabel('P (GPa)')
    ax.legend()

    ax2.set_ylabel('r (km)')
    ax2.set_xlabel('T (K)')
    if title is None:
        title = str(times[ii] * p.sec2Gyr) + ' Gyr'
    plt.suptitle(title)

    # if c is not None:
    #     print('in core: p_cmb', c.p_cmb*1e-9, 'GPa | rho_cen', c.rho_cen, '| R_c', c.R_c*1e-3, 'km | T_cmb', c.T_cmb, 'K | T_cen', c.T_cen, 'K')
    # plt.show()


def plot_TP(p=None, c=None, T_ic=None, **kwargs):
    if p is None:
        p = np.linspace(136e9, 365e9)
    Tp = c.adiabat(p, **kwargs)
    Tm = c.melt_curve(p, **kwargs)

    plt.figure()
    plt.plot(p, Tp, c='k', ls='-', label='adiabat')
    plt.plot(p, Tm, c='k', ls='--', label='solidus')
    plt.axhline(c.T_cmb, lw=0.5, c='r', label='CMB temperature or pressure')
    plt.axvline(c.p_cmb, lw=0.5, c='r')
    if T_ic is not None:
        plt.axhline(T_ic, lw=0.5, c='b', label='test ICB temperature')
    plt.xlabel('pressure (Pa)')
    plt.ylabel('temperature (K)')
    plt.legend()
    plt.show()


def plot_density_profile(r=None, c=None, **kwargs):
    import matplotlib.pyplot as plt

    if r is None:
        r = c.rv
    plt.figure()
    c.rho_ = c.density_profile(r=r)
    c.rho_ic = c.density_profile(r=c.R_ic)


def plot_mass_profile(r=None, c=None, **kwargs):
    # update profiles

    import matplotlib.pyplot as plt

    plt.figure()
    if r is None:
        r = np.linspace(0, c.R_c)
    mass = c.mass_profile(r)
    plt.plot(mass / p.M_E, r * 1e-3, c='k', lw=0.5, label='analytic profile')
    plt.axhline(c.R_c * 1e-3, c='k', lw=0.5, ls='--')
    plt.axvline(mass[-1] / p.M_E, c='k', lw=1, ls='--', label='analytic core mass')
    plt.axvline(c.CMF * c.M_p / p.M_E, c='r', lw=0.5, label='required mass from CMF')
    plt.xlabel('mass (M_E)')
    plt.ylabel('radius (km)')
    print('central density:', c.rho_cen, 'kg/m3')
    print('mass error:', mass[-1] - c.CMF * c.M_p)
    plt.legend()
    plt.show()


def plot_psi_profile(r=None, c=None, **kwargs):
    # update profiles

    import matplotlib.pyplot as plt

    plt.figure()
    if r is None:
        r = np.linspace(0, c.R_c)
    psi = c.gravpotential_profile(r)
    psi_icb = c.gravpotential_profile(c.R_ic)
    print('psi_icb * M_oc', psi_icb * c.M_oc)
    plt.plot(psi, r * 1e-3, c='k', lw=0.5, label='analytic profile')
    plt.axhline(c.R_c * 1e-3, c='k', lw=0.5, ls='-')
    plt.axvline(c.R_c * 1e-3, c='k', lw=0.5, ls='-')
    plt.axvline(psi_icb, c='k', lw=1, ls='--', label='ICB')
    plt.axhline(c.R_ic * 1e-3, c='k', lw=1, ls='--')
    plt.xlabel('gravitational potential (J)')
    plt.ylabel('radius (km)')
    plt.legend()
    # plt.show()


def plot_flux_dependence(x_vec=None, x_name='CMF', core_params=None,
                         flux_colours=['xkcd:pink', 'xkcd:chartreuse', 'xkcd:peach', 'xkcd:cerulean', 'k'],
                         **kwargs):
    # for fixed cooling rate how do other fluxes change w/ x

    if core_params is None:
        # testing using Nimmo+ 2004 parameters - based on Table 4
        core_params = core_params_default
        core_params.update(test_params_default)

    fig, axes = plt.subplots(1, 1)

    if x_vec is None:
        x_vec = np.linspace(0.1, 0.5, num=8)

    for ii, x in enumerate(x_vec):
        print('\n\n---------------------------->', ii, '| x', x)
        core_params.update({x_name: x})

        # calculate heat fluxes
        c = TerrestrialCore(**core_params)
        c.solve(plot=False, **kwargs)

        if np.size(c.Q_s) == 1:
            Q_s = c.Q_s
            Q_l = c.Q_l
            Q_g = c.Q_g
            Q_r = c.Q_r
        else:
            Q_s = c.Q_s[-1]
            Q_l = c.Q_l[-1]
            Q_g = c.Q_g[-1]
            Q_r = c.Q_r[-1]
        total_flux = Q_s + Q_l + Q_g + Q_r

        fluxes = np.array([Q_s, Q_l, Q_g, Q_r, total_flux]) * 1e-12
        labels = ['specific heat', 'latent heat', 'grav. heat', 'rad. heat', 'total']
        for jj in range(len(fluxes)):
            if ii == 0:
                label = labels[jj]
            else:
                label = None
            axes.plot(x, fluxes[jj], marker='.', c=flux_colours[jj], label=label)
        # print('R_ic', c.R_ic)

    plt.xlabel(x_name)
    plt.ylabel('flux (TW)')
    axes.legend()
    plt.show()


def plot_flux_evol(c=None, t=None, pl=None, **kwargs):
    fig, axes = plt.subplots(3, 1)
    axes[0].plot(t * p.sec2Gyr, c.T_cmb, 'k-', label='$T_{cmb}$')
    axes[1].plot(t * p.sec2Gyr, c.Q_cmb * 1e-12, label='$Q_{cmb}$')
    axes[1].plot(t * p.sec2Gyr, c.Q_r * 1e-12, label='$Q_{r}$')
    axes[1].plot(t * p.sec2Gyr, c.Q_s * 1e-12, label='$Q_{s}$')
    axes[1].plot(t * p.sec2Gyr, c.Q_l * 1e-12, label='$Q_{l}$')
    axes[1].plot(t * p.sec2Gyr, c.Q_g * 1e-12, label='$Q_{g}$')
    # axes[2].plot(t * p.sec2Gyr, c.dTdt / p.sec2Gyr, 'k-')
    axes[2].plot(t * p.sec2Gyr, c.R_ic * 1e-3, 'k-')
    axes[2].axhline(c.R_c * 1e-3, c='k', ls='--', lw=0.5)
    axes[2].set_ylabel('$R_{IC}$ (km)')
    axes[1].set_ylabel('$Q$ (TW)')
    axes[-1].set_xlabel('t (Gyr)')

    if pl is not None:
        axes[0].plot(t * p.sec2Gyr, pl.T_m, c='xkcd:grey', ls='-', lw=0.5, label='$T_m$')
        axes[0].plot(t * p.sec2Gyr, pl.Tp, c='xkcd:grey', ls='--', lw=0.5, label='$T_m^\prime$')
        axes[0].legend()
        axes[1].plot(t * p.sec2Gyr, pl.Q_ubl * 1e-12, ls='--', label='$Q_{ubl}$')
        fig.suptitle('core + mantle evolution')
        axes[0].set_ylabel('$T$ (K)')
    else:
        fig.suptitle('core evolution')
        axes[0].set_ylabel('$T_{cmb}$ (K)')
    axes[1].legend()
    # plt.show()
