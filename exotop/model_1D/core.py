import numpy as np
from scipy.special import erf
from mpmath import nsum, exp, inf, sinh
from model_1D import parameters as p
from model_1D import thermal as th
from model_1D import noack_interior as intr
import matplotlib.pyplot as plt
from scipy import integrate

# testing using Nimmo+ 2004 parameters - based on Table 4
core_params_default = {
    # constants
    'alpha_c': 1.35e-5,  # thermal expansion in 1/K
    'K_0': 634112127376.8201,  # compressibility at 0 pressure, back-calculated from Table 1
    'rho_0': 7019,  # density at 0 pressure
    'C_p': 840,
    'L_h': 750e3,  # latent heat melting iron in J/kg
    'k_c': 50,  # core thermal conductivity in W/m/K
    'gamma_melt': 9e-9,  # melting temp gradient at 350 GPa (~ICB) in K/Pa, Nimmo p 367
    'gamma_ad': 7.5e-9,  # adiabatic gradient at ICB in K/Pa, Nimmo p 368

    # radioisotope defaults using O'Neill+ 2020 SSR equation (1)
    'U_part':0.0655,
    'Th_part': 0.0588,
    'K_part':0.1558,
    'x_Eu': 1,  # concentration of r-process elements wrt solar (i.e. Eu, U, Th)

    # test variations?
    'chi_0': 4.2e-2,  # initial concentraton of light elements (mostly O?), Nimmo 2004 p 368
    'delta_rho': 590,  # density contrast at ICB - TODO: depends on chi_0 - see Alfe+ 2002
    'theta': 0.11,
    # +- 0.015, can use melting curve eqn to estimate knowing mag of temperature reduction from contaminants, mostly S & Si
    # 'R_c': 3480e3,  # outer core radius
    'rho_cen': 12500,  # central density in kg/m3 - initial guess
}

test_params_default = {
    # 'p_cmb': 136e9,  # for testing - couple later
    # 'T_cmb': 4155,
    'Q_cmb': 9e12,
    'T_a1': 3.5e-12,  # for core adiabat, fitted parameters
    'T_a2': -1.8e-24,
    'h_rad_c': 1.5e-12,
    # 'test_cooling_rate': -33 * p.sec2Gyr,
    # 'R_ic': 1220e3,
    # 'T_ic': 5581,
    'C_r': -13.45e3,
}


######      CLASS DEF     ######
class TerrestrialCore():

    def __init__(self, cmf_in=None, full_setup=True, r_res=200, **kwargs):

        # define default attributes
        default_attr = {
            # constants
            'alpha_c': 1.35e-5,  # thermal expansion in 1/K - should this avg value change a la noack parameterisation
            'K_0': 634112127376.8201,  # compressibility at 0 pressure, back-calculated from Table 1
            'rho_0': 7019,  # density at 0 pressure
            'C_p': 840,  # should this avg value change a la noack parameterisation
            'L_h': 750e3,  # latent heat melting iron in J/kg
            'k_c': 50,  # core thermal conductivity in W/m/K
            'gamma_melt': 9e-9,  # melting temp gradient at 350 GPa (~ICB) in K/Pa, Nimmo p 367
            'gamma_ad': 7.5e-9,  # adiabatic gradient at ICB in K/Pa, Nimmo p 368
            'delta': 0.5,  #liquid-solid light element partitioning

            # test variations?
            'chi_0': 4.2e-2,  # initial concentraton of light elements (mostly O?), Nimmo 2004 p 368
            'delta_rho': 590,  # density contrast at ICB - TODO: depends on chi_0 - see Alfe+ 2002
            'theta': 0.11,
            # +- 0.015, can use melting curve eqn to estimate knowing mag of temperature reduction from contaminants, mostly S & Si
            'R_c': 3480e3,  # outer core radius
            'rho_cen': 12500,  # central density in kg/m3
        }

        # add input parameters, use default if not given
        default_attr.update(kwargs)
        self.__dict__.update((k, v) for k, v in default_attr.items())

        # get initial conditions at CMB
        pre = intr.InteriorStructure()
        pre.init_structure(which='hot', CMF=cmf_in)
        self.__dict__.update(**pre.__dict__)

        self.rv_core = np.linspace(0, self.R_c, num=r_res)

        if full_setup:
            # add derived parameters - that won't evolve
            # find central density to match input CMF
            self.L = np.sqrt((3 * self.K_0 * np.log(self.rho_cen / self.rho_0 + 1)) / (
                    2 * np.pi * p.G * self.rho_0 * self.rho_cen))  # initial guess
            self.rho_cen = self.get_central_density(m1=self.CMF * p.M_E, **kwargs)

            self.L = np.sqrt((3 * self.K_0 * np.log(self.rho_cen / self.rho_0 + 1)) / (
                    2 * np.pi * p.G * self.rho_0 * self.rho_cen))  # length scale
            self.D = np.sqrt(3 * self.C_p / (2 * np.pi * self.alpha_c * self.rho_cen * p.G))  # another length scale
            self.M_c = mass_profile(self.R_c, rho_cen=self.rho_cen, L=self.L)  # mass of core

            self.pv_core = pressure_profile(r=self.rv_core, R_c=self.R_c, p_cmb=self.p_cmb, rho_cen=self.rho_cen,
                                            L=self.L)

            # initial conditions
            self.evolve_profiles(**kwargs)


    def evolve_profiles(self, **kwargs):

        if np.size(self.T_cmb) > 1:
            print('using full soluton for T_cmb with shape', np.shape(self.T_cmb))
            self.T_cmb = np.expand_dims(self.T_cmb, axis=0)
            print(' new T_cmb shape', np.shape(self.T_cmb))
            # self.T_cen = np.tile(self.T_cmb,(len(self.rv_core),1)) / np.exp(-self.R_c ** 2 / self.D ** 2)
            # print('T_cen shape ', np.shape(self.T_cen))
        # else:
        self.T_cen = self.T_cmb / np.exp(-self.R_c ** 2 / self.D ** 2)

        # update profiles in core
        try:
            self.Tv_core = temperature_profile(self.rv_core, T_cen=self.T_cen, D=self.D)
        except ValueError:
            print('catching valueerror - updating profiles with full time evol')
            self.Tv_core = temperature_profile(r=self.rv_core, T_cen=self.T_cen.T, D=self.D)
            print('  Tv_core solution shape', np.shape(self.Tv_core))

        # update inner core boundary
        self.R_ic, self.T_ic = get_inner_core(**self.__dict__)

        self.M_oc = mass_profile(r=self.R_c, rmin=self.R_ic, **self.__dict__)  # mass of outer core
        self.rho_ic = density_profile(r=self.R_ic, **self.__dict__)
        self.g_ic = gravity_profile(r=self.R_ic, **self.__dict__)
        self.psi_ic = gravpotential_profile(r=self.R_ic, **self.__dict__)

        # update compositions
        if self.C_r is None:
            self.C_r = 1 / (self.gamma_melt - self.gamma_ad) * self.T_ic / self.T_c / (
                        self.rho_ic * self.g_ic)  # Gubbins 2003 eq 36 relating dR_i/dt to dT/dt
        self.chi = self.chi_0 / (1 - (self.R_ic / self.R_c) ** 3 + self.delta * (self.R_ic / self.R_c) ** 3)


    def get_central_density(self, m1=None, tol=0.0001, max_iter=1000, **kwargs):
        """ iterate rho_cen until mass profile produces the desired m1 value from CMF parameterisation. tol in kg"""
        rho_cen0 = self.rho_cen  # initial guess
        m0 = mass_profile(r=self.R_c, L=self.L, rho_cen=rho_cen0)  # current mass of core
        n = 0
        # print('initial error:', abs(m1 - m0)/m1)
        while (tol < abs(m1 - m0) / m1) and (n < max_iter):
            n = n + 1
            rho_cen1 = rho_cen0 * m1 / m0
            L = np.sqrt((3 * self.K_0 * np.log(rho_cen1 / self.rho_0 + 1)) / (
                    2 * np.pi * p.G * self.rho_0 * rho_cen1))  # length scale
            # D = np.sqrt(3 * c.C_p / (2 * np.pi * c.alpha_c * rho_cen1 * p.G))  # another length scale
            m0 = mass_profile(self.R_c, rho_cen=rho_cen1, L=L)  # mass of core
            # print('guess', n, ':', rho_cen1, 'kg/m3', 'rho1/rho0:', rho_cen1/rho_cen0, 'm1/m0', m1/m0)
            rho_cen0 = rho_cen1

        return rho_cen0


""" adapted from Nimmo+ 2004 G J Int"""

def K_0(L=7272e3, rho_0=7019, rho_cen=12500):
    # compressibility at 0 pressure -> let this be fixed at Nimmo value (todo: iterate?)
    return L**2 * 2 * np.pi * p.G * rho_0 * rho_cen / (3 * np.log(rho_cen/rho_0 + 1))


def adiabat(p, T_cmb=None, p_cmb=None, T_a1=None, T_a2=None, **kwargs):
    """ adiabatic temperature as a function of pressure determined by NImmo from fitting adiabatic T(r) profile and relating z to pressure
    therefore constants T_a1 etc may vary with different central densities rho_cen and T_cen ?? test
    T_c, p_c : temperature and pressure at CMB"""

    return T_cmb * (1 + T_a1 * p + T_a2 * p ** 2) / (1 + T_a1 * p_cmb + T_a2 * p_cmb ** 2)


def melt_curve(p, theta=None, T_m0=1695, T_m1=10.9e-12, T_m2=-8.0e-24, chi=None, which='nimmo', **kwargs):
    """ theta depends especially on concentration of S and Si (reduces T at ICB), but other data are for pure Fe
    uncertainty on Tm1 and Tm2 is about 1.1e-12 Pa^-1 and 1.5e-24 Pa^-2
    stixrude parameterisation depends directly on concentration of light elements (can assume constant if mostly S/Si contributing)"""
    if which == 'stixrude':
        # core melt temperature parameterisation from Stixrude (2014)
        return 6500 * (p*1e-9 / 340)**0.515 / (1 - np.log(1 - chi))
    elif which == 'nimmo':
        return T_m0 * (1 - theta) * (1 + T_m1 * p + T_m2 * p ** 2)


def get_inner_core(T_cmb=None, p_cmb=None, rv_core=None, pv_core=None, T_cen=None, **kwargs):
    # todo: fit to adiabat equation for different T, p at outer core radius

    Tm = melt_curve(pv_core, **kwargs)
    # solve adiabatic T(p) given T_c
    if np.size(T_cmb) == 1:
        Tp = adiabat(pv_core, T_cmb=T_cmb, p_cmb=p_cmb, **kwargs)
        idx = np.argwhere(np.diff(np.sign(Tp - Tm))).flatten()
        if np.size(Tp[idx]) > 0:
            T_ic = Tp[idx].item()
            R_ic = rv_core[idx].item()
        else:
            # no solid inner core
            T_ic = T_cen
            R_ic = 0

    else:
        print('\ngetting inner core w full solution, T_cmb', np.shape(T_cmb))
        # using full solution for T_cmb
        T_ic = np.zeros(np.shape(T_cmb.T))
        R_ic = np.zeros(np.shape(T_cmb.T))
        for ii, ti in enumerate(T_cmb.T):
            Tp = adiabat(pv_core, T_cmb=ti, p_cmb=p_cmb, **kwargs)
            idx = np.argwhere(np.diff(np.sign(Tp - Tm))).flatten()

            # try:
            if np.size(Tp[idx]) > 0:
                T_ic[ii] = Tp[idx].item()
                R_ic[ii] = rv_core[idx].item()
            else:
                # no solid inner core
                T_ic[ii] = T_cen.T[ii][0]
                R_ic[ii] = 0
            # except ValueError as e:
            #     plt.figure()
            #     plt.plot(Tp, pv_core, 'k-', label='adiabat')
            #     plt.plot(Tm, pv_core, 'k--', label='solidus')
            #     plt.legend()
            #     plt.title('timestep'+str(ii)+'$T_cmb$'+str(ti))
            #     plt.xlabel('T (K)')
            #     plt.ylabel('P (Pa)')
            #     plt.show()
            #     raise(e)

    return R_ic, T_ic

def temperature_profile(r, T_cen=None, D=None, **kwargs):
    return T_cen * np.exp(-r ** 2 / D ** 2)


def pressure_profile(r, R_c=None, p_cmb=None, rho_cen=None, L=None, **kwargs):
    """ p_c : pressure at CMB """
    # should this evolve?? no??
    # print('r', r, 'R_c', R_c, 'p_c', p_c, 'rho_cen', rho_cen, 'L', L)
    I = lambda x: (((3 * x ** 2 / 10) - (L ** 2 / 5)) * np.exp(-x ** 2 / L ** 2))
    return p_cmb + 4 * np.pi * p.G * rho_cen ** 2 / 3 * (I(R_c) - I(r))


def mass_profile(r, rmin=0, rho_cen=None, L=None, **kwargs):
    I = lambda x: (-(L ** 2 / 2) * x * np.exp(-(x ** 2) / L ** 2) + L ** 3 / 4 * np.sqrt(np.pi) * erf(x / L))
    m = 4 * np.pi * rho_cen * (I(r) - I(rmin))
    return m


def gravity_profile(r, rho_cen=None, L=None, **kwargs):
    g = 4 * np.pi / 3 * p.G * rho_cen * r * (1 - (3 * r ** 2 / (5 * L ** 2)))
    return g


def density_profile(r, rho_cen=None, L=None, **kwargs):
    rho = rho_cen * np.exp(-r ** 2 / L ** 2)
    return rho


def gravpotential_profile(r, R_c=None, rho_cen=None, L=None, **kwargs):
    # gravitational potential energy, relative to 0 potential at the cmb
    I = lambda x: 2 / 3 * np.pi * p.G * rho_cen * x ** 2 * (1 - (3 * x ** 2 / (10 * L ** 2)))
    psi = I(r) - I(R_c)
    return psi


def get_cmb_pressure(T_c, **kwargs):
    # depends on thermodynamic properties of mantle
    # todo
    return None


def grav_heat(L=None, R_c=None, R_ic=None, rho_cen=None, rho_ic=None, chi=None, delta_rho=None, psi_ic=None, C_r=None,
              M_oc=None, **kwargs):
    if np.size(R_ic) == 1 and R_ic == 0:
        return 0
    else:
        # gravitational heat - change in energy associated with the separation of the light component
        C_c = 4 * np.pi * R_ic ** 2 * rho_ic * chi / M_oc
        beta_c = 1 / chi * delta_rho / rho_ic  # compositional expansion coefficient - think rho here is rho_ic but maybe rho at cmb? NImmo eq 21

        # analytic nimmo way
        C_g2 = 3 * L ** 2 / 16 - R_c ** 2 / 2 * (1 - (3 * R_c ** 2 / (10 * L ** 2)))
        f_g = lambda x: ((3 / 20 * x ** 5 - L ** 2 / 8 * x ** 3 - L ** 2 * C_g2 * x) * np.exp(
            -x ** 2 / L ** 2) + C_g2 / 2 * L ** 3 * np.sqrt(np.pi) * erf(x / L))
        I_g = 8 * np.pi ** 2 * rho_cen ** 2 * p.G / 3 * (f_g(R_c) - f_g(R_ic))
        Q_g_bar = (I_g - M_oc * psi_ic) * beta_c * C_c * C_r

        if np.size(R_ic) > 1:
            Q_g_bar = np.squeeze(Q_g_bar)

        return Q_g_bar


def latent_heat(R_ic=None, L_h=None, rho_ic=None, C_r=None, **kwargs):
    # latent heat
    Q_L_bar = 4 * np.pi * R_ic ** 2 * L_h * rho_ic * C_r
    if np.size(R_ic) > 1:
        Q_L_bar = np.squeeze(Q_L_bar)
    return Q_L_bar


def specific_heat(L=None, D=None, T_cen=None, rho_cen=None, R_c=None, C_p=None, T_cmb=None, **kwargs):
    # specific heat
    A = np.sqrt((1 / L ** 2 + 1 / D ** 2) ** -1)
    I_s = 4 * np.pi * T_cen * rho_cen * (
            -A ** 2 * R_c / 2 * np.exp(-R_c ** 2 / A ** 2) + A ** 3 * np.sqrt(np.pi) / 4 * erf(R_c / A))
    Q_s_bar = -C_p / T_cmb * I_s  # Q_s_bar = Q_s / dT_c/dtrho_ic
    if np.size(T_cmb) > 1:
        Q_s_bar = np.squeeze(Q_s_bar)
    return Q_s_bar


def radioactive_heat(L=None, D=None, rho_cen=None, T_cen=None, R_c=None, M_c=None, t=None, #h_rad_c=None,
                     tau_i=np.array([1250, 4468, 703.8, 14050]),  # half life in Myr
                     h_i=np.array([28.761e-6, 94.946e-6, 568.402e-6, 26.368e-6]),  # heat production in W/kg
                        c_i=np.array([30.4e-9, 22.7e-9, 0.16e-9, 85e-9]),  # BSE concentration in kg/kg
                     part_i=np.array([0.1158, 0.0655, 0.0655, 0.0588]),
                    age=4.5, x_Eu=1,
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
    h_rad_c = th.h_rad(t, c_i*part_i, h_i, tau_i, age, x_Eu=x_Eu)

    Q_r = M_c * h_rad_c
    # print('   in rad heat: M_c', M_c, 'h_rad_c', h_rad_c )
    return Q_r


def core_gradient(c, test_cooling_rate=None, t=None, **kwargs):
    c.Q_r = radioactive_heat(t=t, **c.__dict__)
    c.Q_s_bar = specific_heat(**c.__dict__)
    c.Q_l_bar = latent_heat(**c.__dict__)
    c.Q_g_bar = grav_heat(**c.__dict__)
    if test_cooling_rate is not None:
        print('Q_s =', c.Q_s_bar * test_cooling_rate * 1e-12, 'TW')
        print('Q_l =', c.Q_l_bar * test_cooling_rate * 1e-12, 'TW')
        print('Q_g =', c.Q_g_bar * test_cooling_rate * 1e-12, 'TW')
        print('Q_r =', c.Q_r * 1e-12, 'TW')

    dTdt = (c.Q_r - c.Q_cmb) / -(c.Q_l_bar + c.Q_g_bar + c.Q_s_bar)  # note typo in NImmo paper where Q bar should be -
    c.Q_l = c.Q_l_bar * dTdt
    c.Q_g = c.Q_g_bar * dTdt
    c.Q_s = c.Q_s_bar * dTdt

    print('dTdt numerator', (c.Q_r - c.Q_cmb), '/ denomenator', c.Q_l_bar + c.Q_g_bar + c.Q_s_bar)
    return dTdt


def core_heat_flow(c,  **kwargs):
    """
    L_h: latent heating
    delta_rho: density difference across inner core boundary due to presence of light element in outer core (measured for earth??)
    gamma_melt: dT/dp of melting for core
    gamma_ad: adiabatic dT/dp fore core
    chi_0: initial concentration of light element in outer core (expelled from inner core during freezing)
    delta: mass partition coefficient of light element between solid and liquid - results should be insensitive
    """

    c.evolve_profiles(**kwargs)

    # calculate heat fluxes in governing equation, Nimmo+ 2004 eq 30, where Q_bar = Q / dT_c/dt
    c.dTdt = core_gradient(c, **kwargs)
    return c


def print_core(c):
    print('\n----parameters')
    print('R_c = ', c.R_c*1e-3, 'km')
    print('L = ', c.L*1e-3, 'km')
    print('D = ', c.D*1e-3, 'km')
    print('p_cmb =', c.p_cmb*1e-9, 'GPa')
    print('C_r = ', c.C_r)
    # print('p_cen =', c.p_cen*1e-9, 'GPa')
    print('rho_cen =', c.rho_cen, 'kg/m3')
    print('\n----timestep')
    print('T_cmb = ', c.T_cmb, 'K')
    print('R_ic = ', c.R_ic*1e-3, 'km')
    print('T_ic =', c.T_ic, 'K')
    print('chi = ', c.chi)


def plot_TP(p=None, T_cmb=None, p_cmb=None, T_ic=None, **kwargs):
    import matplotlib.pyplot as plt

    if p is None:
        p = np.linspace(136e9, 365e9)
    Tp = adiabat(p, T_cmb=T_cmb, p_cmb=p_cmb, **kwargs)
    Tm = melt_curve(p, **kwargs)

    plt.figure()
    plt.plot(p, Tp, c='k', ls='-', label='adiabat')
    plt.plot(p, Tm, c='k', ls='--', label='solidus')
    plt.axhline(T_cmb, lw=0.5, c='r', label='CMB temperature or pressure')
    plt.axvline(p_cmb, lw=0.5, c='r')
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
    c.rho_ = density_profile(r=r, **c.__dict__)
    c.rho_ic = density_profile(r=c.R_ic, **c.__dict__)


def plot_mass_profile(r=None, c=None, **kwargs):

    # update profiles

    import matplotlib.pyplot as plt

    plt.figure()
    if r is None:
        r = np.linspace(0, c.R_c)
    mass = mass_profile(r, **c.__dict__)
    plt.plot(mass/p.M_E, r*1e-3, c='k', lw=0.5, label='analytic profile')
    plt.axhline(c.R_c*1e-3, c='k', lw=0.5, ls='--')
    plt.axvline(mass[-1]/p.M_E, c='k', lw=1, ls='--', label='analytic core mass')
    plt.axvline(c.CMF*c.M_p/p.M_E, c='r', lw=0.5, label='required mass from CMF')
    plt.xlabel('mass (M_E)')
    plt.ylabel('radius (km)')
    print('central density:', c.rho_cen, 'kg/m3')
    print('mass error:', mass[-1] - c.CMF*c.M_p)
    plt.legend()
    plt.show()


def plot_psi_profile(r=None, c=None, **kwargs):

    # update profiles

    import matplotlib.pyplot as plt

    plt.figure()
    if r is None:
        r = np.linspace(0, c.R_c)
    psi = gravpotential_profile(r, **c.__dict__)
    psi_icb = gravpotential_profile(c.R_ic, **c.__dict__)
    print('psi_icb * M_oc', psi_icb * c.M_oc)
    plt.plot(psi, r*1e-3, c='k', lw=0.5, label='analytic profile')
    plt.axhline(c.R_c*1e-3, c='k', lw=0.5, ls='-')
    plt.axvline(c.R_c * 1e-3, c='k', lw=0.5, ls='-')
    plt.axvline(psi_icb, c='k', lw=1, ls='--', label='ICB')
    plt.axhline(c.R_ic*1e-3, c='k', lw=1, ls='--')
    plt.xlabel('gravitational potential (J)')
    plt.ylabel('radius (km)')
    plt.legend()
    # plt.show()


def plot_flux_dependence(cmf_vec=None, core_params=None, **kwargs):
    # for fixed cooling rate how do other fluxes change w cmf

    if core_params is None:
        # testing using Nimmo+ 2004 parameters - based on Table 4
        core_params = core_params_default
        core_params.update(test_params_default)

    fig, axes = plt.subplots(5, 1)

    if cmf_vec is None:
        cmf_vec = np.linspace(0.1, 0.5, num=8)

    for ii, cmf in enumerate(cmf_vec):
        print('\n\n---------------------------->', ii, '| cmf', cmf)
        core_params.update({'cmf_in': cmf})

        # calculate heat fluxes
        c = TerrestrialCore(**core_params)
        c = core_heat_flow(c, test_cooling_rate=c.test_cooling_rate)

        c.Q_r = radioactive_heat(**c.__dict__)
        c.Q_s_bar = specific_heat(**c.__dict__)
        c.Q_l_bar = latent_heat(**c.__dict__)
        c.Q_g_bar = grav_heat(**c.__dict__)
        total_flux = c.Q_r + (c.Q_s_bar + c.Q_l_bar + c.Q_g_bar)*c.test_cooling_rate

        axes[0].plot(cmf, c.Q_s_bar * c.test_cooling_rate * 1e-12, marker='.', c='k')
        axes[1].plot(cmf, c.Q_l_bar * c.test_cooling_rate * 1e-12, marker='.', c='k')
        axes[2].plot(cmf, c.Q_g_bar * c.test_cooling_rate * 1e-12, marker='.', c='k')
        axes[3].plot(cmf, c.Q_r * 1e-12, marker='.', c='k')
        axes[4].plot(cmf, total_flux * 1e-12, marker='.', c='k')

        print('R_ic', c.R_ic)

    plt.xlabel('CMF')
    # plt.ylabel('flux (TW)')
    axes[0].set_ylabel('specific heat (TW)')
    axes[1].set_ylabel('latent heat (TW)')
    axes[2].set_ylabel('grav. heat (TW)')
    axes[3].set_ylabel('rad. heat (TW)')
    axes[4].set_ylabel('total (TW)')
    axes[0].legend()
    plt.show()


def LHS(t, y, c=None, **kwargs):
    """ ODE equation to solve, LHS = 0 """

    c.T_cmb = y
    c.evolve_profiles(**kwargs)

    dTdt = core_gradient(c, t=t, test_cooling_rate=None)
    return dTdt


def solve_core(c, t0=0, tf=4.5, plot=True, **kwargs):
    T_c0 = c.T_cmb  # initialise
    f = integrate.solve_ivp(fun=lambda t, y: LHS(t, y, **dict(c=c, **kwargs)),
                            t_span=(t0 * 1e9 * p.years2sec, tf * 1e9 * p.years2sec), y0=[T_c0],
                            max_step=100e6 * p.years2sec,
                            method='RK45',
                        )

    print('>>>>finished integration >>>', 'T_cmb solution', np.shape(c.T_cmb))
    # update profiles with T(t) solution
    c.T_cmb = f.y[0]
    c.evolve_profiles()
    c.dTdt = core_gradient(c, t=f.t)
    c.Q_r = radioactive_heat(t=f.t, **c.__dict__)
    c.Q_s = specific_heat(**c.__dict__) * c.dTdt
    c.Q_l = latent_heat(**c.__dict__) * c.dTdt
    c.Q_g = grav_heat(**c.__dict__) * c.dTdt

    c.T_cmb = np.squeeze(c.T_cmb)  # fix
    c.T_cen = np.squeeze(c.T_cen)
    if plot:
        fig, axes = plt.subplots(3, 1)
        axes[0].plot(f.t*p.sec2Gyr, c.T_cmb, 'k-')
        axes[0].set_ylabel('$T_{cmb}$ (K)')
        axes[1].plot(f.t*p.sec2Gyr, c.Q_r * 1e-12, label='Q_r')
        axes[1].plot(f.t * p.sec2Gyr, c.Q_s * 1e-12, label='Q_s')
        axes[1].plot(f.t * p.sec2Gyr, c.Q_l * 1e-12, label='Q_l')
        axes[1].plot(f.t * p.sec2Gyr, c.Q_g * 1e-12, label='Q_g')
        axes[1].set_ylabel('$Q$ (TW)')
        axes[1].legend()
        axes[2].plot(f.t * p.sec2Gyr, c.dTdt / p.sec2Gyr)
        axes[2].set_ylabel('dT/dt (K/Gyr)')
        axes[-1].set_xlabel('t (Gyr)')
    return c


T_cmb0 = 4800
c = TerrestrialCore(cmf_in=0.3, T_cmb=T_cmb0, **core_params_default, **test_params_default)
solve_core(c)

# plot_psi_profile(r=None, c=c)
# plot_psi_profile(r=None, c=TerrestrialCore(cmf_in=0.35, **core_params_default, **test_params_default))

plt.show()
# plot_flux_dependence()

# L_test = np.sqrt((3 * core_params_default['K_0'] * np.log(core_params_default['rho_cen'] / core_params_default['rho_0'] + 1)) / (
#                     2 * np.pi * p.G * core_params_default['rho_0'] * core_params_default['rho_cen']))
# M_c = mass_profile(3480e3, rmin=0, rho_cen=12500, L=L_test)
# print('predicted cmf:', M_c/p.M_E)