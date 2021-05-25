import numpy as np
from . import parameters
from . import astroenvironment as ast
from . import geometry as geom
from . import thermal

# todo: update attributes function

######      CLASS DEF     ######
class TerrestrialPlanet():

    def __init__(self, **kwargs):
        """ Initialize with kwargs, TODO: defaults? """
        # define default attributes
        default_attr = dict(
                            # bulk property defaults
                            M_p = parameters.M_E,
                            sma = 1, # semimajor axis in au
                            L = 1, # stellar luminosity in solar units
                            Alb = 0, # planetary albedo 
                            Th_0 = parameters.Th_0, 
                            R_p0 = None, # fixed R_p input value, if none then comes from mass
                            R_c0 = None, # fixed R_c input value " "
                            T_s = None, # if none then calculated from energy balance
                            CMF = 0.3, # needs a value to avoid errors
                            rho_c = 7200, # Density of iron core in kg m^-3 
                            rho_m = 3500, # Density of silicate mantle in kg m^-3 rho_lith = 2800,
                            Ra_crit_u = 450, # critical Rayleigh number (660 in Driscoll & Bercovici 2014)
                            beta_u = 1/3, # defaults to 1/3
                            beta_c = 1/3, # defaults to 1/3
                            # what pressure should you take densities at?
                            
                            # thermodynamic defaults
                            c_m = 1200, # specific heat capacity from Dorn, Noack & Rozal 2018 in J kg−1 K−1 
                            c_c = 840, # speific heat capacity iron core Thiriet+ 2019
                            alpha_m = 2e-5, # thermal expansivity of silicate mantle in K^-1
                            k_m = 4, # thermal conductivity of silicate mantle in W m^−1 K^−1
                            k_lm = 4, # thermal conductivity lower mantle in W m^−1 K^−1, 10 from Driscoll & Bercovici
                            
                            # radioisotope defaults
                            lambda_n = parameters.lambda_n, # astrophysically constant
                            p_n = parameters.p_n,
                            K_0 = parameters.K_0,
                            U_0_235 = parameters.U_0_235,
                            U_0_238 = parameters.U_0_238,
                            X_K = 250, # initial abundance of K in wt ppm in Treatise on Geophysics
                            X_U = 2e-2, # initial abundane of U in wt ppm ""
                            X_Th = 7e-2, # initial abundance of Th in wt ppm ""
                            H_f = 4.6e-12, # radiogenic heating in W/kg at 4.5 Gyr from Javoy (1999), CI chondrites
                            H_0 = 22.65771894e-12,  # radiogenic heating at t0 based on above
                            
                            # viscosity defaults
                            a_rh=2.44, # for beta=1/3 from Thiriet+ (2019)
                            eta_pre=None,
                            eta_0 = 1e21, # reference eta from Thiriet+ (2019)
                            T_ref = 1600, # reference T from Thiriet+ (2019)
                            Ea=300e3, # activation energy in J, K&W (1993) dry olivine
                            V_rh=6e-6, # activation volume in m^3, K&W (1993)  dry olivine
                            mu=80e9, # shear modulus in Pa, K&W (1993)  dry olivine
                            A_rh=8.7e15, # pre-exponential factor in s^-1, K&W (1993)  dry olivine
                            h_rh=2.07e-3, # grain size in m, K&W (1993)  dry olivine
                            B_rh=0.5e-9, # Burgers vector, K&W (1993)  dry olivine
                            m_rh=2.5, # grain size exponent, K&W (1993)  dry olivine
                           )  
        
        
        # add input parameters, use default if not given
        default_attr.update(kwargs) 
        self.__dict__.update((k,v) for k,v in default_attr.items())
        # add derived parameters
        self.init_derived(**kwargs)
        # placeholders for thermal evol
        self.T_m = None
        self.T_c = None
        self.D_l = None
        self.t = None

            
    def init_derived(self, **kwargs):
        """ Parameters derived from input parameters """
        if 'ident' not in kwargs.keys():
            try:
                self.ident = '%.2f'%(self.M_p/parameters.M_E)+' M$_E$, CMF='+'%.1f'%(self.CMF) # run id
            except KeyError:
                self.ident = '%.2f'%(self.M_p/parameters.M_E)+' M$_E$' # run id
        if self.R_p0 is None:
            self.R_p = ast.radius_zeng(self.M_p, self.CMF)*parameters.R_E # in m
        else:
            self.R_p = self.R_p0

        if self.CMF is not None:
            self.CRF = self.CMF ** 2
            self.R_c = self.R_p * self.CRF
            self.M_m = self.M_p * (
                    1 - self.CMF)  # mass of mantle, updated immediately in thermal code including lid dynamics
            self.M_c = self.M_p * self.CMF  # M_p - M_m
        elif self.R_c0 is not None:
            self.R_c = self.R_c0
            self.CRF = self.R_c0 / self.R_p0
            self.M_c = 4 / 3 * np.pi * self.R_c ** 3 * self.rho_c
            self.CMF = self.M_c / self.M_p
            self.M_m = 4 / 3 * np.pi * (self.R_p ** 3 - self.R_c ** 3) * self.rho_m  # M_p - M_c
        else:
            self.M_m = self.M_p * (
                        1 - self.CMF)  # mass of mantle, updated immediately in thermal code including lid dynamics
            self.CRF = self.CMF ** 0.5  # Zeng & Jacobsen 2017
            self.M_c = self.M_p * self.CMF  # M_p - M_m
            self.R_c = self.R_p * self.CRF

        self.SA_p = geom.SA(R=self.R_p)
        self.SA_c = geom.SA(R=self.R_c) # core surface area 
        self.g_sfc = ast.grav(self.M_p, self.R_p)
        if self.CMF>0:
            self.g_cmb = ast.grav(self.M_c, self.R_c)
        else:
            self.g_cmb = 0
        self.kappa_m = thermal.thermal_diffusivity(self.k_m, self.rho_m, self.c_m)

        if self.T_s is None:
            self.q_out = ast.q_sfc_outgoing(R_p=self.R_p, SA_p=self.SA_p, L=self.L, Alb=self.Alb, sma=self.sma)
            self.T_s = ast.T_sfc(self.q_out)

        # radiogenic element abundance rel. to U
        self.c_n = [a*b for a,b in zip([self.U_0_238, self.U_0_235, self.Th_0, self.K_0],[1, 1, self.X_Th/self.X_U, 
                                                                                         self.X_K/self.X_U])]
#         self.c_n = np.array([self.U_0_238, self.U_0_235, self.Th_0, self.K_0])*np.array([1, 1, self.X_Th/self.X_U, 
#                                                                                          self.X_K/self.X_U])


    def set_attrs(self, **kwargs):
        self.__dict__.update((k,v) for k,v in kwargs.items())
    
    def copy(self, **kwargs):
        self_args = self.__dict__
        self_args.update(kwargs)
        new = TerrestrialPlanet(self_args)
        return new

    def nondimensionalise(self, **kwargs):
        # set natural scales
        T = np.maximum(self.T_c, self.T_m) - self.T_s  # fucked for cases where T_c < T_m
        T = self.T_c - self.T_s
        L = self.R_p - self.R_c
        try:
            T0 = T[0]
        except IndexError:
            T0 = T  # steady state?

        self.T_m_prime = self.T_m/T0
        self.T_c_prime = self.T_c/T0
        self.dT_m_prime = self.dT_m/T0
        self.dT_rh_prime = self.dT_rh/T0
        self.T_l_prime = self.T_l/T0

        self.d_m_prime = self.d_m/L
        self.D_l_prime = self.D_l/L
        self.delta_rh_prime = self.delta_rh/L

        self.alpha_m_prime = 1
        self.rho_m_prime = 1
        self.g_sfc_prime = 1
        self.k_m_prime = 1
        self.kappa_m_prime = 1
        self.c_m_prime = 1
        self.c_c_prime = 1

        self.eta_m_prime = 1/self.Ra_i * np.exp(-np.log(self.delta_eta) * self.T_m_prime)
        self.dyn_top_rms_prime = self.dyn_top_rms/(self.alpha_m*T0*L)
        self.heuristic_h_prime = self.heuristic_h/(self.alpha_m*T0*L)

        self.T_scale = T
        self.L_scale = L
        # incomplete