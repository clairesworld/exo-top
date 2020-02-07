import numpy as np
from scipy import interpolate
import thermal as therm
import topography as top
import astroenvironment as ast
import structure as struct 
import parameters as p

# todo: update attributes function

######      CLASS DEF     ######
class TerrestrialPlanet():

    def __init__(self, **kwargs):
        """ Initialize with kwargs, TODO: defaults? """
        # define default attributes
        default_attr = dict(lambda_n = p.lambda_n, 
                            p_n = p.p_n,
                            K_0 = p.K_0,
                            U_0_235 = p.U_0_235,
                            U_0_238 = p.U_0_238, 
                            Th_0 = p.Th_0, 
                            R_p0 = None, # fixed R_p input value
                            R_c0 = None, # fixed R_c input value
                            T_s = None,
                            CMF = 0, # needs a value to avoid errors
                            L = 1, # stellar luminosity in solar units
                            Alb = 0, # planetary albedo
                            sma = 1, # semimajor axis in au
                            Ra_crit_u = 660, # critical Rayleigh number (in Driscoll & Bercovici 2014)
                            rho_c = 8000, # Density of iron core in kg m^-3 
                            rho_m = 3300, # Density of silicate mantle in kg m^-3 rho_lith = 2800,
                            # what pressure should you take these densities at?
                            c_m = 1200, # specific heat capacity from Dorn, Noack & Rozal 2018 in J kg−1 K−1 
                            c_c = 530, # speific heat capacity iron core Nimmo+ 1997
                            alpha_m = 2e-5, # thermal expansivity of silicate mantle in K^-1
                            k_m = 4, # thermal conductivity of silicate mantle in W m^−1 K^−1
                            k_lm = 10, # thermal conductivity lower mantle in W m^−1 K^−1 from Driscoll & Bercovici
                            X_K = 250, # initial abundance of K in wt ppm in Treatise on Geophysics
                            X_U = 2e-2, # initial abundane of U in wt ppm ""
                            X_Th = 7e-2, # initial abundance of Th in wt ppm ""
                           )  
        
        # add input parameters, use default if not given
        default_attr.update(kwargs) 
        self.__dict__.update((k,v) for k,v in default_attr.items())
        # add derived parameters
        self.init_derived(**kwargs) 

            
    def init_derived(self, **kwargs):
        """ Parameters derived from input parameters """
        try:
            self.ident = '%.2f'%(kwargs['M_p']/p.M_E)+' M$_E$, CMF='+'%.1f'%(kwargs['CMF']) # run id
        except KeyError:
            self.ident = '%.2f'%(kwargs['M_p']/p.M_E)+' M$_E$' # run id
        if self.R_p0 is None:
            self.R_p = struct.radius_zeng(self.M_p, self.CMF)*p.R_E # in m
        else:
            self.R_p = self.R_p0
        if self.R_c0 is None:
            self.M_m = self.M_p*(1 - self.CMF) # mass of mantle
            self.CRF = self.CMF**0.5 # Zeng & Jacobsen 2017
            self.M_c = self.M_p*self.CMF # M_p - M_m
            self.R_c = self.R_p*self.CRF
            #R_c = struct.radius_seager(M_p*CMF, CMF=0, m1=4.34*M_E, r1=2.23*R_E) # EoS for iron... is this consistent?
        else:
            self.R_c = self.R_c0
            self.CRF = self.R_c0/self.R_p0
            self.M_c = 4/3*np.pi*self.R_c**3 * self.rho_c
            self.CMF = self.M_c/self.M_p
            self.M_m = 4/3*np.pi*(self.R_p**3 - self.R_c**3)*self.rho_m  #M_p - M_c
        self.SA_p = struct.SA(R=self.R_p)
        self.SA_c = struct.SA(R=self.R_c) # core surface area 
        self.g_sfc = struct.grav(self.M_p, self.R_p)
        if self.CMF>0:
            self.g_cmb = struct.grav(self.M_c, self.R_c)
        else:
            self.g_cmb = 0
        self.kappa_m = therm.thermal_diffusivity(self.k_m, self.rho_m, self.c_m)

        if self.T_s is None:
            self.q_out = ast.q_sfc_outgoing(R_p=self.R_p, SA_p=self.SA_p, **kwargs)
            self.T_s = ast.T_sfc(self.q_out)

        # radiogenic element abundance rel. to U
        self.c_n = [a*b for a,b in zip([self.U_0_238, self.U_0_235, self.Th_0, self.K_0],[1, 1, self.X_Th/self.X_U, 
                                                                                         self.X_K/self.X_U])]
#         self.c_n = np.array([self.U_0_238, self.U_0_235, self.Th_0, self.K_0])*np.array([1, 1, self.X_Th/self.X_U, 
#                                                                                          self.X_K/self.X_U])


    def set_attrs(self, **kwargs):
        # todo??
        pass