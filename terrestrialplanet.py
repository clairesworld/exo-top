import numpy as np
from scipy import interpolate
import thermal as therm
import topography as top
import astroenv as ast
import structure as st 



######      CLASS DEF     ######
class TerrestrialPlanet():
    
    ###### PHYSICAL CONSTANTS ######
    global M_E, R_E, L_sun, years2sec, sec2Gyr, AU2m, R_b, sb
    M_E = 5.972e24 # earth mass in kg
    R_E = 6371e3 # earth radius in m
    L_sun =  3.9e26 # solar luminosity in W
    years2sec = 31557600
    sec2Gyr = 1e-9/years2sec
    AU2m = 1.5e11
    R_b = 8.3144598 # universal gas constant in J mol −1 K −1
    sb = 5.67e-8 # Stefan Boltzmann constant in W m^-2 K^-4

    def __init__(self, **kwargs):
        """ Initialize with kwargs, TODO: defaults? """
        # define default attributes
        # [238U, 235U, 232Th, 40K] order for radioisotopes
        print(M_E)
        default_attr = dict(lambda_n = np.array([0.15541417, 
                                                 0.98458406, 0.04951051, 0.55011681])*1e-9/years2sec, # decay constant
                            p_n = [95.13e-6, 568.47e-6, 26.3e-6, 28.47e-6], # heating rate in W kg^-1 from Treatise, Dye (2012) 
                            K_0 = 0.0117e-2, # ratio of 40-K to total K at time 0 (in Treatise on Geophysics)
                            U_0_235 = 0.0072, # ratio of 235-U to total U at time 0 (in Treatise on Geophysics)
                            U_0_238 = 0.9927, # ratio of 238-U to total U at time 0 (in Treatise on Geophysics)
                            Th_0 = 1, # ratio of 232-Th to total Th at time 0 (in Treatise on Geophysics)
                            R_p0 = None, # fixed R_p input value
                            R_c0 = None, # fixed R_c input value
                            T_s = None,
                            CMF = 0,
                            complexity = 3,
                            # TODO: give everything a default?
                            ident = '%.2f'%(kwargs['M_p']/M_E)+' M$_E$, CMF='+'%.1f'%(kwargs['CMF']) # run id
                           )  
        default_attr.update(kwargs) # add input parameters, use default if not given
        self.__dict__.update((k,v) for k,v in default_attr.items())
        
        self.init_derived() # add derived parameters
#         for key, value in kwargs.items():
#             setattr(self, key, value)

            
    def init_derived(self, **kwargs):
        """ Parameters derived from input parameters """
        if self.R_p0 is None:
            self.R_p = st.radius_zeng(self.M_p, self.CMF)*R_E # in m
        else:
            self.R_p = self.R_p0
        if self.R_c0 is None:
            self.M_m = self.M_p*(1 - self.CMF) # mass of mantle
            self.CRF = self.CMF**0.5 # Zeng & Jacobsen 2017
            self.M_c = self.M_p*self.CMF # M_p - M_m
            self.R_c = self.R_p*self.CRF
            #R_c = radius_seager(M_p*CMF, CMF=0, m1=4.34*M_E, r1=2.23*R_E) # EoS for iron... is this consistent?
        else:
            self.R_c = self.R_c0
            self.CRF = self.R_c0/self.R_p0
            self.M_c = 4/3*np.pi*self.R_c**3 * self.rho_c
            self.CMF = self.M_c/self.M_p
            self.M_m = 4/3*np.pi*(self.R_p**3 - self.R_c**3)*self.rho_m  #M_p - M_c
        self.SA_p = st.SA(R=self.R_p)
        self.SA_c = st.SA(R=self.R_c) # core surface area 
        self.g_sfc = st.grav(self.M_p, self.R_p)
        if self.CMF>0:
            self.g_cmb = st.grav(self.M_c, self.R_c)
        else:
            self.g_cmb = 0
        self.kappa_m = therm.thermal_diffusivity(self.k_m, self.rho_m, self.c_m)

        if self.T_s is None:
            self.q_out = ast.q_sfc_outgoing(R_p=self.R_p, SA_p=self.SA_p, **kwargs)
            self.T_s = ast.T_sfc(self.q_out)

        # radiogenic element abundance rel. to U
        self.c_n = np.array([self.U_0_238, self.U_0_235, self.Th_0, self.K_0])*np.array([1, 1, self.X_Th/self.X_U, 
                                                                                         self.X_K/self.X_U])

