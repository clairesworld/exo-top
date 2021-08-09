from model_1D import parameters as p

baseline_in = dict(
    ident='baseline',  # must match dict name ****_in
    M_p=5.972e24,
    sma=1,
    Alb=0,
    CMF=0.3,
    Ra_crit_u=450,
    rho_m=3500,
    rho_c=7200,
    c_m=1142,
    c_c=840,
    beta_u=1 / 3,
    k_m=4,
    alpha_m=2.5e-5,
    H_0=22.65771894e-12,  # not used -radiogenic heating in W/kg at t0  based on below
    H_f=4.6e-12,  # not usde - radiogenic heating in W/kg at 4.5 Gyr from Javoy (1999) BSE estimate from CI chondrites
    T_s=273,
    eta_pre=1.6e11,  # arrhenius pre-exponential factor - dry olivine
    Ea=300e3,
    x_Eu=1,  # at 0.7 this is similar to older Treatise values

    #     # viscosity
    #     nu_0 = 0.27e17, # use for constant visc option
    #     a_rh=2.44, # for beta=1/3 from Thiriet+ (2019)
    eta_0=1e21,  # reference eta from Thiriet+ (2019)
    T_ref=1600,  # reference T from Thiriet+ (2019)
    #     Ea=300e3, # activation energy in J, K&W (1993) dry olivine
    #     V_rh=6e-6, # activation volume in m^3, K&W (1993)  dry olivine
    #     mu=80e9, # shear modulus in Pa, K&W (1993)  dry olivine
    #     A_rh=8.7e15, # pre-exponential factor in s^-1, K&W (1993)  dry olivine
    #     h_rh=2.07e-3, # grain size in m, K&W (1993)  dry olivine
    #     B_rh=0.5e-9, # Burgers vector, K&W (1993)  dry olivine
    #     m_rh=2.5, # grain size exponent, K&W (1993)  dry olivine
)
baseline_run = dict(T_m0=1750, T_c0=2250, D_l0=150e3, tf=4.5, visc_type='KW')  # model params

Earthbaseline_in = dict(
    ident='Earthbaseline',  # must match dict name ****_in
    M_p=5.972e24,
    sma=1,
    Alb=0,
    CMF=0.3,
    Ra_crit_u=450,
    rho_m=3300,
    rho_c=7200,
    c_m=1200,
    c_c=7800,
    beta_u=1 / 3,
    k_m=3,
    alpha_m=3e-5,
    H_0=22.65771894e-12,  # radiogenic heating in W/kg at t0  based on below
    H_f=4.6e-12,  # radiogenic heating in W/kg at 4.5 Gyr from Javoy (1999) BSE estimate from CI chondrites

    #     # viscosity
    #     nu_0 = 0.27e17, # use for constant visc option
    #     a_rh=2.44, # for beta=1/3 from Thiriet+ (2019)
    #     eta_0 = 1e21, # reference eta from Thiriet+ (2019)
    #     T_ref = 1600, # reference T from Thiriet+ (2019)
    #     Ea=300e3, # activation energy in J, K&W (1993) dry olivine
    #     V_rh=6e-6, # activation volume in m^3, K&W (1993)  dry olivine
    #     mu=80e9, # shear modulus in Pa, K&W (1993)  dry olivine
    #     A_rh=8.7e15, # pre-exponential factor in s^-1, K&W (1993)  dry olivine
    #     h_rh=2.07e-3, # grain size in m, K&W (1993)  dry olivine
    #     B_rh=0.5e-9, # Burgers vector, K&W (1993)  dry olivine
    #     m_rh=2.5, # grain size exponent, K&W (1993)  dry olivine
)
Earthbaseline_run = dict(T_m0=1750, T_c0=2250, D_l0=150e3, tf=4.5, visc_type='KW',
                         backwards_cooling=False)  # model params

Earth_Tosi_in = dict(
    ident='Earth_Tosi',
    Alb=None,
    H_f=23e-12,  # final radiogenic heating in W/kg
    L=None,
    Ra_crit_u=450,
    R_p0=6370e3,
    R_c0=3840e3,
    alpha_m=2e-5,  # thermal expansivity
    k_m=4,  # silicate thermal conductivity
    rho_c=7200,  # core density
    rho_m=3500,  # mantle density
    rho_lith=None,
    c_m=1100,  # specific heat mantle
    c_c=800,  # specific heat for core in J/K/kg
    k_lm=4,  # thermal conduvtivity lower mantle
    beta_u=1 / 3,  # defaults to 1/3
    beta_c=1 / 3,  # defaults to 1/3
    a_rh=2.9,  # for beta=1/3
    Ea=3.35e5 + 2.5e9 * 4e-6,  # activation energy in J for viscosity law
    eta_pre=6.127e10 / 500,  # reference dynamic viscosity in Pa s
    T_ref=1600,  # viscosity law reference temperature in K
    T_s=293,  # fixed surface temp in K
    M_p=p.M_E,  # only used for gravity in this case
    sma=None,
)
Earth_Tosi_run = dict(T_m0=1700, T_c0=1900, D_l0=100e3, tf=4.5, visc_type='KW', rad_type='H0', backwards_cooling=True)

Earth_Foley_in = dict(
    f_melt0=None,
    ident='Earth_Foley',  # must match dict name ****_in
    cc_slope=1.5e-10,  # Pa-1
    rho_lith=3300,
    rho_melt=2800,
    L_melt=600e3,  # latnt heat of melting
    gamma_melt=3.1e-8,  # K Pa^-1
    gamma_m=1e-8,
    M_p=5.972e24,
    CMF=0.3,
    rho_m=4000,
    rho_c=7200,
    c_m=1250,
    c_c=840,
    beta_u=1 / 3,
    k_m=5,
    alpha_m=3e-5,
    H_0=22.65771894e-12,
    H_f=4.6e-12,  # not usde - radiogenic heating in W/kg at 4.5 Gyr from Javoy (1999) BSE estimate from CI chondrites
    T_s=273,
    Ea=300e3,
    a_rh=2.5,  # for beta=1/3
    eta_pre=4e10,
)
Earth_Foley_run = dict(T_m0=2000, T_c0=2000, D_l0=100e3, tf=6, visc_type='KW', rad_type='H0',
                       q_type='theta')  # model params

Venusbaseline_in = dict(
    ident='Venusbaseline',  # must match dict name ****_in
    M_p=4.867e24,
    CMF=0.3,
    Ra_crit_u=450,
    rho_m=3300,
    rho_c=7200,
    c_m=1200,
    c_c=7800,
    beta_u=1 / 3,
    a_rh=2.44,
    k_m=3,
    alpha_m=3e-5,
    T_s=730,
    Ea=300e3,
    # nu_0=0.27e17,  # use for constant visc case
    # H_0=22.65771894e-12,  # radiogenic heating in W/kg at t0 (Mars value from Thiriet - UPDATE)
    # H_f=4.6e-12  # radiogenic heating in W/kg at 4.5 Gyr from Javoy (1999) BSE estimate from CI chondrites
)
Venusbaseline_run = dict(T_m0=1750, T_c0=2250, D_l0=150e3, tf=4.5, visc_type='KW',
                         backwards_cooling=False)  # model params

VenusHuang15_in = dict(
    ident='VenusHuang15',  # must match dict name ****_in
    M_p=4.867e24,
    R_p0=6050e3,
    R_c0=3330e3,
    Ra_crit_u=450,
    rho_m=3300,
    rho_c=7200,
    c_m=1200,
    beta_u=1 / 3,
    a_rh=2.44,
    k_m=3.206,  # from thermal expansivity
    alpha_m=2e-5,
    T_s=730,
    Ea=300e3,
    # use 150 kJ/mol but with r-dependence... maybe assume 300 and set prefactor so you get same Ra for same <T>
    H_0=22.65771894e-12,  # radiogenic heating in W/kg at t0 (Mars value from Thiriet - UPDATE)
    H_f=4.6e-12  # radiogenic heating in W/kg at 4.5 Gyr from Javoy (1999) BSE estimate from CI chondrites
)
# include olivie-spinel phase change, spinel-perovskite phase change
VenusHuang15_run = dict(T_m0=1750, T_c0=2250, D_l0=200e3, tf=4.5, visc_type='KW',
                        backwards_cooling=False)  # model params

# VenusGolle_in = dict(
#     ident = 'VenusGolle', # must match dict name ****_in
#     M_p = 4.867e24,
#     R_p0 = 6050e3,
#     R_c0 = 3085e3,
#     D_l_const = 300,
#     Ra_crit_u = 450,
#     rho_m = 3300,
#     rho_c = 7200,
#     c_m = ,
#     beta_u = 1/3,
#     a_rh = 2.44 ,
#     k_m = ,
#     alpha_m =,
#     T_s = 730,
#     Ea = 300e3,
#     H_0 = 0,
#     Ra_F_const = 3e8,
#     nu_0 = 1e21/3300,
# )


VenusNimmo_V20_in = dict(
    ident='VenusNimmo_V20',  # must match dict name ****_in
    M_p=4.867e24,
    R_p0=6050e3,
    d_m_const=705e3,
    D_l_const=176e3,
    q_core_const=15e-3,
    Ra_crit_u=450,
    rho_m=3300,
    rho_c=7200,
    c_m=1200,  # from k, kappa
    beta_u=1 / 3,
    a_rh=2.44,
    k_m=3.17,
    alpha_m=3.28e-5,
    T_s=730,
    H_0=0,
    Ra_F_const=1.58e7,
    nu_0=0.27e17,
)
VenusNimmo_V20_run = dict(T_m0=1364 + 273, T_c0=2250, D_l0=176e3, tf=4.5, visc_type='constant',
                          backwards_cooling=False)  # model params

VenusNimmo_V10_in = dict(
    ident='VenusNimmo_V10',  # must match dict name ****_in
    M_p=4.867e24,
    R_p0=6050e3,
    d_m_const=705e3,
    D_l_const=176e3,
    q_core_const=15e-3,
    Ra_crit_u=450,
    rho_m=3300,
    rho_c=7200,
    c_m=1200,  # from k, kappa
    beta_u=1 / 3,
    a_rh=2.44,
    k_m=3.17,
    alpha_m=3.28e-5,
    T_s=730,
    H_0=0,
    Ra_F_const=7.86e6,
    nu_0=0.54e17,
)
VenusNimmo_V10_run = dict(T_m0=1367 + 273, T_c0=2250, D_l0=176e3, tf=4.5, visc_type='constant',
                          backwards_cooling=False)  # model params

VenusNimmo_V1_in = dict(
    ident='VenusNimmo_V10',  # must match dict name ****_in
    M_p=4.867e24,
    R_p0=6050e3,
    d_m_const=705e3,
    D_l_const=176e3,
    q_core_const=15e-3,
    Ra_crit_u=450,
    rho_m=3300,
    rho_c=7200,
    c_m=1200,  # from k, kappa
    beta_u=1 / 3,
    a_rh=2.44,
    k_m=3.17,
    alpha_m=3.28e-5,
    T_s=730,
    H_0=0,
    Ra_F_const=3.95e6,
    nu_0=0.29e18,
)
VenusNimmo_V1_run = dict(T_m0=1210 + 273, T_c0=2250, D_l0=176e3, tf=4.5, visc_type='constant',
                         backwards_cooling=False)  # model params

VenusKH92_in = dict(
    ident='VenusKH92',  # must match dict name ****_in
    M_p=4.867e24,
    R_p0=6050,
    D_l_const=130e3,
    d_m_const=2800e3,
    dT_const=1000,
    Ra_crit_u=450,
    rho_m=3300,
    c_m=909,
    beta_u=1 / 3,
    a_rh=2.44,
    k_m=3,
    alpha_m=3e-5,
    T_s=730,
    Ra_const=1e6,
    nu_0=1.93e22,  # calculate from Ra and dT
    H_0=0,
)

Mars_Breuer_in = dict(
    ident='Mars_Breuer',
    Alb=None,
    H_0=4e-12,  # final radiogenic heating in W/kg
    X_K=305,  # initial abundance of K in wt ppm
    X_U=16e-3,  # initial abundane of U in wt ppm
    X_Th=56e-3,  # initial abundance of Th in wt ppm
    L=None,
    Ra_crit_u=450,
    R_p0=3400e3,
    R_c0=1700e3,
    alpha_m=2e-5,  # thermal expansivity
    k_m=4,  # silicate thermal conductivity
    rho_c=7200,  # core density
    rho_m=3500,  # mantle density
    rho_lith=None,
    c_m=1142,  # <----??? TODO: check if you need constant volume c_p
    c_c=840,  # specific heat for core in J/K/kg
    k_lm=4,  # thermal conduvtivity lower mantle
    beta_u=1 / 3,  # defaults to 1/3
    beta_c=1 / 3,  # defaults to 1/3
    a_rh=2.21,  # for beta=1/3
    Ea=300e3,  # activation energy in J for viscosity law
    eta_0=1e21,  # reference dynamic viscosity in Pa s
    T_ref=1600,  # viscosity law reference temperature in K
    T_s=220,  # fixed surface temp in K
    M_p=6.39e23,  # only used for gravity in this case
    sma=None,
)
Mars_Breuer_run = dict(T_m0=1900, T_c0=2200, D_l0=100e3, tf=4.5, visc_type='KW', backwards_cooling=False)

Mars_Hauck_in = dict(
    ident='Mars_Hauck',
    Alb=None,
    H_0=5e-11,  # final radiogenic heating in W/kg
    X_K=305,  # initial abundance of K in wt ppm
    X_U=16e-3,  # initial abundane of U in wt ppm
    X_Th=56e-3,  # initial abundance of Th in wt ppm
    L=None,
    Ra_crit_u=450,
    R_p0=3390e3,
    R_c0=1550e3,
    alpha_m=2e-5,  # thermal expansivity
    k_m=4,  # silicate thermal conductivity
    rho_c=7200,  # core density
    rho_m=3527,  # mantle density
    rho_lith=None,
    c_m=1149,  # <----??? TODO: check if you need constant volume c_p
    c_c=571,  # specific heat for core in J/K/kg
    k_lm=4,  # thermal conduvtivity lower mantle
    beta_u=1 / 3,  # defaults to 1/3
    beta_c=1 / 3,  # defaults to 1/3
    A_rh=5.3e15,
    B_rh=5e-10,
    mu=80e9,
    h_rh=1e-3,
    m_rh=1,
    a_rh=2.21,  # for beta=1/3
    Ea=240e3,  # activation energy in J for viscosity law
    T_s=220,  # fixed surface temp in K
    M_p=6.39e23,  # only used for gravity in this case
    sma=None,
)
Mars_Hauck_run = dict(T_m0=1723, T_c0=2000, D_l0=100e3, tf=4.5, visc_type='KW', backwards_cooling=False)

Venus_Driscoll_in = dict(  # table 3
    ident='Venus_Driscoll',
    H_0=2.98e-12,  # final radiogenic heating in W/kg
    L=None,
    Ra_crit_u=660,
    R_p0=6371e3,
    R_c0=3480e3,
    rho_m=4800,
    rho_c=11900,
    k_m=4.2,
    k_lm=10,
    alpha_m=3e-5,
    c_m=1265,
    c_c=840,
    Ea=300e3,
    nu_0=7e7,
    T_s=736,  # fixed surface temp in K
    M_p=4.867e24,  # only used for gravity in this case
    sma=None,
)
Venus_Driscoll_run = dict(T_m0=3000, T_c0=5000, D_l0=300e3, tf=4.5, visc_type='KW', backwards_cooling=False)

Venus_in = dict(
    ident='Venus',  # must match dict name ****_in
    M_p=4.867e24,
    R_p0=6051.8e3,
    R_c0=2900e3,
    Ra_crit_u=450,
    rho_m=3500,
    rho_c=7200,
    c_m=1142,
    c_c=840,
    beta_u=0.335,  # defaults to 1/3
    beta_c=1 / 3,  # defaults to 1/3
    a_rh=2.54,
    k_m=4,
    k_lm=4,
    alpha_m=2.5e-5,
    T_s=467 + 273,
    # viscosity 
    Ea=300e3,
    eta_0=1e21,
    T_ref=1600,
    H_0=23e-12,  # radiogenic heating in W/kg at t0 (Mars value from Thiriet - UPDATE)
    H_f=4.6e-12  # radiogenic heating in W/kg at 4.5 Gyr from Javoy (1999) BSE estimate from CI chondrites
)
Venus_run = dict(
    T_m0=1750, T_c0=2250, D_l0=600e3, tf=4.5, visc_type='Thi', backwards_cooling=False)  # model params

Earth_in = dict(
    ident='Earth',
    CMF=0.3,
    sma=1,
    L=1,
    Alb=0,
    Ra_crit_u=450,
    rho_m=3500,
    rho_c=7200,
    c_m=1142,
    c_c=840,
    beta_u=0.335,  # defaults to 1/3
    beta_c=1 / 3,  # defaults to 1/3
    a_rh=2.54,
    k_m=4,
    k_lm=4,
    alpha_m=2.5e-5,
    # viscosity 
    Ea=300e3,
    eta_0=1e21,
    T_ref=1600,
    H_0=23e-12,  # radiogenic heating in W/kg at t0 (Mars value from Thiriet - UPDATE)
    H_f=4.6e-12  # radiogenic heating in W/kg at 4.5 Gyr from Javoy (1999) BSE estimate from CI chondrites
)
Earth_run = dict(T_m0=1750, T_c0=2000, D_l0=137e3, tf=4.5, visc_type='Thi', backwards_cooling=False)

Moon1_in = dict(
    M_p=7.34767309e22,
    R_p0=1740e3,
    R_c0=390e3,
    Ra_crit_u=450,
    rho_m=3300,
    rho_c=7200,
    c_m=1142,
    c_c=840,
    beta_u=0.346,  # defaults to 1/3
    beta_c=1 / 3,  # defaults to 1/3
    a_rh=2.44,
    k_m=4,
    k_lm=4,
    alpha_m=2.5e-5,
    T_s=250,
    Ea=300e3,
    eta_0=1e21,
    T_ref=1600,
    X_K=83,  # initial abundance of K in wt ppm
    X_U=33e-3,  # initial abundane of U in wt ppm
    X_Th=125e-3,  # initial abundance of Th in wt ppm
    H_f=7e-12,  # radiogenic heating in W/kg at t_f
    ident='Moon1'
)
Moon1_run = dict(T_m0=1750, T_c0=2000, D_l0=445e3, tf=4.5, visc_type='KW', backwards_cooling=False)

Mercury1_in = dict(
    M_p=3.285e23,
    R_p0=2440e3,
    R_c0=2010e3,
    Ra_crit_u=450,
    rho_m=3500,
    rho_c=7200,
    c_m=1142,
    c_c=840,
    beta_u=0.335,  # defaults to 1/3
    beta_c=1 / 3,  # defaults to 1/3
    a_rh=2.54,
    k_m=4,
    k_lm=4,
    alpha_m=2.5e-5,
    T_s=440,
    # viscosity 
    Ea=300e3,
    eta_0=1e21,
    T_ref=1600,
    X_K=400,  # initial abundance of K in wt ppm
    X_U=28e-3,  # initial abundane of U in wt ppm
    X_Th=50e-3,  # initial abundance of Th in wt ppm
    H_f=5e-12,  # radiogenic heating in W/kg at t_f
    ident='Mercury1'
)
Mercury1_run = dict(T_m0=1750, T_c0=2000, D_l0=137e3, tf=4.5, visc_type='KW', backwards_cooling=False)

Mars1_in = dict(
    Alb=None,
    X_K=305,  # initial abundance of K in wt ppm
    X_U=16e-3,  # initial abundane of U in wt ppm
    X_Th=56e-3,  # initial abundance of Th in wt ppm
    L=None,
    Ra_crit_u=450,
    R_p0=3400e3,
    R_c0=1700e3,
    g_sfc=3.7,
    alpha_m=2.5e-5,  # thermal expansivity
    k_m=4,  # silicate thermal conductivity
    CMF=0.24,  # not used
    rho_c=7200,  # core density
    rho_m=3500,  # mantle density
    c_m=1142,  # <----???
    c_c=840,  # specific heat for core in J/K/kg
    k_lm=4,  # thermal conduvtivity lower mantle
    beta_u=1 / 3,  # defaults to 1/3
    beta_c=1 / 3,  # defaults to 1/3
    a_rh=2.44,  # for beta=1/3
    Ea=300e3,  # activation energy in J for viscosity law
    eta_0=1e21,  # reference dynamic viscosity in Pa s
    T_ref=1600,  # viscosity law reference temperature in K
    T_s=250,  # fixed surface temp in K
    M_p=6.39e23,  # only used for gravity in this case
    sma=None,
    ident='Mars1',
    H_0=23e-12,
    H_f=4e-12,  # radiogenic heating in W/kg at t_f
    # backwards_cooling=True,
)
Mars1_run = dict(T_m0=1750, T_c0=2250, D_l0=300e3, tf=4.5, visc_type='Thi', rad_type='H0', backwards_cooling=False,
                 lid_heating=False)

TRAPPIST1e_in = dict(
    ident='TRAPPIST1e',  # must match dict name ****_in
    M_p=0.692 * p.M_E,
    R_p0=0.788 * p.R_E,
    T_s=273,
    CMF=0.246,
    FeMg=0.48,
    rho_m=3300,
    k_m=3,
    alpha_m=3e-5,
    H_0=22.65771894e-12,  # radiogenic heating in W/kg at t0  based on below
    H_f=4.6e-12,  # radiogenic heating in W/kg at 4.5 Gyr from Javoy (1999) BSE estimate from CI chondrites
)
TRAPPIST1e_run = dict(T_m0=1750, T_c0=2250, D_l0=150e3, tf=6, visc_type='KW', backwards_cooling=False)  # model params

TRAPPIST1f_in = dict(
    ident='TRAPPIST1f',  # must match dict name ****_in
    M_p=1.039 * p.M_E,
    R_p0=1.045 * p.R_E,
    T_s=273,
    CMF=0.201,
    FeMg=0.45,
    rho_m=3300,
    k_m=3,
    alpha_m=3e-5,
    H_0=22.65771894e-12,  # radiogenic heating in W/kg at t0  based on below
    H_f=4.6e-12,  # radiogenic heating in W/kg at 4.5 Gyr from Javoy (1999) BSE estimate from CI chondrites
)
TRAPPIST1f_run = dict(T_m0=1750, T_c0=2250, D_l0=150e3, tf=6, visc_type='KW', backwards_cooling=False)  # model params

TRAPPIST1g_in = dict(
    ident='TRAPPIST1g',  # must match dict name ****_in
    M_p=1.321 * p.M_E,
    R_p0=1.129 * p.R_E,
    T_s=273,
    CMF=0.161,
    FeMg=0.34,
    rho_m=3300,
    k_m=3,
    alpha_m=3e-5,
    H_0=22.65771894e-12,  # radiogenic heating in W/kg at t0  based on below
    H_f=4.6e-12,  # radiogenic heating in W/kg at 4.5 Gyr from Javoy (1999) BSE estimate from CI chondrites
)
TRAPPIST1g_run = dict(T_m0=1750, T_c0=2250, D_l0=150e3, tf=6, visc_type='KW', backwards_cooling=False)  # model params

TRAPPIST1h_in = dict(
    ident='TRAPPIST1h',  # must match dict name ****_in
    M_p=0.326 * p.M_E,
    R_p0=0.755 * p.R_E,
    T_s=273,
    CMF=0.165,
    FeMg=0.35,
    rho_m=3300,
    k_m=3,
    alpha_m=3e-5,
    H_0=22.65771894e-12,  # radiogenic heating in W/kg at t0  based on below
    H_f=4.6e-12,  # radiogenic heating in W/kg at 4.5 Gyr from Javoy (1999) BSE estimate from CI chondrites
)
TRAPPIST1h_run = dict(T_m0=1750, T_c0=2250, D_l0=150e3, tf=6, visc_type='KW', backwards_cooling=False)  # model params
