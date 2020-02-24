Venus_in = dict(
    ident = 'Venus', # must match dict name ****_in
    M_p = 4.867e24,
    R_p0 = 6051.8e3,
    R_c0 = 2900e3,
    Ra_crit_u = 450,
    rho_m = 3500,
    rho_c = 7200,
    c_m = 1142,
    c_c = 840,
    beta_u = 0.335, # defaults to 1/3
    beta_c = None, # defaults to 1/3
    a_rh = 2.54 ,
    k_m = 4,
    k_lm = 4,
    alpha_m = 2.5e-5,
    T_s = 467+273,
    # viscosity 
    Ea = 300e3,
    eta_0 = 1e21,
    T_ref = 1600,
    H_0 = 4.6e-12 # radiogenic heating in W/kg at 4.5 Gyr from Javoy (1999) BSE estimate from CI chondrites
)
Venus_run = dict(
    T_m0=1750, T_c0=2250, D_l0=600e3, tf=4.5, visc_type='Thi', complexity=3)  # model params

Earth_in = dict(
    ident = 'Earth',
    CMF = 0.3,
    sma = 1,
    L=1,
    Alb=0,
    Ra_crit_u = 450,
    rho_m = 3500,
    rho_c = 7200,
    c_m = 1142,
    c_c = 840,
    beta_u = 0.335, # defaults to 1/3
    beta_c = None, # defaults to 1/3
    a_rh = 2.54,
    k_m = 4,
    k_lm = 4,
    alpha_m = 2.5e-5,
    # viscosity 
    Ea = 300e3,
    eta_0 = 1e21,
    T_ref = 1600,
    H_0 = 4.6e-12 # radiogenic heating in W/kg at 4.5 Gyr from Javoy (1999) BSE estimate from CI chondrites
)
Earth_run = dict(T_m0=1750, T_c0=2000, D_l0=137e3, tf=4.5, visc_type='Thi', complexity=3)

Moon1_in = dict(
    M_p = 7.34767309e22,
    R_p0 = 1740e3,
    R_c0 = 390e3,
    Ra_crit_u = 450,
    rho_m = 3300,
    rho_c = 7200,
    c_m = 1142,
    c_c = 840,
    beta_u = 0.346, # defaults to 1/3
    beta_c = None, # defaults to 1/3
    a_rh = 2.44,
    k_m = 4,
    k_lm = 4,
    alpha_m = 2.5e-5,
    T_s = 250,
    Ea = 300e3,
    eta_0 = 1e21,
    T_ref = 1600,
    X_K = 83, # initial abundance of K in wt ppm 
    X_U = 33e-3, # initial abundane of U in wt ppm 
    X_Th = 125e-3, # initial abundance of Th in wt ppm 
    H_0 = 7e-12, # radiogenic heating in W/kg at t_f
    ident = 'Moon1'
    )
Moon1_run = dict(T_m0=1750, T_c0=2000, D_l0=445e3, tf=4.5, visc_type='Thi', complexity=3)

Mercury1_in = dict(
    M_p = 3.285e23,
    R_p0 = 2440e3,
    R_c0 = 2010e3,
    Ra_crit_u = 450,
    rho_m = 3500,
    rho_c = 7200,
    c_m = 1142,
    c_c = 840,
    beta_u = 0.335, # defaults to 1/3
    beta_c = None, # defaults to 1/3
    a_rh = 2.54, 
    k_m = 4,
    k_lm = 4,
    alpha_m = 2.5e-5,
    T_s = 440,
    # viscosity 
    Ea = 300e3,
    eta_0 = 1e21,
    T_ref = 1600,
    X_K = 400, # initial abundance of K in wt ppm 
    X_U = 28e-3, # initial abundane of U in wt ppm 
    X_Th = 50e-3, # initial abundance of Th in wt ppm 
    H_0 = 5e-12, # radiogenic heating in W/kg at t_f
    ident='Mercury1'
    )
Mercury1_run = dict(T_m0=1750, T_c0=2000, D_l0=137e3, tf=4.5, visc_type='Thi', complexity=3)

Mars1_in = dict(
     Alb=None, 
     H_0=4e-12, # final radiogenic heating in W/kg
     X_K=305, # initial abundance of K in wt ppm
     X_U=16e-3, # initial abundane of U in wt ppm 
     X_Th=56e-3, # initial abundance of Th in wt ppm 
     L=None, 
     Ra_crit_u=450, 
     R_p0=3400e3, 
     R_c0=1700e3,
     alpha_m=2.5e-5, # thermal expansivity
     k_m=4, # silicate thermal conductivity
     CMF=0.24, # not used
     rho_c=7200, # core density
     rho_m=3500, # mantle density 
     rho_lith=None, 
     c_m=1142, #<----??? TODO: check if you need constant volume c_p
     c_c=840, # specific heat for core in J/K/kg
     k_lm=4, # thermal conduvtivity lower mantle
     beta_u=None, # defaults to 1/3
     beta_c=None, # defaults to 1/3 
     a_rh=2.44, # for beta=1/3 
     Ea=300e3, # activation energy in J for viscosity law
     eta_0=1e21, # reference dynamic viscosity in Pa s
     T_ref=1600, # viscosity law reference temperature in K
     T_s=250, # fixed surface temp in K
     M_p=6.39e23, # only used for gravity in this case
     sma=None, 
     ident='Mars1'
     )
Mars1_run =  dict(T_m0=1750, T_c0=2250, D_l0=300e3, tf=4.5, visc_type='Thi', complexity=3)