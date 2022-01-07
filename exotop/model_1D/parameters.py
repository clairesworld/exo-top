###### PHYSICAL CONSTANTS ######
M_E = 5.972e24  # earth mass in kg
R_E = 6371e3  # earth radius in m
TO = 1.35e9 * 1000 ** 3  # earth ocean volume
L_sun = 3.9e26  # solar luminosity in W
G = 6.67408e-11
years2sec = 31557600
sec2Gyr = 1e-9 / years2sec
AU2m = 1.5e11
R_b = 8.3144598  # universal gas constant in J mol −1 K −1
sb = 5.67e-8  # Stefan Boltzmann constant in W m^-2 K^-4

# other solar system
M_Mars = 6.39e23
R_Mars = 3389.5e3
M_Venus = 4.867e24
R_Venus = 6051.8e3

# radioisotope data from Treatise on Geophys, Dye (2012)
# order: [238U, 235U, 232Th, 40K]
lambda_n = [a * 1e-9 / years2sec for a in [0.15541417, 0.98458406, 0.04951051, 0.55011681]]  # decay constant in s^-1
p_n = [95.13e-6, 568.47e-6, 26.3e-6, 28.47e-6]  # heating rate in W kg^-1
U_0_238 = 0.9927  # ratio of 238-U to total U at time 0 (natural abundance)
U_0_235 = 0.0072  # ratio of 235-U to total U at time 0
Th_0 = 1  # ratio of 232-Th to total Th at time 0
K_0 = 0.0117e-2  # ratio of 40-K to total K at time 0, think these ratios are by weight but double check

# topography scaling relationship
beta_h = [9.265808299836122, -0.5662122195385831, -1.466326823555549, 0.07316787161626584]
cov_beta_h = [[ 7.78051245e+00, -4.34305471e-01, -9.95308352e-01,  5.55550360e-02],
 [-4.34305471e-01,  2.47533957e-02,  5.55549078e-02, -3.16618176e-03],
 [-9.95308352e-01,  5.55549078e-02,  1.27403161e-01, -7.11091023e-03],
 [ 5.55550360e-02, -3.16618176e-03 ,-7.11091023e-03,  4.05240762e-04]]



