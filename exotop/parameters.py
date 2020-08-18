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

# radioisotope data from Treatise on Geophys, Dye (2012)
# order: [238U, 235U, 232Th, 40K] 
global lambda_n, p_n, K_0, U_0_235, U_0_238, Th_0
lambda_n = [a*1e-9/years2sec for a in [0.15541417, 0.98458406, 0.04951051, 0.55011681]] # decay constant in s^-1
p_n = [95.13e-6, 568.47e-6, 26.3e-6, 28.47e-6] # heating rate in W kg^-1 
K_0 = 0.0117e-2 # ratio of 40-K to total K at time 0, think these ratios are by weight but double check
U_0_235 = 0.0072 # ratio of 235-U to total U at time 0
U_0_238 = 0.9927 # ratio of 238-U to total U at time 0 
Th_0 = 1 # ratio of 232-Th to total Th at time 0 
