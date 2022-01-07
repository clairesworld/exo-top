""" https://zenodo.org/record/5560138 """
import numpy as np
from scipy.optimize import fsolve, bisect
from scipy.special import erf
from model_1D import thermal as th
import matplotlib.pyplot as plt

def strength_profile(g, Tm, Ts, rho, R_p=None, R_l=None, T_l=None, a0=None, p=1/2, q=2, lith_age=320, wet=False, plot=True):
    # rho=3500  # crust/mantle (average) density

    # dislocation creep parameters
    n_dis=3.5;
    A_0=1.1e5; ## Stress should be in MPa
    E=520000; # J/mol
    V=10e-6; # m^3/mol (could be 5-20e-6)
    R=8.314; # gas constant

    G0=77.4e9; # Pa
    G_prime=1.61;

    if p == 1/2 and q == 2:
        # For p=1/2 and q=2; from Jain et al 2017 case a
        E_lt=452*1000; # J/mol
        V_lt=8e-6; # could be negative, but small number
        A_lt=10**-7.78;
        sigma_p=7.32e9; # Pa

    elif p == 1 and q == 2:
        # Low temperature plasticity based on Jain et al 2017 for p=1 and q=2; case a
        E_lt=225*1000; # J/mol
        V_lt=3e-6; # could be negative, but small number
        A_lt=10**-7.69;
        sigma_p=6.20e9; # Pa

    e_dis=1e-13; # 1/s  # strain rate
    k=4  #5; # Thermal conductivity

    # g=unifrnd(5,40); # gravity
    if R_l is None:
        age=lith_age  #unifrnd(5,320); # Million years lith age
        # Ts=unifrnd(200,1000);
        # Tm=unifrnd(1500,2000);
        delta=2*np.sqrt(31.5*age);
        z = np.linspace(0,delta*1000,5000); # Depth to 200 km in units of meters

        # Error function temperature profile
        T = Ts + (Tm - Ts) * erf((z / 1000) / (2 * np.sqrt(31.5 * age)));

        # Heat flux
        qs = k * (Tm - Ts) * (1 / 1000) * (1 / (np.pi * 31.5 * age)) ** (1 / 2);

    else:
        z=np.linspace(0, R_p - R_l)  #np.linspace(0,delta*1000,5000); # Depth to 200 km in units of meters
        T = th.sph_conduction(r=R_p - z, k_m=k, T_l=T_l, T_s=Ts, R_p=R_p, R_l=R_l, a0=a0)
        qs = th.sph_flux(r=R_p - z, a0=a0, k_m=k, T_l=T_l, T_s=Ts, R_p=R_p, R_l=R_l)
        print('max z', np.max(z)*1e-3, 'km')

    # Pressure profile
    pressure=rho*g*z

    # Re-doing fault strength like in Kohlstedt 1995
    v_stress_grad=34*(g/9.8); #MPa/km
    if wet:
        P_p_grad= 1000*g*1000/1e6; # MPa/km; set to 0 if no pore pressure
    else:
        P_p_grad = 0
    fault_strength1=4.9*(v_stress_grad*z/1000-P_p_grad*z/1000)+P_p_grad*z/1000-v_stress_grad*z/1000
    fault_strength2=3.1*(v_stress_grad*z/1000-P_p_grad*z/1000)+210+P_p_grad*z/1000-v_stress_grad*z/1000
    fault_strength=np.minimum(fault_strength1,fault_strength2);

    # compare classic Byerlee law instead
    sigma_s = np.zeros_like(pressure)
    for ii, pp in enumerate(pressure):
        if pp < 200e6:
            sigma_s[ii] = 0.85 * pp
        else:
            sigma_s[ii] = 50 + 0.6 * pp
    # fault_strength = sigma_s*1e-6  # in MPa
    # print('Byerlee', sigma_s*1e-6, 'MPa')

    flow_strength = np.zeros_like(fault_strength)
    for i in range(len(z)):
        sigma_p1=sigma_p*(1+G_prime*pressure[i]/G0);

        def duct_func(x, p_Pa):
            return A_lt * (x / 1e6) ** 2 * np.exp(
                -((E_lt + p_Pa * V_lt) / (R * T[i])) * (1 - (x / sigma_p1) ** p) ** q) + A_0 * np.exp(-(E + p_Pa * V) / (R * T[i])) * (x / 1e6) ** n_dis - e_dis

        zero = bisect(duct_func, 0, sigma_p1, args=(pressure[i]), xtol=2e-12, rtol=8.881784197001252e-16, maxiter=100, full_output=False, disp=True)
        # zero = fsolve(duct_func,[0, sigma_p1], args=pressure[i])
        flow_strength[i]=zero  # in MPa I think

    # print('\nviscous', flow_strength*1e-6, 'MPa')
    strength=np.minimum(fault_strength, flow_strength/1e6)

    i_bdt= np.argmax(strength) #find(strength==max(strength));
    z_bdt=z[i_bdt]/1e3;  # depth of bdt in km
    sigma_bdt = strength[i_bdt]  # strength at bdt

    if plot:
        plt.plot(strength, z*1e-3, 'k-', label='min strength')
        plt.plot(flow_strength*1e-6, z*1e-3, 'r--', label='viscous strength')
        plt.plot(sigma_s*1e-6, z * 1e-3, 'b--', label="friction: Byerlee's law", lw=1)
        plt.plot(fault_strength, z * 1e-3, 'g--', label='friction: Kohlstedt+ 1995', lw=0.5)
        # plt.plot(fault_strength2, z * 1e-3, 'g--', lw=1)
        plt.axvline(sigma_bdt, c='k', alpha=0.2, lw=0.5)
        plt.xlabel('Strength (MPa)')
        plt.ylabel('Depth (km)')
        plt.legend()
        plt.gca().invert_yaxis()
        plt.xlim([0, 3000])
    return z, strength*1e6   # in m, Pa


def z_brittle(z, strength):
    i_bdt= np.argmax(strength) #find(strength==max(strength));
    z_bdt=z[i_bdt]/1e3;  # depth of bdt in km
    sigma_bdt = strength[i_bdt]  # strength at bdt
    return z_bdt, sigma_bdt


def max_topo(Y, g, rho, C=1/3):
    return (C ** -1 * Y) / (rho * g)


# test_depth, test_sigma =  brittle_ductile_transition(3.804447e+01, 1.584532e+03, 2.077603e+02, 3500, R_p=None, R_l=None, T_l=None, a0=None, p=1/2, q=2, lith_age=7.959980e+01)
# print('BDT depth',test_depth*1e-3, 'km', 'sigma', test_sigma*1e-6, 'MPa')
# plt.show()


z, Y = strength_profile(9.8, Tm=1650, Ts=300, rho=2700, R_p=None, R_l=None, T_l=None, a0=None, p=1/2, q=2, lith_age=320, wet=True, plot=True)

Y_max

### iterate to find h