#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# %%
"""
Created on Tue Mar 29 16:58:54 2022

@author: evgeni grishin
"""
import pandas as pd
import numpy as np
import math
import matplotlib
import matplotlib.colors as colors
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from matplotlib import rc
from scipy.optimize import fsolve
from astropy import constants as const
from scipy.io import savemat
# rc('text', usetex=False)
# %%
# CONSTANTS
c          = const.c.cgs.value                # speed of light (c) in cgs
msun       = const.M_sun.cgs.value            # solar mass (Msun) in cgs
G          = const.G.cgs.value                # gravitational constant (G) in cgs
au         = const.au.cgs.value               # astronomical unit (AU) in cgs
rsun       = const.R_sun.cgs.value            # solar radius (Rsun) in cgs
rsun_to_au = const.R_sun.value/const.au.value # solar radius in au
day_to_s   = 86400                            # day in s


# %%
# FUNCTIONS
def period_ratio(m1, m2, m3, a1, a2):
    # m1, m2 and m3 are the masses of the primary, secondary and tertiary, respectively
    # a1 and a2 are the semi-major axes of the inner and outer orbit, respectively
    # given that we calculate a ratio, the values of the variables can be in any units, as long as they are consistent with one another
    m_bin = m1 + m2
    m_tot = m_bin + m3
    return (a1/a2)**1.5 * (m_tot/m_bin)**0.5

# second timescale - 1st order (no dependence on inclination)
def t_sec(m1, m2, m3, a1, a2, e2):
    # [m1]=[m2]=[m3]=Msun
    # [a1]=[a2]=au
    # [t_sec]=s
    m1 = m1 * msun
    m2 = m2 * msun
    m3 = m3 * msun
    a1 = a1 * au
    a2 = a2 * au

    mu_in = m1*m2/(m1+m2)
    m_bin = (m1+m2) 
    L1 = mu_in * (G * m_bin * a1)**0.5
    C = G * mu_in * m3/8/a2/(1-e2**2)**1.5 * (a1/a2)**2
    
    return L1 / 6 / C
    
# how bad the double averaging is (e.g. Grishin, Perets and Fragione 2018; Luo et al., 2016)
def epsilon_sa(m1, m2, m3, a1, a2, e2):
    # e2 is the eccentricity of the outer orbit
    # [e2]=adim
    m_bin = m1+m2
    m_tot = m_bin + m3
    return (a1/a2/(1-e2**2))**1.5*m3/(m_tot*m_bin)**0.5

# angular momentum ratio (e.g. Haim and Katz 2018)
def eta(m1, m2, m3, a1, a2, e2):
    m_bin = m1+m2
    m_tot = m_bin + m3
    mu_in = m1*m2/m_bin
    mu_out = m3*m_bin/m_tot
    return mu_in/mu_out*(m_bin*a1/m_tot/a2/(1-e2**2))**0.5

# importance of GR corrections (e.g. Liu, Munoz and Lai 2015)
def epsilon_gr(m1, m2, m3, a1, a2, e2):
    # [m1]=[m2]=[m3]=Msun
    # [a1]=[a2]=au
    # [e2]=adim
    m_bin = m1+m2
    rg_over_a1= G*m_bin*msun/c**2/(a1*au)
    return 3 * m_bin*((1-e2**2)**1.5)/m3*(a2/a1)**3 * rg_over_a1


# importance of tidal corrections (e.g. Liu, Muñoz and Lai 2015)
def epsilon_tide(m1, m2, m3, a1, a2, e2, r1, r2, k1, k2):
    # [m1]=[m2]=[m3]=Msun
    # [a1]=[a2]=au
    # r1 and r2 are the radius of the primary and the secondary, respectively
    # [r1]=[r2]=au
    # k1 and k2 are the apsidal motion constants of the primary and the secondary, respectively
    # [k1]=[k2]=adim
    m_bin = m1+m2
    b2 = a2 * (1-e2**2)**0.5
    return 15 * m_bin * b2**3 * (m1**2 * k1 * r1**5 + m2**2 * k2 * r2**5 ) / a1**8 / m1 / m2/ m3

# octupole term (e.g. Naoz 2016 and many others)
def epsilon_oct(m1, m2, a1, a2, e2):
    return (m1-m2)/(m1+m2)* (a1/a2) * e2/(1-e2**2)

# ALEJANDRO: include a function to check for stability?
# m_mardling = [140*(1 + (x/2.8/a1)**2.5 ) for x in a2]
# def f_tory(q):
#     return 10**(-0.6+0.04*q) * q **(0.32+0.1*q)

# %%
# calculate the maximal eccentricity according to Mangipudi et al. 2022 and the appendix with tides
def get_maximal_eccentricity(m1, m2, m3, a1, a2, e2, eps_SA1, eps_GR1, eps_Tide1, cos_inc, debug):
    # [m1]=[m2]=[m3]=Msun
    # [a1]=[a2]=au
    # [e2]=adim
    # [r1]=[r2]=au
    # [k1]=[k2]=adim
    # [cos_inc]=adim
    cosi_0    = cos_inc
    ETA       = eta(m1, m2, m3, a1, a2, e2)
    
    if debug:
        print("ETA=",ETA)
        print("eps_GR1=",eps_GR1)
        print("eps_SA1=",eps_SA1)  
        print("eps_Tide1=",eps_Tide1)          
    
    # in the next function, what does F means?
    def F(e_m=0.5,cosi_0=cosi_0, eta=ETA):
        return (1 - e_m**2)*(-3 + eta*cosi_0 + e_m**2*eta**2/4) + 5*(cosi_0 + eta*e_m**2/2)**2 - eps_SA1*(3/8)*( 
            -((1 - e_m**2)*cosi_0**3/e_m**2) - (eta/2)*(1 - e_m**2) + (cosi_0 + e_m**2*eta/2)**3*(1/e_m**2 - 16) 
            - 9*(1 - e_m**2)*(cosi_0 + e_m**2*eta/2)) - (eps_GR1*8/3)*((1 - e_m**2)/e_m**2)*(1 - (1 - e_m**2)**(-1/2)) \
            - eps_Tide1*(8/45)*(1 - (1 + 3*e_m**2 + (3/8)*e_m**4)/( (1 - e_m**2)**(9/2) )  )

    # What does a mean? Isn't it supposed to be e_max or e_max^2 or something?
    a = fsolve(F,0.999)
    
    def dH(e_m, eta = ETA, eps_SA = eps_SA1, cosi_0 = cosi_0):
        return ( ((1 + (cosi_0 + eta/2)*eta)/(np.sqrt(3/5) + np.sqrt(1 - e_m**2)*eta))*(9/8)*eps_SA*e_m  )* \
            ( np.sqrt(1 - e_m**2) - ((2*e_m**2 - 1)/2)*(((1 + (cosi_0 + eta/2)*eta)/(np.sqrt(3/5) + np.sqrt(1 - e_m**2)*eta))*(9/8)*eps_SA))
    da = dH(e_m=0.5)
    
    if a < 0:
        a=0
    if a+da >= 1:
        return 1
    return a+da


# %%
def calculate_grids(NN, Ni, period, m1, m2, r1, r2, r1_core, r2_core, k1, k2, incl_SA, incl_GR, incl_Tides, debug):
    # NN x NN is the size of the outer grid of masses and separations
    # Ni is the number of uniform realisation in cos(inclination)
    # period is the orbital period of the inner binary
    # [period]=days
    # [m1]=[m2]=Msun
    # [r1]=[r2]=Rsun
    # [r1_core]=[r2_core]=Rsun    
    # [k1]=[k2]=adim
    # incl_SA is a boolean that if true calculates epsilon_SA and otherwise sets it to zero
    # incl_GR is a boolean that if true calculates incl_GR and otherwise sets it to zero
    # incl_Tide is a boolean that if true calculates incl_Tide and otherwise sets it to zero
    e2 = 0.0
    a1 = (G*(m1+m2)*msun * (period * day_to_s)**2/4/np.pi/np.pi)**0.3333/au # [a1]=au
    r1 = r1 * rsun / au; r2 = r2 * rsun / au # [r1]=[r2]=au
    r1_core = r1_core * rsun / au; r2_core = r2_core * rsun / au # [r1]=[r2]=au    
    m3 = np.logspace(0,np.log10(200),NN) # [m3]=Msun
#     a2 = np.logspace(np.log10(a1) + np.log10(2), np.log10(a1)+ np.log10(300), NN)    
    # ALEJANDRO: what is with the np.log10(2) in the previous, commented line?
    a2 = np.logspace(np.log10(0.1), np.log10(10), NN)    
    eccs = np.zeros(shape=[NN,NN])
    # ALEJANDRO: what is rps supposed to stand for? I know is RPeriapsiS?
    rps = np.zeros(shape=[NN,NN])
    fraction = np.zeros(shape=[NN,NN])
    cos_inc_for_max_ecc = np.zeros(shape=[NN,NN])
    cos_incs=np.linspace(-1,1,Ni)
    tertiary_mass = np.zeros(shape=[NN,NN])
    outer_separation = np.zeros(shape=[NN,NN])
    secular_timescale = np.zeros(shape=[NN,NN])
    eps_SA = np.zeros(shape=[NN,NN])    
    eps_GR = np.zeros(shape=[NN,NN])
    eps_Tides = np.zeros(shape=[NN,NN])
    ETA = np.zeros(shape=[NN,NN])
    
    if debug:
        print("r1 and r2:",r1,r2)
        print("r1_core and r2_core:",r1_core,r2_core)
        print("a1, a2min, a2max:",a1,np.min(a2),np.max(a2))
        print("m3min and m3max:",np.min(m3),np.max(m3))

    for i in range(0,NN):
        if debug:
            print (i)
        for j in range(0,NN):
            tertiary_mass[i,j] = m3[i]
            outer_separation[i,j] = a2[j]
            secular_timescale[i,j] = t_sec(m1, m2, m3[i], a1, a2[j], e2)
            ETA[i,j] = eta(m1, m2, m3[i], a1, a2[j], e2)
            
            if incl_SA:
                eps_SA[i,j] = epsilon_sa(m1, m2, m3[i], a1, a2[j], e2)
            
            if incl_GR:
                eps_GR[i,j] = epsilon_gr(m1, m2, m3[i], a1, a2[j], e2)
            
            if incl_Tides:
                eps_Tides[i,j] = epsilon_tide(m1, m2, m3[i], a1, a2[j], e2, r1_core, r2_core, k1, k2)
            
            # Evgeni, I am not sure we need the following lines here
            # eccs[i,j] = get_maximal_eccentricity(m1, m2, m3[i], a1, a2[j], e2, eps_SA[i,j], eps_GR[i,j], eps_Tides[i,j], cos_inc, False) 
            # rps[i,j] = a1 * (1 - eccs[i,j]) * au / rsun # ALEJANDRO: Why does rps needs to be in Rsun?
            # What we want is to get the maximum eccentricity, minimum periapsis, and inclination for each m3/aout combo
            # I think this is better achieve within the next for loop

            counter = 0
            ec = 0.0
            rpp = a1
            cos_inc = 0.0
            for k in range(0,Ni):
                ec_temp = get_maximal_eccentricity(m1, m2, m3[i], a1, a2[j], e2, eps_SA[i,j], eps_GR[i,j], eps_Tides[i,j], cos_incs[k], False)
                rpp_temp = a1 * (1 - ec) 
                cos_inc_temp = cos_incs[k]
                
                # Here we check what is the maximum eccentricity in the cos(i) parameter space
                if ec_temp>ec:
                    ec = ec_temp
                    rpp = rpp_temp
                    cos_inc = cos_inc_temp
                # The next condition is for L2 filling and only valid for equal-mass binaries
                # the system enters a contact phase if the radius reaches the second lagrangian point
                # L2 = 1.32*R_{RL}, which for an equal-mass binary results in
                # L2 = 0.5002r_p = 0.5002*a*(1-e) (we check for contact phase at periastron)
                if 0.5002*rpp_temp/r1 < 1:
                    counter=counter+1

            eccs[i,j] = ec
            # Double check rps is ploting the correct quantity
            rps[i,j] = rpp
            cos_inc_for_max_ecc[i,j] = cos_inc
            fraction[i,j] = counter/Ni       

    return tertiary_mass, outer_separation, eccs, rps, cos_inc_for_max_ecc, fraction, secular_timescale, ETA, eps_SA, eps_GR, eps_Tides


# %%
# CALCULATE
# NN=40; Ni=30; 
NN=100; 
Ni=30; 
metallicity = 0.0001
is_CHE = 1
# period_days=1; m1=55; m2=55; r1=7; r2=7; k1=0.0; k2=0.0;
period_days=2; m1=55; m2=55; r1=7; r2=7; r1_core=7; r2_core=7; k1=0.024; k2=0.024;
eps_SA = True
eps_GR = True
eps_Tides = True

string_to_save = "triple_Z="+str(metallicity)+"_CHE="+str(is_CHE)+"_M1=M2="+str(m2)+"_Porb="+str(period_days)

if eps_SA:
    string_to_save+="_SA"
    
if eps_GR:
    string_to_save+="_GR"   
    
if eps_Tides:
    string_to_save+="_Tides"   
    
string_to_save+=".mat"
print(string_to_save)

m3, a2, eccs, rps, cos_inc, fraction, tau_ZKL, ETA, eps_SA, eps_GR, eps_Tides = calculate_grids(NN, Ni, period_days, m1, m2, r1, r2, r1_core, r2_core, k1, k2, eps_SA, eps_GR, eps_Tides, False)

mdic = {"m3":m3, "a2":a2, "eccs":eccs, "rps":rps, "cos_inc":cos_inc, "fraction":fraction, "tau_ZKL":tau_ZKL, "ETA":ETA, "eps_SA":eps_SA, "eps_GR":eps_GR, "eps_Tides":eps_Tides}
savemat(string_to_save, mdic)

plot_flag = True
if plot_flag:
    # PLOT
    plt.rc('text', usetex=False)
    matplotlib.rcParams.update({'font.size': 20})
    plt.contourf(m3, a2, eccs)
    plt.xlabel(r'$m_3 / \rm M_\odot$')
    plt.ylabel(r'$a_2 / \rm au$')
    plt.xscale('log')
    plt.yscale('log')
    plt.text(2.03, 1., r'$e_{\rm max}$')
    plt.colorbar()

    plt.figure(2)
    plt.contourf(m3, a2, fraction) #[2,2.5,3,3.5,4,4.5])
    plt.xlabel(r'$m_3 / \rm M_\odot$')
    plt.ylabel(r'$a_2 / \rm au$')
    plt.xscale('log')
    plt.yscale('log')
    plt.text(2.02, 1., 'fraction')        
    plt.colorbar()

# %%
NN=100; 
Ni=30; 
period_days=2; m1=55; m2=55; r1=7; r2=7; r1_core=7; r2_core=7; k1=0.024; k2=0.024;
a1 = (G*(m1+m2)*msun * (period_days * day_to_s)**2/4/np.pi/np.pi)**0.3333/au # [a1]=au
m3 = 1.0
a2 = 2.0
e2 = 0.0

# eps_SA1 = epsilon_sa(m1, m2, m3, a1, a2, e2)
eps_SA1 = 0
eps_Tide1 = 0
debug = False

# eps_GR1 = epsilon_gr(m1, m2, m3, a1, a2, e2)
eps_GR1 = 2.05

print(eps_SA1)
print(eps_GR1)
print(get_maximal_eccentricity(m1, m2, m3, a1, a2, e2, eps_SA1, eps_GR1, eps_Tide1, 0, debug))
print(get_maximal_eccentricity(m1, m2, m3, a1, a2, e2, eps_SA1, eps_GR1, eps_Tide1, -1, debug))
print(get_maximal_eccentricity(m1, m2, m3, a1, a2, e2, eps_SA1, eps_GR1, eps_Tide1, +1, debug))

# %%
