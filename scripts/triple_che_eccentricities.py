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
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from matplotlib import rc
from scipy.optimize import fsolve
import matplotlib.colors as colors
rc('text', usetex=False)
from astropy import constants as const
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
def get_maximal_eccentricity(m1, m2, m3, a1, a2, e2, r1, r2, k1, k2, cos_inc, debug):
    # [m1]=[m2]=[m3]=Msun
    # [a1]=[a2]=au
    # [e2]=adim
    # [r1]=[r2]=au
    # [k1]=[k2]=adim
    # [cos_inc]=adim
    cosi_0    = cos_inc
    ETA       = eta(m1, m2, m3, a1, a2, e2)
    eps_GR1   = epsilon_gr(m1, m2, m3, a1, a2, e2)
    eps_SA1   = epsilon_sa(m1, m2, m3, a1, a2, e2)
    eps_Tide1 = epsilon_tide(m1, m2, m3, a1, a2, e2, r1, r2, k1, k2)
    
    if debug:
        print("ETA=",ETA)
        print("eps_GR1=",eps_GR1)
        print("eps_SA1=",eps_SA1)  
        print("eps_Tide1=",eps_Tide1)          
    
    def F(e_m=0.5,cosi_0=cosi_0, eta=ETA):
        return (1 - e_m**2)*(-3 + eta*cosi_0 + e_m**2*eta**2/4) + 5*(cosi_0 + eta*e_m**2/2)**2 - eps_SA1*(3/8)*( 
            -((1 - e_m**2)*cosi_0**3/e_m**2) - (eta/2)*(1 - e_m**2) + (cosi_0 + e_m**2*eta/2)**3*(1/e_m**2 - 16) 
            - 9*(1 - e_m**2)*(cosi_0 + e_m**2*eta/2)) - (eps_GR1*8/3)*((1 - e_m**2)/e_m**2)*(1 - (1 - e_m**2)**(-1/2)) \
            - eps_Tide1*(8/45)*(1 - (1 + 3*e_m**2 + (3/8)*e_m**4)/( (1 - e_m**2)**(9/2) )  )
    a = fsolve(F,0.999)
    
    def dH(e_m, eta = ETA, eps_SA = eps_SA1, cosi_0 = cosi_0):
        return ( ((1 + (cosi_0 + eta/2)*eta)/(np.sqrt(3/5) + np.sqrt(1 - e_m**2)*eta))*(9/8)*eps_SA*e_m  )* \
            ( np.sqrt(1 - e_m**2) - ((2*e_m**2 - 1)/2)*(((1 + (cosi_0 + eta/2)*eta)/(np.sqrt(3/5) + np.sqrt(1 - e_m**2)*eta))*(9/8)*eps_SA   )  )
    da = dH(e_m=0.5)
    
    if a < 0:
        a=0
    if a+da >= 1:
        return 1
    return a+da


# %%
# # Testing data. Probably should be deleted when everything is working.
# NN=20; Ni=15; period_days=0.8; m1=45; m2=45; r1=7; r2=7; k1=0.02; k2=0.02;
# ain = 0.07559;
# print(get_maximal_eccentricity(m1, m2, 45, ain, 0.3, 0, r1*rsun_to_au, r2*rsun_to_au, k1, k2, 0, 1))
# print(get_maximal_eccentricity(m1, m2, 45, ain, 0.3, 0, r1*rsun_to_au, r2*rsun_to_au, k1, k2, -1, 1))
# print(get_maximal_eccentricity(m1, m2, 45, ain, 0.3, 0, r1*rsun_to_au, r2*rsun_to_au, k1, k2, 1, 1))
# print()
# print(get_maximal_eccentricity(m1, m2, 90, ain, 0.3, 0, r1*rsun_to_au, r2*rsun_to_au, k1, k2, 0, 1))
# print(get_maximal_eccentricity(m1, m2, 90, ain, 0.3, 0, r1*rsun_to_au, r2*rsun_to_au, k1, k2, -1, 1))
# print(get_maximal_eccentricity(m1, m2, 90, ain, 0.3, 0, r1*rsun_to_au, r2*rsun_to_au, k1, k2, 1, 1))
# print()
# print(get_maximal_eccentricity(m1, m2, 90, ain, 5, 0, r1*rsun_to_au, r2*rsun_to_au, k1, k2, 0, 1))
# print(get_maximal_eccentricity(m1, m2, 90, ain, 5, 0, r1*rsun_to_au, r2*rsun_to_au, k1, k2, -1, 1))
# print(get_maximal_eccentricity(m1, m2, 90, ain, 5, 0, r1*rsun_to_au, r2*rsun_to_au, k1, k2, 1, 1))
# print()
# print(get_maximal_eccentricity(m1, m2, 90, ain, 10, 0, r1*rsun_to_au, r2*rsun_to_au, k1, k2, 0, 1))
# print(get_maximal_eccentricity(m1, m2, 90, ain, 10, 0, r1*rsun_to_au, r2*rsun_to_au, k1, k2, -1, 1))
# print(get_maximal_eccentricity(m1, m2, 90, ain, 10, 0, r1*rsun_to_au, r2*rsun_to_au, k1, k2, 1, 1))

# %%
# ALEJANDRO: this function I still need to check/clean
# plot e_max for a particular choice of inclination cosinc - set e_maxplot = True
# plot fractions of collisoins - uniform in cos inc with Ni inclinations
# set fraction_plot = True
# the outer grid of masses and separations is of size NN x NN

#the system is colliding if the pericentre r_p < 1.3 * (r1+r2)
# ALEJANDRO: double check this previous condition
# the system enters a contact phase if the radius reaches the second lagrangian point
# L2 = 1.32*R_{RL}, which for an equal-mass binary results in
# L2 = 0.5002r_p = 0.5002*a*(1-e) (we check for contact phase at periastron)
# therefore, the condition is r1=r2=r >= 0.5002r_p
# 0.5002r_p/r1 < 1

def plot_grids(NN, Ni, period, m1, m2, r1, r2, k1, k2, cos_inc, debug):
    # NN and Ni are ...
    # period is the orbital period of the inner binary
    # [period]=days
    # [m1]=[m2]=Msun
    # [r1]=[r2]=Rsun
    # [k1]=[k2]=adim
    # [cos_inc]=adim
    # AVG: explain all things here
    # AVG: probably better a script to first calculate the grid and another to then plot
    a1 = (G*(m1+m2)*msun * (period * day_to_s)**2/4/np.pi/np.pi)**0.3333/au # [a1]=au?
    r1 = r1 * rsun / au; r2 = r2 * rsun / au # [r1]=[r2]=au
    m3 = np.logspace(0,2,NN) # [m3]=Msun
#     a2 = np.logspace(np.log10(a1) + np.log10(2), np.log10(a1)+ np.log10(300), NN)    
    # ALEJANDRO: what is with the np.log10(2)?
    a2 = np.logspace(np.log10(0.1), np.log10(10), NN)    
    eccs = np.zeros(shape=[NN,NN])
    # ALEJANDRO: what is rps?
    rps = np.zeros(shape=[NN,NN])
    fraction = np.zeros(shape=[NN,NN])
    cos_incs=np.linspace(-1,1,Ni)
    
    tertiary_mass = np.zeros(shape=[NN,NN])
    outer_separation = np.zeros(shape=[NN,NN])
    
    if debug:
        print("r1 and r2:",r1,r2)
        print("a1, a2min, a2max:",a1,np.min(a2),np.max(a2))
        print("m3min and m3max:",np.min(m3),np.max(m3))

    for i in range(0,NN):
        if debug:
            print (i)
        for j in range(0,NN):
            # ALEJANDRO: what is the next lines for? I don't get why we need to user define a cos_inc 
            # if later we will explore the cos_inc parameter space 
            tertiary_mass[i,j] = m3[i]
            outer_separation[i,j] = a2[j]
            eccs[i,j] = get_maximal_eccentricity(m1, m2, m3[i], a1, a2[j], 0, r1 ,r2, k1, k2, cos_inc, False) 
            rps[i,j] = a1 * (1 - eccs[i,j]) * au / rsun # ALEJANDRO: Why does rps needs to be in Rsun?
            counter=0
            for k in range(0,Ni):
                ec = get_maximal_eccentricity(m1, m2, m3[i], a1, a2[j], 0, r1 ,r2, k1, k2, cos_incs[k], False)
                rpp = a1 * (1 - ec) 
#                 print("rpp=",rpp)
#                 if rpp/1.3/(r1+r2)<=1:
                if 0.5002*rpp/r1 < 1:
#                 if r1 >= 0.5002*rpp or r2 >= 0.5002*rpp:
                    counter=counter+1

            fraction[i,j] = counter/Ni       

    return tertiary_mass, outer_separation, eccs, rps, fraction


# %%
# CALCULATE
NN=40; Ni=30; 
period_days=1; m1=45; m2=45; r1=7; r2=7; k1=0.02; k2=0.02;
m3, a2, eccs, rps, fraction = plot_grids(NN, Ni, period_days, m1, m2, r1, r2, k1, k2, 0.0, False)

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
#         plt.subplots_adjust(left=0.18, bottom=0.15, right=0.97, top=0.93)

plt.figure(2)
plt.contourf(m3, a2, fraction) #[2,2.5,3,3.5,4,4.5])
plt.xlabel(r'$m_3 / \rm M_\odot$')
plt.ylabel(r'$a_2 / \rm au$')
plt.xscale('log')
plt.yscale('log')
plt.text(2.02, 1., 'fraction')        
plt.colorbar()
#         plt.subplots_adjust(left=0.18, bottom=0.15, right=0.97, top=0.93)

# %%
# Plot
plt.rc('text', usetex=False)
matplotlib.rcParams.update({'font.size': 20})
plt.contourf(m3, a2, eccs)
plt.xlabel(r'$m_3 / \rm M_\odot$')
plt.ylabel(r'$a_2 / \rm au$')
plt.xscale('log')
plt.yscale('log')
plt.text(2.03, 1., r'$e_{\rm max}$')
plt.colorbar()
#         plt.subplots_adjust(left=0.18, bottom=0.15, right=0.97, top=0.93)

plt.figure(2)
plt.contourf(m3, a2, fraction) #[2,2.5,3,3.5,4,4.5])
plt.xlabel(r'$m_3 / \rm M_\odot$')
plt.ylabel(r'$a_2 / \rm au$')
plt.xscale('log')
plt.yscale('log')
plt.text(2.02, 1., 'fraction')        
plt.colorbar()
#         plt.subplots_adjust(left=0.18, bottom=0.15, right=0.97, top=0.93)

# %%
