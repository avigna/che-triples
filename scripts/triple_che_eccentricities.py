#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# %%
"""
Created on Tue Mar 29 16:58:54 2022

@author: evgeni
"""
import pandas as pd
import numpy as np
import math
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from matplotlib import rc
from scipy.optimize import fsolve
import matplotlib.colors as colors
rc('text', usetex=False)
from astropy import constants as const
# %%
# %%
#define constants
# c=2.99792458e10
# msun=2e33
# G=6.674e-8
# au=1.5e13
# rsun = 6.9e10

c=const.c.cgs.value
msun=const.M_sun.cgs.value
G=const.G.cgs.value
au=const.au.cgs.value
rsun=const.R_sun.cgs.value

# define some useful quantities

#this is self explanatory - I hope
def period_ratio(m1,m2,m3,a1,a2):
    m_bin = m1 + m2
    m_tot = m_bin + m3
    return (a1/a2)**1.5 * (m_tot/m_bin)**0.5

# how bad the double averaging is (e.g. Grishin, Perets and Fragione 2018; Luo et al., 2016)
def eps_sa(m1,m2,m3,a1,a2,e2):
    m_bin = m1+m2
    m_tot = m_bin + m3
    return (a1/a2/(1-e2**2))**1.5*m3/(m_tot*m_bin)**0.5

# angular momentum ratio (e.g. Haim and Katz 2018)
def eta(m1,m2,m3,a1,a2,e2):
    m_bin = m1+m2
    m_tot = m_bin + m3
    mu_in = m1*m2/m_bin
    mu_out = m3*m_bin/m_tot
    return mu_in/mu_out*(m_bin*a1/m_tot/a2/(1-e2**2))**0.5

#importance of GR corrections (e.g. Liu, Munoz and Lai 2015)
def eps_gr(m1,m2,m3,a1,a2,e2):
    m_bin = m1+m2
    rg_over_a1= G*m_bin*msun/c**2/(a1*au)
 #   print (rg_over_a1)
    return 3 * m_bin*(1-e2**2)**1.5/m3*(a2/a1)**3 * rg_over_a1

#importance of tidal corrections (e.g. Liu, Munoz and Lai 2015)
def eps_tide(m1,m2,m3,a1,a2,e2, r1, r2, k1, k2):
    m_bin = m1+m2
    b2 = a2 * (1-e2**2)**0.5
    return 15 * m_bin * b2**3 * (m1**2 * k1 * r1**5 + m2**2 * k2 * r2**5 ) / a1**8 / m1 / m2/ m3

#octupole term (e.g. Naoz 2016 and many others)
def eps_oct(m1,m2,a1,a2,e2):
    return (m1-m2)/(m1+m2)* (a1/a2) * e2/(1-e2**2)

# %%
# calculate the maximal eccentricity according to Mangipudi et al., 2022 and the appendix with tides
def get_ecc(m1, m2, m3, a1, a2, e2, r1, r2, k1, k2, cos_inc):
    cosi_0 = cos_inc
    ETA = eta(m1,m2,m3,a1,a2,e2)
    eps_GR1 = eps_gr(m1, m2, m3, a1, a2, e2)
    eps_SA1 = eps_sa(m1, m2, m3, a1, a2, e2)
    eps_Tide1 = eps_tide(m1, m2, m3, a1, a2, e2, r1, r2, k1, k2)
    
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

print (get_ecc(50, 50, 200, 0.2, 1, 0, 7/214, 7/214, 0.001, 0,0 ))

# %%
import matplotlib
#plot e_max for a particular choice of inclination cosinc - set e_maxplot = True
#plot fractions of collisoins - uniform in cos inc with Ni inclinations
# set fraction_plot = True
#the outer grid of masses and separations is of size NN x NN

#the system is colliding if the pericentre r_p < 1.3 * (r1+r2)
# AVG: double check this previous condition
def plot_grids(NN, Ni,period,m1,m2,r1,r2,k1,k2,cosinc,col):
    # AVG: explain all things here
    # AVG: probably better a script to first calculate the grid and another to then plot
    e_maxplot = True
    fraction_plot = True

#     m1=70; m2=70; 
    # Change 86400 to a constant
    a1 = (G*(m1+m2)*msun * (period * 86400)**2/4/np.pi/np.pi)**0.3333/au
    r1 = r1 * rsun / au; r2 = r2 * rsun / au
    m3 = np.logspace(0,4,NN)
    a2 = np.logspace(np.log10(a1) + np.log10(2), np.log10(a1)+ np.log10(300), NN)
    eccs = np.zeros(shape=[NN,NN])
    rps = np.zeros(shape=[NN,NN])
    fraction = np.zeros(shape=[NN,NN])
    cos_incs=np.linspace(-1,1,Ni)

    for i in range(0,NN):
        print (i)
        for j in range(0,NN):
            eccs[i,j] = get_ecc(m1,m2,m3[i],a1,a2[j],0,r1 ,r2,k1,k2, cosinc) 
            rps[i,j] = a1 * (1 - eccs[i,j]) * au / rsun
            counter=0
            for k in range(0,Ni):
               ec = get_ecc(m1,m2,m3[i],a1,a2[j],0,r1 ,r2,k1,k2, cos_incs[k])
               rpp = a1 * (1 - ec) 

               if rpp/1.3/(r1+r2)<=1:
                   counter=counter+1

            fraction[i,j] = counter/Ni       
                

    plt.rc('text', usetex=False)
    matplotlib.rcParams.update({'font.size': 20})
    if e_maxplot:
        plt.contourf(np.log10(a2), np.log10(m3), eccs, levels=np.linspace(0,1,21))
        plt.text(1., 4.03, r'$e_{\rm max}$')
        plt.xlabel(r'$\log a_2 / \rm au$')
        plt.ylabel(r'$\log m_2 / \rm M_\odot$')
        plt.colorbar()
        plt.subplots_adjust(left=0.18, bottom=0.15, right=0.97, top=0.93)

    #   plt.contour((a2), np.log10(m3), rps/(r1+r2) * rsun/au, [1.3], colors=col) #[2,2.5,3,3.5,4,4.5])
    if fraction_plot:
        plt.figure(2)
        plt.contourf(np.log10(a2), np.log10(m3), fraction) #[2,2.5,3,3.5,4,4.5])
        plt.text(1., 4.02, 'fraction')
        plt.xlabel(r'$\log a_2 / \rm au$')
        plt.ylabel(r'$\log m_2 / \rm M_\odot$')
        plt.colorbar()
        plt.subplots_adjust(left=0.18, bottom=0.15, right=0.97, top=0.93)


    #  plt.colorbar(format="%.2f")
 #   plt.text(4,150 + 60*cosinc, str(cosinc), color=col)
    return eccs, rps, fraction

# #%%
# NN=10
# a1=0.16
# m3 = np.logspace(-0,4,NN)
# a2 = np.logspace(np.log10(a1) + np.log10(2), np.log10(a1)+ np.log10(300), NN)
# m3_h=[140/24*(x/a1)**3 for x in a2]
# m_mardling = [140*(1 + (x/2.8/a1)**2.5 ) for x in a2]
# def f_tory(q):
#     return 10**(-0.6+0.04*q) * q **(0.32+0.1*q)

# q_m = [140/x for x in m3]
# ff = [f_tory(qq) for qq in q_m]
# a2_tory = [max(-1, a1/f) for f in ff]
# #plt.plot([np.log10(x) for x in np.logspace(-6,0,100)], [(f_tory(x)) for x in np.logspace(-6,0,100)])
# #plt.ylim([0,3])
# #plt.xlim([-1,2])

# #%%

# plt.contour(np.log10(a2), np.log10(m3), (res[0]+0*1e-2),levels=np.linspace(0.01,0.81,9), cmap='RdBu') #[2,2.5,3,3.5,4,4.5])
# plt.colorbar()
# plt.plot(np.log10(a2),np.log10(m3_h), linewidth=3, color='gray')
# plt.plot(np.log10(a2_tory), np.log10(m3), linewidth=3)
# plt.plot(np.log10(a2),np.log10(m_mardling), linewidth=3, color='k')
# plt.axvline(np.log10(2.2*a1), color='red', linewidth=3)
# #plt.plot(np.log10(a2), np.log10(m3_h), color='k', linewidth='3')
# plt.ylim([0.,4])
# plt.xlim([-0.5,1.7])
# plt.xlabel(r'$\log_{10}\ a_2 / \rm au$')
# plt.ylabel(r'$log_{10}\ m_3 / \rm M_\odot$')
# plt.text(1.5,4.1, r'$f_{\rm contact}$')
# plt.subplots_adjust(left=0.15, bottom=0.16, right=0.97, top=0.92)

# %%
NN=40; Ni=30; r1=10; r2=10; m1=45; m2=45; k1=0.02; k2=0.02;
res = plot_grids(NN,Ni,2,m1,m2,r1,r2,k1,k2,0.5, 'red')

# %%
