#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# %%
"""
Created on Tue Mar 29 16:58:54 2022

@author: evgeni grishin
"""
import numpy as np
from matplotlib import pyplot as plt

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

# importance of GR corrections (e.g. Liu, Muñoz and Lai 2015)
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

def fsolve_bisection(func, a, b, args=(), tolerance=1e-6, max_iterations=100):
    cosi_0, ETA, eps_SA1, eps_GR1, eps_Tide1 = args
    """
    Bisection method for finding a root of a real-valued function.

    Parameters:
        func (function): The function for which to find the root.
        a (float): The lower bound of the interval.
        b (float): The upper bound of the interval.
        args (tuple): Extra arguments to pass to the function.
        tolerance (float): The desired accuracy of the solution.
        max_iterations (int): Maximum number of iterations.

    Returns:
        float: Approximation of the root.
    """
    if func(a, *args) * func(b, *args) >= 0:
        print ('no solution for eps_gr = ', eps_GR1, 'eps_tide = ', eps_Tide1, 'cos_0 = ', cosi_0)
       # raise ValueError("Function does not change sign over the interval [a, b].")
        return 0

    iteration = 0
    while (b - a) / 2 > tolerance and iteration < max_iterations:
        c = (a + b) / 2
        if func(c, *args) == 0:
            break
        elif func(c, *args) * func(a, *args) < 0:
            b = c
        else:
            a = c
        iteration += 1

    return (a + b) / 2

#Evgeni: 06.05.2024
# function to solve for. The first three terms are from Mingipudi. The last term is new and based on Liu, Muñoz and Lai 2015
# Hasn't been tested with N-body where eps_tide and eps_sa are relevant.
def F(e_m , cosi_0, ETA, eps_SA1, eps_GR1, eps_Tide1):
    zlk_term = (1 - e_m**2)*(-3 + ETA*cosi_0 + e_m**2*ETA**2/4) + 5*(cosi_0 + ETA*e_m**2/2)**2
    sa_term = -eps_SA1 * (3/8) *(-((1 - e_m**2)*cosi_0**3/e_m**2) - (ETA/2)*(1 - e_m**2) + (cosi_0 + e_m**2*ETA/2)**3*(1/e_m**2-16)- 9*(1 - e_m**2)*(cosi_0 + e_m**2*ETA/2)) 
    gr_term = - (eps_GR1*8/3)*((1 - e_m**2)/(e_m**2+1e-8))*(1 - (1 - e_m**2)**(-1/2))
    tide_term = - eps_Tide1*(8/45)*(1 - (1 + 3*e_m**2 + (3/8)*e_m**4)/( (1 - e_m**2)**(9/2) )  )
        
    return  zlk_term + sa_term + gr_term + tide_term

# %%
# calculate the maximal eccentricity according to Mangipudi et al. 2022 and the appendix with tides
#28.03.2024 - Evgeni: added a brute force calculation of a grid search. It runs slower and less accurate
#However it can be more stable, since fsolve struggles with the large dynamical range
# The initial grid choice also depends on epg_GR
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
        
    #add the dashed line from eq 23 of Anderson, Storch and Lai 2017:
    eps_srf = 4 * (eps_GR1 + eps_Tide1) / 3 #eq 18 with epr_rot= 0
    cos_I_crit = - 2 / (ETA + 1e-12) * (1 + eps_srf / 2)
    if cosi_0 <= cos_I_crit:
        return 0      
       
    # in the next function, what does F means? 
    args = (cosi_0, ETA, eps_SA1, eps_GR1, eps_Tide1)
    e_max = fsolve_bisection(F, a=1e-8, b=1-1e-8, args = args)
    
    delta_e = 0
    # the fluctuating terms in eq 24 of Mangipudi et al., 2022, calculated only for eps_sa > 0.01
    if eps_SA1 > 0.01:
        j_min = (1 - e_max**2)**0.5
        K2 = j_min * cosi_0 + j_min**2 * ETA / 2
        C = 9 / 8 * (1 + np.fabs(K2)* ETA) / (0.6**0.5 + j_min * ETA)
        delta_e = C * eps_SA1 * e_max * (j_min - e_max**2 * C * eps_SA1 / 2 )

    e_tot = e_max + delta_e
    if e_tot < 0:
        e_tot = 0
    if e_tot >= 1:
        e_tot = 1
    
    return e_tot

# %%
# TEST
# Alejandro: are we sending this to a test file/function?
# def test_eccentricity_with_gr(N):
#     m1=m2=1
#     m3=10
#     a1=1
#     a2=10
#     e2=0
#     eps_SA1=0.
#     eps_Tide1=0
#     eps_GR1=np.logspace(-1.5, 1.5,N)#np.logspace(-2,1,201)
#     cos_inc=0
#     debug = False
#     em=[]
#     em_analytic2 = [(1 - 64 / 81 * x**2) for x in eps_GR1]
#     for i in range(0,len(eps_GR1)):
#         em_to_append = get_maximal_eccentricity(m1, m2, m3, a1, a2, e2, eps_SA1, eps_GR1[i], eps_Tide1, cos_inc, debug)
#         em.append(em_to_append)#%%
#     plt.figure(figsize=(8,6))
#     plt.plot(eps_GR1, em, linewidth=3, label='numerical')
#     plt.plot(eps_GR1, [max(1e-20, x)**0.5 for x in em_analytic2],linewidth=3, label=r'$\epsilon_{\rm GR} \ll 1\ \rm limit$')
#     plt.xlabel(r'$\log \epsilon_{\rm GR}$')
#     plt.ylabel(r'$\log e_{\rm max}$')
#     plt.ylim(2e-5,2)
#     plt.legend()

#     plt.yscale('log')
#     plt.xscale('log')
#     plt.subplots_adjust(left=0.15, bottom=0.15, right=0.98, top=0.92,)

# test_eccentricity_with_gr(1000)

# %%
def calculate_e_grid(N,Ni, m1, m2, p, r1, r2, rc1, rc2, k1, k2, eps_SA_flag, eps_GR_flag, eps_Tides_flag, debug):
    a1_cubed = G * msun * (m1 + m2) * (p * day_to_s)**2 / 4 / np.pi**2
    a1_in_au = a1_cubed**(1/3) / au
    r1 = r1 * rsun / au
    r2 = r2 * rsun / au

    m3 = np.logspace(0,2,N)
    a2 = np.logspace(np.log10(2 * a1_in_au), np.log10(100 * a1_in_au), N)

    e2 = 0
    eps_SA1 = 0
    eps_GR1 = 0
    eps_Tide1 = 0       
    eps_SA = np.zeros([N,N])
    eps_GR = np.zeros([N,N])
    eps_Tide = np.zeros([N,N])      
    ecc_grid = np.zeros([N,N])
    etas = np.zeros([N,N])
    f_merger = np.zeros([N,N])
    tau_sec = np.zeros([N,N])
    cos_inc_max_ecc = np.zeros(shape=[N,N])
        
    if Ni==-1:
        # Here we should probably make it an user-defined choice which angle to pick
        cos_inc = 0
    else:
        cos_incs=np.linspace(-1,1,Ni)
        print(cos_incs)
        ec_max_temp = 0.0
        cos_inc_temp = []
    
    for i in range(0,len(m3)):
        for j in range (0, len(a2)):
            
            etas[i][j] = eta(m1, m2, m3[i], a1_in_au, a2[j], e2)
            tau_sec[i][j] = t_sec(m1, m2, m3[i], a1_in_au, a2[j], e2)
            
            if eps_SA_flag:
                eps_SA1 = epsilon_sa(m1, m2, m3[i], a1_in_au, a2[j], e2)
                eps_SA[i][j] = eps_SA1
            
            if eps_GR_flag:
                eps_GR1 = epsilon_gr(m1, m2, m3[i], a1_in_au, a2[j], e2)
                eps_GR[i][j] = eps_GR1                
            
            if eps_Tides_flag:
                eps_Tide1 = epsilon_tide(m1, m2, m3[i], a1_in_au, a2[j], e2, r1, r2, k1, k2)
                eps_Tide[i][j] = eps_Tide1        
        
            if Ni==-1:
                ecc_grid[i][j] = get_maximal_eccentricity(m1, m2, m3[i], a1_in_au, a2[j], e2, eps_SA1, eps_GR1, eps_Tide1, cos_inc, debug)
            else:
                ec = 0.0
                N_merger = 0
                for val in cos_incs:
                    e_temp = get_maximal_eccentricity(m1, m2, m3[i], a1_in_au, a2[j], e2, eps_SA1, eps_GR1, eps_Tide1, val, debug)

                    # Assess if a merger has occured following Marchant et al. (2016), which states that 
                    # the condition for merger is: R >= 1.32*f(1)*a*(1-e) \approx 0.5*a*(1-e) for q=1 
                    if r1>=(0.5*a1_in_au*(1-e_temp)):
                        N_merger = N_merger+1                   
                    
                    # Here we check what is the maximum eccentricity in the cos(i) parameter space
                    if e_temp>ec:
                        ec = e_temp
                        cos_inc_temp = val
                    
                ecc_grid[i][j] = ec
                cos_inc_max_ecc[i][j] = cos_inc_temp
                print(N_merger)
                f_merger[i][j] = N_merger/Ni
            
    p2 = [2*np.pi* ((x*au)**3 / G / msun / (m1+m2))**0.5/day_to_s  for x in a2]

    return cos_inc_max_ecc, ecc_grid, eps_SA, eps_GR, eps_Tide, etas, f_merger, m3, p2, tau_sec


# %%
# CALCULATE
NN=100; 
Ni=51; # Make sure the list includes i=0

# CHE
is_CHE = 1; metallicity = 0.0001; period_days=1; m1=55; m2=55; r1=7; r2=7; r1_core=3; r2_core=3; k1=0.024; k2=0.024;

# # TIC 470710327
# is_CHE = 0; metallicity = 0.142; period_days=1.1; m1=6; m2=6; r1=2.8; r2=2.8; r1_core=2.8; r2_core=2.8; k1=0.014; k2=0.014;

debug_flag = False
eps_SA_flag = True
eps_GR_flag = True
eps_Tides_flag = True

string_to_save = "triple_Z="+str(metallicity)+"_CHE="+str(is_CHE)+"_M1=M2="+str(m2)+"_Porb="+str(period_days)

if eps_SA_flag:
    string_to_save+="_SA"
    
if eps_GR_flag:
    string_to_save+="_GR"   
    
if eps_Tides_flag:
    string_to_save+="_Tides"   
    
string_to_save+=".mat"
print(string_to_save)

cos_inc_max_ecc, ecc_grid, eps_SA, eps_GR, eps_Tide, etas, f_merger, m3, p2, tau_sec = calculate_e_grid(NN, Ni, m1, m2, period_days, r1, r2, r1_core, r2_core, k1, k2, eps_SA_flag, eps_GR_flag, eps_Tides_flag, debug_flag)

# %%
# SAVE
from scipy.io import savemat

print(string_to_save)

mdic = {"cos_inc":cos_inc_max_ecc, "eccs":ecc_grid, "eps_SA":eps_SA, "eps_GR":eps_GR, "eps_Tide":eps_Tide, "eta":etas, "f_merger":f_merger, "m3":m3, "p2":p2, "tau_sec":tau_sec}
savemat(string_to_save, mdic)

# %%
plot_flag = True

# PLOT
if plot_flag:
    fig, ax = plt.subplots(figsize=(12,4.5))
    plt.subplot(121)
    CS1 = plt.contourf(m3, p2, np.transpose(ecc_grid), levels = np.linspace(0,1,21), cmap = 'viridis')
    CS2 = plt.contour(m3, p2, np.transpose(ecc_grid),  [0.1, 0.3, 0.5, 0.7, 0.9], cmap='Greys', linestyles='dashed', linewidths=3)

    plt.xlabel(r'$m_3\ [M_\odot]$')
    plt.ylabel(r'$P_{\rm out}\ \rm [day]$')
    plt.axvline(50, color='red', linestyle='dashed', linewidth=3)

    plt.xscale('log')
    plt.yscale('log')
    clb = plt.colorbar(CS1)
    clb.set_label(r'$e_{\rm max}$', labelpad=-150, y=1.12, rotation=0, size = 26)
    ax.clabel(CS2, inline=True, fontsize=10)

    plt.subplot(122)
    CS1 = plt.contourf(m3, p2, np.transpose(np.log10(etas)), cmap = 'viridis', levels=np.linspace(-2.5,1.5,41))
    CS2 = plt.contour(m3, p2, np.transpose(np.log10(etas)),  [0], linestyles='dashed', linewidths=3)
    plt.xlabel(r'$m_3\ [M_\odot]$')
    plt.ylabel(r'$P_{\rm out}\ \rm [day]$')

    plt.xscale('log')
    plt.yscale('log')
    clb = plt.colorbar(CS1)
    clb.set_label(r'$\log_{\rm 10} \eta$', labelpad=-150, y=1.12, rotation=0, size = 26)
    ax.clabel(CS2, inline=True, fontsize=10)

    plt.subplots_adjust(top = 0.88, bottom=0.16, left=0.1, right= 0.97, hspace=0.25)


# %%
plt.scatter(cos_inc_max_ecc,etas)
plt.xlabel(r'$cos(i)|_{e_{max}}$')
plt.ylabel(r'$\eta$')

# %%
f_merger

# %%

# %%
