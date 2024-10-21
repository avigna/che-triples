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
from scipy.io import savemat

import compact_triple_analysis as cta
# %%
# From MESA simulations
# Z_SMC
idx_Z_SMC = [39, 184, 1340, 2286]
age_yr_Z_SMC = [8.392999e+03, 2.697549e+06, 5.310639e+06, 5.683500e+06]
apsidal_constant_k2_Z_SMC = [0.018226, 0.011376, 0.006834, 0.001726]
period_days_Z_SMC = [1.100177, 1.087004, 2.220218, 4.525685]
radius_Rsun_Z_SMC = [9.038323, 10.171930, 2.016422, 0.448425]
mass_conv_core_Z_SMC = [37.751202, 38.110672, 19.155549, 0.000000]
mass_Msun_Z_SMC = [54.997851, 51.413071, 26.583868, 18.619564]

# 0.1 Z_SMC
idx_0_1_Z_SMC = [47, 1784, 2132, 2504]
age_yr_0_1_Z_SMC = [5.329669e+03, 5.120677e+06, 5.307573e+06, 5316369.558564772]
apsidal_constant_k2_0_1_Z_SMC = [0.023533, 0.008881, 0.004897, 0.0007136587012452427]
period_days_0_1_Z_SMC = [1.099654, 1.450562, 1.668519, 1.6824947859764587]
radius_Rsun_0_1_Z_SMC = [7.577242, 2.401224, 1.939817, 0.8534730543749688]
mass_conv_core_Msun_0_1_Z_SMC = [37.805430, 41.412994, 38.487511, 0.0]
mass_Msun_0_1_Z_SMC = [54.999836, 43.976249, 43.793192, 43.79319230382309]

# %%
# DEFINE SYSTEM
debug_flag = False
NN=250; 
Ni=181; # Make sure the list includes i=0

#-----------------------------------------------------------------------------------#
# Z=Z_SMC: 55+55 Msun circular binary with a 1.1 d orbital period (at ZAMS)
# # Evaluated at ZAMS
# chosen_metallicity = 0.0035;
# chosen_apsidal_constant_k2 = apsidal_constant_k2_Z_SMC[0];
# chosen_mass_Msun = mass_Msun_Z_SMC[0];
# chosen_period_days = period_days_Z_SMC[0];
# chosen_radius_Rsun = radius_Rsun_Z_SMC[0];
# flag_CHE = 1;

# # Evaluated at HeMS
# chosen_metallicity = 0.0035;
# chosen_apsidal_constant_k2 = apsidal_constant_k2_Z_SMC[2];
# chosen_mass_Msun = mass_Msun_Z_SMC[2];
# chosen_period_days = period_days_Z_SMC[2];
# chosen_radius_Rsun = radius_Rsun_Z_SMC[2];
# flag_CHE = 1;

# # Evaluated at the end of the simulation
# chosen_metallicity = 0.0035;
# chosen_apsidal_constant_k2 = apsidal_constant_k2_Z_SMC[2];
# chosen_mass_Msun = mass_Msun_Z_SMC[3];
# chosen_period_days = period_days_Z_SMC[3];
# chosen_radius_Rsun = radius_Rsun_Z_SMC[3];
# flag_CHE = 1;

# # Evaluated at the end of the simulation, assuming BBH formation
# chosen_metallicity = 0.0035;
# chosen_apsidal_constant_k2 = 0.0;
# chosen_mass_Msun = mass_Msun_Z_SMC[3];
# chosen_period_days = period_days_Z_SMC[3];
# chosen_radius_Rsun = radius_Rsun_Z_SMC[3];
# flag_CHE = 1;
#-----------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------#
# Z=0.1Z_SMC: 55+55 Msun circular binary with a 1.1 d orbital period (at ZAMS)
# Evaluated at ZAMS
chosen_metallicity = 0.00035;
chosen_apsidal_constant_k2 = apsidal_constant_k2_0_1_Z_SMC[0];
chosen_mass_Msun = mass_Msun_0_1_Z_SMC[0];
chosen_period_days = period_days_0_1_Z_SMC[0];
chosen_radius_Rsun = radius_Rsun_0_1_Z_SMC[0];
flag_CHE = 1;

# # Evaluated at HeMS
# chosen_metallicity = 0.00035;
# chosen_apsidal_constant_k2 = apsidal_constant_k2_0_1_Z_SMC[1];
# chosen_mass_Msun = mass_Msun_0_1_Z_SMC[1];
# chosen_period_days = period_days_0_1_Z_SMC[1];
# chosen_radius_Rsun = radius_Rsun_0_1_Z_SMC[1];
# flag_CHE = 1;

# # Evaluated at the end of the simulation
# chosen_metallicity = 0.00035;
# chosen_apsidal_constant_k2 = apsidal_constant_k2_0_1_Z_SMC[3];
# chosen_mass_Msun = mass_Msun_0_1_Z_SMC[3];
# chosen_period_days = period_days_0_1_Z_SMC[3];
# chosen_radius_Rsun = radius_Rsun_0_1_Z_SMC[3];
# flag_CHE = 1;

# # Evaluated at the end of the simulation, assuming BBH formation
# chosen_metallicity = 0.00035;
# chosen_apsidal_constant_k2 = 0;
# chosen_mass_Msun = mass_Msun_0_1_Z_SMC[3];
# chosen_period_days = period_days_0_1_Z_SMC[3];
# chosen_radius_Rsun = radius_Rsun_0_1_Z_SMC[3];
# flag_CHE = 1;

print("k_2 = ", chosen_apsidal_constant_k2)
print("M/Msun = ", chosen_mass_Msun)
print("Z = ", chosen_metallicity)
print("P_{orb}/d = ", chosen_period_days)
print("R/Rsun = ", chosen_radius_Rsun)
print("CHE: ", flag_CHE)

# %%
is_CHE = flag_CHE; metallicity = chosen_metallicity; period_days=chosen_period_days; m1=chosen_mass_Msun; m2=m1; r1=chosen_radius_Rsun; r2=r1; r1_core=chosen_radius_Rsun; r2_core=r1_core; k1=chosen_apsidal_constant_k2; k2=k1;

# %%
eps_SA_flag = False
eps_GR_flag = False
eps_Tides_flag = False

string_to_save = "triple_Z="+str(metallicity)+"_CHE="+str(is_CHE)+"_M1=M2="+str(m2)+"_Porb="+str(period_days)

if eps_SA_flag:
    string_to_save+="_SA"
    
if eps_GR_flag:
    string_to_save+="_GR"   
    
if eps_Tides_flag:
    string_to_save+="_Tides"   
    
string_to_save+=".mat"
print(string_to_save)

# CALCULATE
cos_inc_max_ecc, ecc_grid, eps_SA, eps_GR, eps_Tide, etas, f_merger, m3, p2, tau_sec = cta.calculate_e_grid(NN, Ni, m1, m2, period_days, r1, r2, r1_core, r2_core, k1, k2, eps_SA_flag, eps_GR_flag, eps_Tides_flag, debug_flag)

# SAVE
print(string_to_save)

mdic = {"cos_inc":cos_inc_max_ecc, "eccs":ecc_grid, "eps_SA":eps_SA, "eps_GR":eps_GR, "eps_Tide":eps_Tide, "eta":etas, "f_merger":f_merger, "m3":m3, "p2":p2, "tau_sec":tau_sec}
savemat(string_to_save, mdic)

# %%
eps_SA_flag = True
eps_GR_flag = False
eps_Tides_flag = False

string_to_save = "triple_Z="+str(metallicity)+"_CHE="+str(is_CHE)+"_M1=M2="+str(m2)+"_Porb="+str(period_days)

if eps_SA_flag:
    string_to_save+="_SA"
    
if eps_GR_flag:
    string_to_save+="_GR"   
    
if eps_Tides_flag:
    string_to_save+="_Tides"   
    
string_to_save+=".mat"
print(string_to_save)

# CALCULATE
cos_inc_max_ecc, ecc_grid, eps_SA, eps_GR, eps_Tide, etas, f_merger, m3, p2, tau_sec = cta.calculate_e_grid(NN, Ni, m1, m2, period_days, r1, r2, r1_core, r2_core, k1, k2, eps_SA_flag, eps_GR_flag, eps_Tides_flag, debug_flag)

# SAVE
print(string_to_save)

mdic = {"cos_inc":cos_inc_max_ecc, "eccs":ecc_grid, "eps_SA":eps_SA, "eps_GR":eps_GR, "eps_Tide":eps_Tide, "eta":etas, "f_merger":f_merger, "m3":m3, "p2":p2, "tau_sec":tau_sec}
savemat(string_to_save, mdic)

# %%
eps_SA_flag = False
eps_GR_flag = True
eps_Tides_flag = False

string_to_save = "triple_Z="+str(metallicity)+"_CHE="+str(is_CHE)+"_M1=M2="+str(m2)+"_Porb="+str(period_days)

if eps_SA_flag:
    string_to_save+="_SA"
    
if eps_GR_flag:
    string_to_save+="_GR"   
    
if eps_Tides_flag:
    string_to_save+="_Tides"   
    
string_to_save+=".mat"
print(string_to_save)

# CALCULATE
cos_inc_max_ecc, ecc_grid, eps_SA, eps_GR, eps_Tide, etas, f_merger, m3, p2, tau_sec = cta.calculate_e_grid(NN, Ni, m1, m2, period_days, r1, r2, r1_core, r2_core, k1, k2, eps_SA_flag, eps_GR_flag, eps_Tides_flag, debug_flag)

# SAVE
print(string_to_save)

mdic = {"cos_inc":cos_inc_max_ecc, "eccs":ecc_grid, "eps_SA":eps_SA, "eps_GR":eps_GR, "eps_Tide":eps_Tide, "eta":etas, "f_merger":f_merger, "m3":m3, "p2":p2, "tau_sec":tau_sec}
savemat(string_to_save, mdic)

# %%
eps_SA_flag = False
eps_GR_flag = False
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

# CALCULATE
cos_inc_max_ecc, ecc_grid, eps_SA, eps_GR, eps_Tide, etas, f_merger, m3, p2, tau_sec = cta.calculate_e_grid(NN, Ni, m1, m2, period_days, r1, r2, r1_core, r2_core, k1, k2, eps_SA_flag, eps_GR_flag, eps_Tides_flag, debug_flag)

# SAVE
print(string_to_save)

mdic = {"cos_inc":cos_inc_max_ecc, "eccs":ecc_grid, "eps_SA":eps_SA, "eps_GR":eps_GR, "eps_Tide":eps_Tide, "eta":etas, "f_merger":f_merger, "m3":m3, "p2":p2, "tau_sec":tau_sec}
savemat(string_to_save, mdic)

# %%
eps_SA_flag = True
eps_GR_flag = True
eps_Tides_flag = False

string_to_save = "triple_Z="+str(metallicity)+"_CHE="+str(is_CHE)+"_M1=M2="+str(m2)+"_Porb="+str(period_days)

if eps_SA_flag:
    string_to_save+="_SA"
    
if eps_GR_flag:
    string_to_save+="_GR"   
    
if eps_Tides_flag:
    string_to_save+="_Tides"   
    
string_to_save+=".mat"
print(string_to_save)

# CALCULATE
cos_inc_max_ecc, ecc_grid, eps_SA, eps_GR, eps_Tide, etas, f_merger, m3, p2, tau_sec = cta.calculate_e_grid(NN, Ni, m1, m2, period_days, r1, r2, r1_core, r2_core, k1, k2, eps_SA_flag, eps_GR_flag, eps_Tides_flag, debug_flag)

# SAVE
print(string_to_save)

mdic = {"cos_inc":cos_inc_max_ecc, "eccs":ecc_grid, "eps_SA":eps_SA, "eps_GR":eps_GR, "eps_Tide":eps_Tide, "eta":etas, "f_merger":f_merger, "m3":m3, "p2":p2, "tau_sec":tau_sec}
savemat(string_to_save, mdic)

# %%
eps_SA_flag = True
eps_GR_flag = False
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

# CALCULATE
cos_inc_max_ecc, ecc_grid, eps_SA, eps_GR, eps_Tide, etas, f_merger, m3, p2, tau_sec = cta.calculate_e_grid(NN, Ni, m1, m2, period_days, r1, r2, r1_core, r2_core, k1, k2, eps_SA_flag, eps_GR_flag, eps_Tides_flag, debug_flag)

# SAVE
print(string_to_save)

mdic = {"cos_inc":cos_inc_max_ecc, "eccs":ecc_grid, "eps_SA":eps_SA, "eps_GR":eps_GR, "eps_Tide":eps_Tide, "eta":etas, "f_merger":f_merger, "m3":m3, "p2":p2, "tau_sec":tau_sec}
savemat(string_to_save, mdic)

# %%
eps_SA_flag = False
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

# CALCULATE
cos_inc_max_ecc, ecc_grid, eps_SA, eps_GR, eps_Tide, etas, f_merger, m3, p2, tau_sec = cta.calculate_e_grid(NN, Ni, m1, m2, period_days, r1, r2, r1_core, r2_core, k1, k2, eps_SA_flag, eps_GR_flag, eps_Tides_flag, debug_flag)

# SAVE
print(string_to_save)

mdic = {"cos_inc":cos_inc_max_ecc, "eccs":ecc_grid, "eps_SA":eps_SA, "eps_GR":eps_GR, "eps_Tide":eps_Tide, "eta":etas, "f_merger":f_merger, "m3":m3, "p2":p2, "tau_sec":tau_sec}
savemat(string_to_save, mdic)

# %%
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

# CALCULATE
cos_inc_max_ecc, ecc_grid, eps_SA, eps_GR, eps_Tide, etas, f_merger, m3, p2, tau_sec = cta.calculate_e_grid(NN, Ni, m1, m2, period_days, r1, r2, r1_core, r2_core, k1, k2, eps_SA_flag, eps_GR_flag, eps_Tides_flag, debug_flag)

# SAVE
print(string_to_save)

mdic = {"cos_inc":cos_inc_max_ecc, "eccs":ecc_grid, "eps_SA":eps_SA, "eps_GR":eps_GR, "eps_Tide":eps_Tide, "eta":etas, "f_merger":f_merger, "m3":m3, "p2":p2, "tau_sec":tau_sec}
savemat(string_to_save, mdic)

# %%
if debug_flag:
    cta.test_eccentricity_with_gr(1000)
