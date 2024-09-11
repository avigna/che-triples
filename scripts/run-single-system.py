#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# %%
import numpy as np
from matplotlib import pyplot as plt
from astropy import constants as const
from scipy.io import savemat

import compact_triple_analysis as cta
# %%
# DEFINE SYSTEM
NN=100; 
Ni=51; # Make sure the list includes i=0

# # CHE
# is_CHE = 1; metallicity = 0.0001; period_days=1; m1=45; m2=m1; r1=7; r2=r1; r1_core=3; r2_core=r1_core; k1=0.022; k2=k1;

# # CHE
# is_CHE = 1; metallicity = 0.0001; period_days=1; m1=55; m2=m1; r1=; r2=r1; r1_core=; r2_core=r1_core; k1=; k2=k1;

# TIC 470710327
is_CHE = 0; metallicity = 0.142; period_days=1.1; m1=6; m2=m1; r1=2.8; r2=r1; r1_core=2.8; r2_core=r1_core; k1=0.014; k2=k1;

debug_flag = False

# %%
if debug_flag:
    cta.test_eccentricity_with_gr(1000)

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
