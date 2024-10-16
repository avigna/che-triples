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
debug_flag = False
NN=200; 
Ni=51; # Make sure the list includes i=0

#-----------------------------------------------------------------------------------#
# Z=Z_SMC: 55+55 Msun circular binary with a 1.1 d orbital period (at ZAMS)
# Evaluated at the end of the simulation, assuming BBH formation
chosen_metallicity = 0.0035;
chosen_apsidal_constant_k2 = 0.0;
chosen_mass_Msun = mass_Msun_Z_SMC[3];
chosen_period_days = period_days_Z_SMC[3];
chosen_radius_Rsun = radius_Rsun_Z_SMC[3];
flag_CHE = 1;
#-----------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------#
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
