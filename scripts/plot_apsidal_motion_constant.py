import random
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import kde
import pandas as pd
import matplotlib.patches as mpatches
from scipy.signal import argrelextrema

debug_flag = True

# +
# history_path = "../data/MESA/00_cat_55_55_1.1d_ZSMC/LOGS1/history.data"
# filename = "apsidal_plot_55_Msun_Z_SMC.pdf"
# filename_full = "apsidal_plot_55_Msun_Z_SMC_full.pdf"

history_path = "../data/MESA/00_cat_55_55_1.1d_0.1ZSMC/LOGS1/history.data"
filename = "apsidal_plot_55_Msun_0_1_Z_SMC.pdf"
filename_full = "apsidal_plot_55_Msun_0_1_Z_SMC_full.pdf"

df = pd.read_csv(history_path,  header=4, skiprows=0, delimiter='\s+',usecols=['age','apsidal_constant_k2','center_h1','center_he4','star_mass','log_R','mass_conv_core','log_Teff','log_L','period_days'])
# -


if debug_flag:
    print(df)

age = df['age']
orbital_period = df['period_days']
center_h1 = df['center_h1']
center_he4 = df['center_he4']
mass = df['star_mass']
mass_conv_core = df['mass_conv_core']
log_Teff = df['log_Teff']
log_L = df['log_L']
log_R = df['log_R']
period_days = df['period_days']
apsidal_constant_k2 = df['apsidal_constant_k2']


# +
n = 5  # number of points to be checked before and after
# Find local peaks
retrieve_min = log_R.iloc[argrelextrema(log_R.values, np.less_equal, order=n)[0]]

idx=np.argmax(apsidal_constant_k2)
# -

# ZAMS_L2_overflow
# L2_overflow
# off_CHE
# PISN
# PPISN
# dir_coll

# +
from matplotlib import rc
from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, \
     grid, savefig, show
rc('text', usetex=True)
rc('font', family='serif')
rc('font', size=16)

plt.figure(figsize=(6,4))
plt.plot(age*10**-6,apsidal_constant_k2,label="k$_2$")
plt.xlabel("t/Myr",fontsize=16)
plt.ylabel("k$_2$",fontsize=16)
plt.xlim([0, max(age*10**-6)])
# plt.legend()
plt.savefig(filename_full, bbox_inches='tight')
plt.show()
# -

plt.plot(age,10**log_R)
plt.scatter(age[retrieve_min.index], 10**log_R[retrieve_min.index], c='r')

if debug_flag:
    print(max(apsidal_constant_k2))
    print(max(mass_conv_core))

# +

print(age[idx])
    print(apsidal_constant_k2[idx])
    print(period_days[idx])
    print(10**log_R[idx])
    print(mass_conv_core[idx])
# -

print(age[retrieve_min.index])
print(apsidal_constant_k2[retrieve_min.index])
print(period_days[retrieve_min.index])
print(10**log_R[retrieve_min.index])
print(mass_conv_core[retrieve_min.index])


