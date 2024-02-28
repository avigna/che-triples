import random
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import kde
import pandas as pd
import matplotlib.patches as mpatches

filename = '../data/MESA/45Msun/history.data' 
df = pd.read_csv(filename,  header=4, skiprows=0, delimiter='\s+',usecols=['age','apsidal_constant_k2','period_days'])


print(df)

age = df['age']
orbital_period = df['period_days']
#P = df['period_days']
#Rl_vs_R = df['rl_relative_overflow_1']
#Rl1 = df['rl_1']
#m1 = df['star_1_mass']
#Jdot = df['Jdot']
#jdot_gr = df['jdot_gr']
apsidal_constant_k2 = df['apsidal_constant_k2']


orbital_period

colors = []

# ZAMS_L2_overflow
# L2_overflow
# off_CHE
# PISN
# PPISN
# dir_coll

#for i in range(len(z)):
 #   colors.append("yellow")

#for i in range(len(z)):
 #   if('L2_overflow' in z[i]):
 #       colors[i]='tomato'
 #   if('ZAMS_L2_overflow' in z[i]):
 #       colors[i]='lightsalmon'
 #   if('dir_coll' in z[i]):
 #       colors[i]='black'
 #   if('off_CHE' in z[i]):
 #       colors[i]='powderblue'
 #   if('PISN' in z[i]):
 #       colors[i]='white'
 #   if('PPISN' in z[i]):
 #       colors[i]='slategrey'

# L2_overflow = mpatches.Patch(color='tomato', label='L2 overflow')
# ZAMS_L2_overflow = mpatches.Patch(color='lightsalmon', label='ZAMS L2 overflow')#dir_coll = mpatches.Patch(color='black', label='Direct collapse')
# PISN = mpatches.Patch(color='white', label='PISN')
# PPISN = mpatches.Patch(color='slategrey', label='PPISN')
# off_CHE = mpatches.Patch(color='powderblue', label='Evol. off CHE')


#plt.plot(age,jdot_gr,label="GR")
plt.plot(age,apsidal_constant_k2,label="k$_2$")


#plt.legend(handles=[L2_overflow,ZAMS_L2_overflow,off_CHE,PISN,PPISN,dir_coll],loc='upper
#left',ncol=2)
#plt.text(2.3,0.5,"log$_{10}$(Z/Z$_{\odot}$)=-2.625",fontsize=10)
plt.title("45 Msun")
plt.xlabel("star age",fontsize=12)
plt.ylabel("k$_2$",fontsize=12)
plt.legend()
plt.savefig("apsidal_plot.png")
plt.show()

