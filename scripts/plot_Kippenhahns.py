#!/usr/bin/env python
# +
import mkipp
import kipp_data
import mesa_data
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch
import numpy as np

# import random
# import numpy as np
# import matplotlib
# import matplotlib.pyplot as plt
# from scipy.stats import kde
import pandas as pd
# import matplotlib.patches as mpatches
# from scipy.signal import argrelextrema

# +
# logs_DIR = ["../data/MESA/00_cat_55_55_1.1d_ZSMC/LOGS1/"]
# history_path = "../data/MESA/00_cat_55_55_1.1d_ZSMC/LOGS1/history.data"
# profile_path = "../data/MESA/00_cat_55_55_1.1d_ZSMC/LOGS1/"
# filename_mass = "../plots/pdf/55_Msun_high_Z/Kippenhahn_Mass_55_Msun_high_Z.pdf"
# filename_radius = "../plots/pdf/55_Msun_high_Z/Kippenhahn_Radius_55_Msun_high_Z.pdf"
# time_of_He_ZAMS = 5.31
# val_skip = 39

logs_DIR = ["../data/MESA/00_cat_55_55_1.1d_0.1ZSMC/LOGS1/"]
history_path = "../data/MESA/00_cat_55_55_1.1d_0.1ZSMC/LOGS1/history.data"
profile_path = "../data/MESA/00_cat_55_55_1.1d_0.1ZSMC/LOGS1/"
filename_mass = "../plots/pdf/55_Msun_low_Z/Kippenhahn_Mass_55_Msun_low_Z.pdf"
filename_radius = "../plots/pdf/55_Msun_low_Z/Kippenhahn_Radius_55_Msun_low_Z.pdf"
time_of_He_ZAMS = 5.12
val_skip = 58

df = pd.read_csv(history_path,  header=4, skiprows=0, delimiter='\s+',usecols=['age','apsidal_constant_k2'])

save_flag = True

# +
from matplotlib import rc
import matplotlib.ticker as ticker
# from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, \
#      grid, savefig, show
rc('text', usetex=True)
rc('font', family='serif')
rc('font', size=16)

#Reading out mixing regions and data, and plotting independently
kipp_args = mkipp.Kipp_Args(logs_dirs = logs_DIR, xaxis="star_age")
fig = plt.figure()
axis = plt.gca()
profile_paths = mesa_data.get_profile_paths(logs_DIR)
#if data is distributed among several history.data files, you can provide them
history_paths = [history_path]
#read profile data
#kipp_data.get_xyz_data returns an object containing
#   xyz_data.xlims : limits of data in x coordinate
#   xyz_data.X     : 2D array of xaxis values of profile data
#   xyz_data.Y     : 2D array of xaxis values of profile data
#   xyz_data.Z     : 2D array of xaxis values of profile data
# the last three can be used as inputs for matplotlib contour or contourf
xyz_data = kipp_data.get_xyz_data(profile_paths, kipp_args)
#read mixing regions 
#kipp_data.get_mixing_zones returns an object containing
#   mixing_zones.zones     : matplotlib Path objects for each mixing zone.
#                            These can be plotted using add_patch(...)
#   mixing_zones.mix_types : Integer array containing the type of each zone
#   mixing_zones.x_coords  : x coordinates for points at the surface
#   mixing_zones.y_coords  : y coordinates for points at the surface
#   mixing_zones.histories : mesa_data history files to access additional data
# the last three can be used as inputs for matplotlib contour or contourf
mixing_zones = kipp_data.get_mixing_zones(history_paths, kipp_args, xlims = xyz_data.xlims)
# just plot convection, overshooting and semiconvection
for i,zone in enumerate(mixing_zones.zones):
    color = ""
    #Convective mixing
    if mixing_zones.mix_types[i] == 1: #convection
        color = '#332288' # purple
    # # ???
    elif mixing_zones.mix_types[i] == 2: #???
        color = '#000000' #         
    #Overshooting 
    elif mixing_zones.mix_types[i] == 3: #overshooting
        color = '#117733' # green       
    #Semiconvective mixing
    elif mixing_zones.mix_types[i] == 4: #semiconvection
        color = '#CC6677'
#     # ???
#     elif mixing_zones.mix_types[i] == 5: #???
#         color = '#000000' #
#     # ???
#     elif mixing_zones.mix_types[i] == 6: #???
#         color = '#000000' #                    
    else:
        continue
    axis.add_patch(PathPatch(zone, color=color, alpha = 0.5, lw = 0))
# if (xyz_data.Z.size > 0):
#     CS = plt.contour(xyz_data.X, xyz_data.Y, xyz_data.Z, [0,4,8], colors='k')
#     plt.clabel(CS, inline=1, fontsize=10)
axis.plot(mixing_zones.x_coords[val_skip:-1], mixing_zones.y_coords[val_skip:-1], 'k', lw=4)
axis.set_xlabel("t/Myr")
axis.set_ylabel("$m/M_{\odot}$")
axis.set_xlim(0,max(mixing_zones.x_coords))
axis.set_ylim(0,max(mixing_zones.y_coords))

plt.axvline(x=time_of_He_ZAMS, color='b', linestyle='--')

# plt.axhline(y=50, color='b', linestyle='-')
# plt.axhline(y=46, color='b', linestyle='--')
# plt.axhline(y=41, color='b', linestyle=':')

ax2 = axis.twinx()  # instantiate a second Axes that shares the same x-axis
# color = '#9F4A54'
color = 'tab:red'
ax2.set_ylabel('k', color=color)  # we already handled the x-label with ax1
ax2.plot(df['age'][val_skip:-1]*10**-6, df['apsidal_constant_k2'][val_skip:-1], color=color, lw=2)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim([0, 0.025])

# ax2.set_yticks([0, 0.01, 0.02])
# yticks = [0, 0.01, 0.02]
# ylabels = [f'{x:1.2f}' for x in yticks]
# ax2.set_yticks(yticks, labels=ylabels)

if save_flag:
    plt.savefig(filename_mass, bbox_inches='tight')
    print(filename_mass)
# +
#Reading out mixing regions and data, and plotting independently
kipp_args = mkipp.Kipp_Args(logs_dirs = logs_DIR, xaxis="star_age", yaxis = "radius")
fig = plt.figure()
axis = plt.gca()
profile_paths = mesa_data.get_profile_paths(logs_DIR)
#if data is distributed among several history.data files, you can provide them
history_paths = [history_path]
#read profile data
#kipp_data.get_xyz_data returns an object containing
#   xyz_data.xlims : limits of data in x coordinate
#   xyz_data.X     : 2D array of xaxis values of profile data
#   xyz_data.Y     : 2D array of xaxis values of profile data
#   xyz_data.Z     : 2D array of xaxis values of profile data
# the last three can be used as inputs for matplotlib contour or contourf
xyz_data = kipp_data.get_xyz_data(profile_paths, kipp_args)
#read mixing regions 
#kipp_data.get_mixing_zones returns an object containing
#   mixing_zones.zones     : matplotlib Path objects for each mixing zone.
#                            These can be plotted using add_patch(...)
#   mixing_zones.mix_types : Integer array containing the type of each zone
#   mixing_zones.x_coords  : x coordinates for points at the surface
#   mixing_zones.y_coords  : y coordinates for points at the surface
#   mixing_zones.histories : mesa_data history files to access additional data
# the last three can be used as inputs for matplotlib contour or contourf
mixing_zones = kipp_data.get_mixing_zones(history_paths, kipp_args, xlims = xyz_data.xlims)
# just plot convection, overshooting and semiconvection
for i,zone in enumerate(mixing_zones.zones):
    color = ""
    #Convective mixing
    if mixing_zones.mix_types[i] == 1: #convection
        color = '#332288' # purple
    # # ???
    elif mixing_zones.mix_types[i] == 2: #???
        color = '#000000' #         
    #Overshooting 
    elif mixing_zones.mix_types[i] == 3: #overshooting
        color = '#117733' # green       
    #Semiconvective mixing
    elif mixing_zones.mix_types[i] == 4: #semiconvection
        color = '#CC6677'
#     # ???
#     elif mixing_zones.mix_types[i] == 5: #???
#         color = '#000000' #
#     # ???
#     elif mixing_zones.mix_types[i] == 6: #???
#         color = '#000000' #        
    else:
        continue
    axis.add_patch(PathPatch(zone, color=color, alpha = 0.5, lw = 0))
# if (xyz_data.Z.size > 0):
#     CS = plt.contour(xyz_data.X, xyz_data.Y, xyz_data.Z, [0,4,8], colors='k')
#     plt.clabel(CS, inline=1, fontsize=10)
axis.plot(mixing_zones.x_coords[val_skip:-1], mixing_zones.y_coords[val_skip:-1], 'k', lw=4)
axis.set_xlabel("t/Myr")
axis.set_ylabel("$r/R_{\odot}$")
axis.set_xlim(0,max(mixing_zones.x_coords))
axis.set_ylim(0,max(mixing_zones.y_coords))
# plt.axhline(y=3.2, color='r', linestyle='-')

plt.axvline(x=time_of_He_ZAMS, color='b', linestyle='--')

# plt.axhline(y=3.5, color='b', linestyle='-')
# plt.axhline(y=7.8, color='b', linestyle='--')
# plt.axhline(y=8.7, color='b', linestyle=':')
# plt.axhline(y=2.4, color='b', linestyle='-')

ax2 = axis.twinx()  # instantiate a second Axes that shares the same x-axis
# color = '#9F4A54'
# color = '#E8AA14'
color = 'tab:red'
ax2.set_ylabel('k', color=color)  # we already handled the x-label with ax1
ax2.plot(df['age'][val_skip:-1]*10**-6, df['apsidal_constant_k2'][val_skip:-1], color=color, lw=2)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim([0, 0.025])

if save_flag:
    plt.savefig(filename_radius, bbox_inches='tight')
    print(filename_radius)