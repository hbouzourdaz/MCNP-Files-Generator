# Implementation of matplotlib function
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib import ticker, cm
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import csv
import re
from itertools import islice

x = [] #X meshgrid
y = [] #Y meshgrid
z = [] #Z meshgrid cell level

# Reading the XY mesh values
file = 'meshtal'
filename = f'{file}'
with open(f'{filename}', 'r') as csvfile:
    for line in islice(csvfile, 13, None):
        plots = csv.reader(csvfile, skipinitialspace=True, delimiter=' ')
    #    print(line)
        for row in plots:
            x.append(float(row[0]))
            y.append(float(row[1]))
        x = np.array(x)
        y = np.array(y)
        print(x, y)
    triang = tri.Triangulation(x, y)

# # reading the filenames
# File_Flux = [] # meshtally flux file data 
# with open('Flux', 'r') as listefile :
#     lines = listefile.readlines()
#     N = len(lines)
#     print(N)
#     for line in lines:
#        file = line.rstrip() #f'{line}'
#        File_Flux.append(file)
#        print(file)
# print('N = ', N)

# print('File_Flux list name : ', File_Flux)




# preparing the canevas
fig, axs = plt.subplots(nrows=1, ncols=1) #, sharex=True, sharey=True) #, figsize = ())

# Z_list = ['Z0', 'Z1', 'Z2', 'Z3', 'Z4', 'Z5']
# print(Z_list[0])
# for file in File_Flux :
#    M = File_Flux.index(file)
#    print('M =', M)

# Getting the Z-value results from the data file
with open(f'{filename}', 'r') as results:
    for line in islice(results, 13, None):
        results = csv.reader(results, skipinitialspace=True, delimiter=' ')
        Res = []
        for row in results:
            Res.append(float(row[3]))
        Z = np.array(Res)
        Rmax = max(Z)
        Rmin = min(Z) #min(Z)
        print(Z)
        # print(file,'=', Z[10:15])
        #print('X =', x[10:15])
        #print('Y =', y[10:15])
    tcf = axs.tricontourf(triang, Z, levels=np.linspace(Rmin, Rmax, 50), cmap=cm.jet)
    axs.set_xlabel(r'X (cm) ')
    axs.set_title(f'{filename}')
    circle1 = plt.Circle((0, 0), 75, linewidth=1.5, linestyle='--', color='r', fill=False)
    circle2 = plt.Circle((0, 0), 85, linewidth=1.5, linestyle='--', color='r', fill=False)
    circle3 = plt.Circle((0, 0), 122, linewidth=1.5, linestyle='--', color='r', fill=False)
    circle4 = plt.Circle((0, 0), 137, linewidth=1.5, linestyle='--', color='r', fill=False)
    axs.add_artist(circle1)
    axs.add_artist(circle2)
    axs.add_artist(circle3)
    axs.add_artist(circle4)
    #axs.axvline(x=-20, ymin=0.39, linewidth=1.5, linestyle='-.', color='r')
    #axs.axvline(x=20, ymin=0.39, linewidth=1.5, linestyle='-.', color='r')
#shared Y axis
axs.set_ylabel(r'Y (cm) ')
axs.set_aspect('equal')

divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "05%", pad="3%")
plt.colorbar(tcf, cax=cax)
#cbar = fig.colorbar(tcf)
# plt.tight_layout() #pad=0.1, w_pad=0.5, h_pad=1.50)
plt.show()
#plt.savefig('figure_1.jpg', dpi=300)
