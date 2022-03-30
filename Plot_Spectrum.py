import numpy as np
import matplotlib.pyplot as plt
import csv
import os 

# Preparing the canevas (figure)
fig, ax = plt.subplots()

File = 'spectrum' #[]

Counter = 0
title = {'U99':1, 'U75':2, 'U67':3, 'U50':4, 'U33':5, 'U25':6}
# read the lines from the file
with open(f'{File}', "r") as f:
    lines = f.readlines()
    f.close()

En = []
Flux99 = []
Flux75 = []
Flux67 = []
Flux50 = []
Flux33 = []
Flux25 = []

with open(f'{File}', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=' ')
    next(plots, None) # skip the header
    for row in plots:
        En.append(float(row[0]))
        Flux99.append(1.802e+15*float(row[1]))
        Flux75.append(1.782e+15*float(row[2]))
        Flux67.append(1.775e+15*float(row[3]))
        Flux50.append(1.760e+15*float(row[4]))
        Flux33.append(1.745e+15*float(row[5]))
        Flux25.append(1.738e+15*float(row[6]))
    En = np.array(En)
    Flux99 = np.array(Flux99)
    Flux75 = np.array(Flux75)
    Flux67 = np.array(Flux67)
    Flux50 = np.array(Flux50)
    Flux33 = np.array(Flux33)
    Flux25 = np.array(Flux25)
    
# Creating the plot
# Tuning the canevas
ax.set(xlabel='Incident Neutron Energy [eV]', ylabel='Flux [n/cm²]', title='Neutron Spectrum')
ax.grid()
ax.loglog(En, Flux99, drawstyle='steps', label='UO2 = 99%')
ax.loglog(En, Flux75, drawstyle='steps', label='UO2 = 75%')
ax.loglog(En, Flux67, drawstyle='steps', label='UO2 = 67%')
ax.loglog(En, Flux50, drawstyle='steps', label='UO2 = 50%')
ax.loglog(En, Flux33, drawstyle='steps', label='UO2 = 33%')
ax.loglog(En, Flux25, drawstyle='steps', label='UO2 = 25%')
ax.legend(loc='lower center', handlelength=3, handleheight=3.5)

# showing and saving the plot
fig.savefig("test.png", DPI=1200)
plt.show()

