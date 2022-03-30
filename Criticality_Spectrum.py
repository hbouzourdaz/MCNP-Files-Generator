  # .................Criticité_Crt ...........................................
# .................. UFAS1 ...................................................
# .................. Python based script .....................................
#................... Author: S.E BENTRIDI, H. BOUZOURDAZ, H. MAKHLOUFI........
# .................. Starting Date: 08/09/2020 ...............................
# .................. present release : 0.3 ...................................
# .................. Last modif. date: 05/06/2021 ............................
import math
import os
import numpy as np
from decimal import *
import re

getcontext().prec = 12

# Initial parameters input
enr_U5 = float(input('Give the value of U5 enrichment in [%]:..... '))
f_UO2 = float(input('Give the value of UO2 volume fraction [%] in fuel pellet: ...... '))
H_core = float(input('Give the active hight of the core [cm]:..... '))


# Constants
NA = 6.022140E23
eV = 1.602176E-19
barn = 1.0E-24
# # Natural U235 enrichment is given by:
# Age = 2.0e+9
# enr_U5 = 100 / (1 + 137.84261 * math.exp(-8.29517579948633e-10 * Age))
# Atomic masses
# Fuel composition Molar masses
M_U5 = 235.04000
M_U8 = 238.05000
M_Pu9 = 239.05216
M_U = 0.01 * (M_U5 * enr_U5 + M_U8 * (100-enr_U5))
M_Th = 232.03806
# Water and Organic Matter Molar masses
M_H1 = 1.00790
M_O16 = 15.99900
M_C12 = 12.01074
M_N14 = 14.00670
# Uraninite and Thorite Molar Masses:
M_UO2 = (M_O16*2)+0.01*(M_U5*enr_U5+M_U8*(100-enr_U5))
M_ThO2 = (M_O16*2)+M_Th
# basical Compound densities [g/cc]
Rho_UO2 = 10.60000  # Uraninite [g/cc]
Rho_ThO2 = 9.86  # Thorite [g/cc]
Rho_Wtr0 = 1.00  # Water at normal P,T conditions [g/cc]
Rho_SiO2 = 2.65000  # Silica or Quartz [g/cc]
# 1MW power conversion constant : Physical Flux to MCNP Flux 
C = 7.62E16
# The Fuel definition
def Fuel(enr_U5, f_UO2):
    global d_UO2, d_ThO2, Rho_Fuel, M_UO2
    M_UO2 = (M_O16*2)+0.01*(M_U5*enr_U5+M_U8*(100-enr_U5))
    d_UO2 = f_UO2 * Rho_UO2 * 0.01
    d_ThO2 = (100-f_UO2) * Rho_ThO2 * 0.01
    Rho_Fuel = d_UO2 + d_ThO2
    print(f"""
--------------------------------------------------------------------------
The Fuel mass density is {Rho_Fuel:.7f} [g/cc]")
with:
Uraninite density fraction of {d_UO2:.7f} [g/CC]
U5 enriched with {enr_U5:.3f} [%]
Thorite density fraction of {d_ThO2:.7f} [g/cc]
--------------------------------------------------------------------------
""")
    return Rho_Fuel

# The isotopes fraction definition
def ZAID_Fuel(enr_U5, f_UO2):
    global d_UO2, d_ThO2, Rho_Fuel, M_UO2
    global d_U5, d_U8, d_Th, d_O16_Fuel
    # Matrix with ZAID generation
    d_U = d_UO2 * (M_U / M_UO2)
    d_U5 = 0.01 * (d_U * enr_U5)
    d_U8 = 0.01 * (d_U * (100 - enr_U5))
    d_Th = d_ThO2 * (M_Th / M_ThO2)
    # Oxygen fraction calculation
    d_O16_UO2 = (d_UO2 * M_O16 * 2) / (M_UO2)
    d_O16_ThO2 = (d_ThO2 * M_O16 * 2) / (M_ThO2)
    d_O16_Fuel = (d_O16_UO2 + d_O16_ThO2)
    return d_U8, d_U5, d_Th, d_O16_Fuel

# Header Definition
def header(enr_U5, f_UO2, H_core):
    global filename, Output
    Rho_Fuel=Fuel(enr_U5, f_UO2)
    # Files naming protocol
    filename = f"U{int(f_UO2)}e{int(enr_U5)}"
    Output = f"U{int(f_UO2)}e{int(enr_U5)}O"
    d_U8, d_U5, d_Th, d_O16_Fuel=ZAID_Fuel(enr_U5, f_UO2)
    # Normalization coefficient for tally volume flux
    N_U5 = (d_U5 * NA * 1E-24 * C)/M_U5
    N_U8 = (d_U8 * NA * 1E-24 * C)/M_U8
    N_Th = (d_Th * NA * 1E-24 * C)/M_Th
    N_Fuel = N_U5 + N_U8 + N_Th
    # Writing the MCNP input file
    with open(f"{filename}", "w+") as header_card:
        header_card.write(
f"""c -------------------------- Head of the input file ---------------------------
c -------------------------- MCNP Simulations input ---------------------------
c -----------------------------------------------------------------------------
c ------------------ U5 Enrichment = {enr_U5:.3f} % -----------------------------------
c ------------------ UO2 volume fraction = {f_UO2:.2f} % ------------------------------------
c ------------------ Giving a Fuel density = {Rho_Fuel:.6f} [g/cc] -----------------
c -----------------------------------------------------------------------------
c ------------------- Cell definition -----------------------------------------
c start file for MCNP project 2021
C NuScale SMR model
C  
1 3 -1.0 -3000 -2000 1000 (60:50:-40) IMP:N=1 
2 2 -7.76 40 -50 -60 (30:20:-10) IMP:N=1    
3 3 -1.0     -30 -20 10 (5:-6:7) IMP:N=1
4 2 -7.76 -6000 -5000 4000 (3000:2000:-1000)  IMP:N=1
5 0  6000:-4000:5000  IMP:N=0
7 0  -5 6 -7 IMP:N=1  FILL=6
6 0 1 -2 3 -4 IMP:N=1 U=6 LAT=1 FILL=-4:3 -4:3 0:0 & 
5 5 5 5 5 5 5 5 &
5 5 5 1 1 5 5 5 &
5 5 1 1 1 1 5 5 &
5 1 1 1 1 1 1 5 &
5 1 1 1 1 1 1 5 &
5 5 1 1 1 1 5 5 &
5 5 5 1 1 5 5 5 &
5 5 5 5 5 5 5 5  
C univers 1 (zone 1)
10 0 100 -200 300 -400 IMP:N=1 U=1 LAT=1 FILL=11
11 3 -1.0 700 IMP:N=1 U=11
12 4 -6.5 800 -700 IMP:N=1 U=11
13 0  900 -800 IMP:N=1 U=11
14 1 -{Rho_Fuel:.7f} -900 IMP:N=1 U=11
C univers 5 (rempli d'eau)
50 3 -1.0 1 -2 3 -4 IMP:N=1 U=5 $LAT=1 FILL=51

C The Surface Block
1 PX  -21.42 $ half of 21.42 cm = 17*1.26 cm
2 PX   0
3 PY  -21.42
4 PY   0
5 CZ   70 $ the closest radial dimension to the fuel assembly
6 PZ  -{H_core/2:.2f}
7 PZ   {H_core/2:.2f}
10 PZ -{H_core/2+28.00:.2f}
20 PZ  {H_core/2+28.00:.2f}
30 CZ  75 $ ~ external radius of zone4 
40 PZ -{H_core/2+38.00:.2f} 
50 PZ  {H_core/2+38.00:.2f}
60 CZ  85 $ for 10-cm of Core's baril made with stainless steel
C primary vessel of the reactor
1000 PZ -{H_core/2+122.40:.2f}
2000 PZ  {H_core/2+122.40:.2f}
3000 CZ  122 $ inner radius of the primary vessel 
4000 PZ -{H_core/2+137.40:.2f}
5000 PZ  {H_core/2+137.40:.2f}
6000 CZ  137 $ outer radis of the primary vessel
C for 7 cells per assembly side of 21.42 cm 
100 PX   0
200 PX   1.2626
300 PY   0
400 PY   1.2626
700 C/Z 0.6313 0.6313  0.4761  $ outer radius of cladding
800 C/Z 0.6313 0.6313  0.41895 $ inner radius of cladding
900 C/Z 0.6313 0.6313  0.41295 $ fuel pin cell radius

C  
C General cards   
C  
C PRDMP 2J 1
C  
C The Materials (fuel at 900 K except Pu-240 not available, D2O and Zircalloy at 300 K)
C  
C Real fuel made from a mixture of UO2 and ThO2
m1 92238.60c {-d_U8:.7f} &
   92235.60c {-d_U5:.7f} &
   90232.60c {-d_Th:.7f} &
   8016.60c  {-d_O16_Fuel:.7f} 
m2 6000.60c 0.08 14000.60c 2.00 24000.42c 19.50 25055.60c 1.50 &
   26000.42c 67.34 28000.42c 9.50
m3 1002.60c 2 8016.60c 1
mt3 hwtr.01t
m4 40000.60c 0.988548 8016.60c 0.00773128 26056.60c 0.00372108
C  Virtual materials used to calculate reaction rates
m235 92235.60c 1 8016.60c 2
m238 92238.60c 1 8016.60c 2
m232 90232.60c 1 8016.60c 2
C The Source  
C  
KSRC -33.4 33.4 0  33.4 33.4 0  -33.4 -33.4 0  33.4 -33.4 0 &
     -22 22 0  22 22 0  -22 -22 0  22 -22 0
C KCODE 4000 1 100 200
KCODE 10000 1 50 250
C Spectrum tallying
c --------------------------- 413 groups energy bining (A. NUTTIN)------------|
e0:n 1.00e-09 1.06e-09 1.12e-09 1.19e-09 1.26e-09 1.33e-09 1.41e-09 &          
1.50e-09 1.58e-09 1.68e-09 1.78e-09 1.88e-09 2.00e-09 2.11e-09 &               
2.24e-09 2.37e-09 2.51e-09 2.66e-09 2.82e-09 2.99e-09 3.16e-09 &               
3.35e-09 3.55e-09 3.76e-09 3.98e-09 4.22e-09 4.47e-09 4.73e-09 &               
5.01e-09 5.31e-09 5.62e-09 5.96e-09 6.31e-09 6.68e-09 7.08e-09 &               
7.50e-09 7.94e-09 8.41e-09 8.91e-09 9.44e-09 1.00e-08 1.06e-08 &               
1.12e-08 1.19e-08 1.26e-08 1.33e-08 1.41e-08 1.50e-08 1.58e-08 &               
1.68e-08 1.78e-08 1.88e-08 2.00e-08 2.11e-08 2.24e-08 2.37e-08 &               
2.51e-08 2.66e-08 2.82e-08 2.99e-08 3.16e-08 3.35e-08 3.55e-08 &               
3.76e-08 3.98e-08 4.22e-08 4.47e-08 4.73e-08 5.01e-08 5.31e-08 &               
5.62e-08 5.96e-08 6.31e-08 6.68e-08 7.08e-08 7.50e-08 7.94e-08 &               
8.41e-08 8.91e-08 9.44e-08 1.00e-07 1.06e-07 1.12e-07 1.19e-07 &               
1.26e-07 1.33e-07 1.41e-07 1.50e-07 1.58e-07 1.68e-07 1.78e-07 &               
1.88e-07 2.00e-07 2.11e-07 2.24e-07 2.37e-07 2.51e-07 2.66e-07 &               
2.82e-07 2.99e-07 3.16e-07 3.35e-07 3.55e-07 3.76e-07 3.98e-07 &               
4.22e-07 4.47e-07 4.73e-07 5.01e-07 5.31e-07 5.62e-07 5.96e-07 &               
6.31e-07 6.68e-07 7.08e-07 7.50e-07 7.94e-07 8.41e-07 8.91e-07 &               
9.44e-07 1.00e-06 1.06e-06 1.12e-06 1.19e-06 1.26e-06 1.33e-06 &               
1.41e-06 1.50e-06 1.58e-06 1.68e-06 1.78e-06 1.88e-06 2.00e-06 &               
2.11e-06 2.24e-06 2.37e-06 2.51e-06 2.66e-06 2.82e-06 2.99e-06 &               
3.16e-06 3.35e-06 3.55e-06 3.76e-06 3.98e-06 4.22e-06 4.47e-06 &               
4.73e-06 5.01e-06 5.31e-06 5.62e-06 5.96e-06 6.31e-06 6.68e-06 &               
7.08e-06 7.50e-06 7.94e-06 8.41e-06 8.91e-06 9.44e-06 1.00e-05 &               
1.06e-05 1.12e-05 1.19e-05 1.26e-05 1.33e-05 1.41e-05 1.50e-05 &               
1.58e-05 1.68e-05 1.78e-05 1.88e-05 2.00e-05 2.11e-05 2.24e-05 &               
2.37e-05 2.51e-05 2.66e-05 2.82e-05 2.99e-05 3.16e-05 3.35e-05 &               
3.55e-05 3.76e-05 3.98e-05 4.22e-05 4.47e-05 4.73e-05 5.01e-05 &               
5.31e-05 5.62e-05 5.96e-05 6.31e-05 6.68e-05 7.08e-05 7.50e-05 &               
7.94e-05 8.41e-05 8.91e-05 9.44e-05 1.00e-04 1.06e-04 1.12e-04 &               
1.19e-04 1.26e-04 1.33e-04 1.41e-04 1.50e-04 1.58e-04 1.68e-04 &               
1.78e-04 1.88e-04 2.00e-04 2.11e-04 2.24e-04 2.37e-04 2.51e-04 &               
2.66e-04 2.82e-04 2.99e-04 3.16e-04 3.35e-04 3.55e-04 3.76e-04 &               
3.98e-04 4.22e-04 4.47e-04 4.73e-04 5.01e-04 5.31e-04 5.62e-04 &               
5.96e-04 6.31e-04 6.68e-04 7.08e-04 7.50e-04 7.94e-04 8.41e-04 &               
8.91e-04 9.44e-04 1.00e-03 1.06e-03 1.12e-03 1.19e-03 1.26e-03 &               
1.33e-03 1.41e-03 1.50e-03 1.58e-03 1.68e-03 1.78e-03 1.88e-03 &               
2.00e-03 2.11e-03 2.24e-03 2.37e-03 2.51e-03 2.66e-03 2.82e-03 &               
2.99e-03 3.16e-03 3.35e-03 3.55e-03 3.76e-03 3.98e-03 4.22e-03 &               
4.47e-03 4.73e-03 5.01e-03 5.31e-03 5.62e-03 5.96e-03 6.31e-03 &               
6.68e-03 7.08e-03 7.50e-03 7.94e-03 8.41e-03 8.91e-03 9.44e-03 &               
1.00e-02 1.06e-02 1.12e-02 1.19e-02 1.26e-02 1.33e-02 1.41e-02 &               
1.50e-02 1.58e-02 1.68e-02 1.78e-02 1.88e-02 2.00e-02 2.11e-02 &               
2.24e-02 2.37e-02 2.51e-02 2.66e-02 2.82e-02 2.99e-02 3.16e-02 &               
3.35e-02 3.55e-02 3.76e-02 3.98e-02 4.22e-02 4.47e-02 4.73e-02 &               
5.01e-02 5.31e-02 5.62e-02 5.96e-02 6.31e-02 6.68e-02 7.08e-02 &               
7.50e-02 7.94e-02 8.41e-02 8.91e-02 9.44e-02 1.00e-01 1.06e-01 &               
1.12e-01 1.19e-01 1.26e-01 1.33e-01 1.41e-01 1.50e-01 1.58e-01 &               
1.68e-01 1.78e-01 1.88e-01 2.00e-01 2.11e-01 2.24e-01 2.37e-01 &               
2.51e-01 2.66e-01 2.82e-01 2.99e-01 3.16e-01 3.35e-01 3.55e-01 &               
3.76e-01 3.98e-01 4.22e-01 4.47e-01 4.73e-01 5.01e-01 5.31e-01 &               
5.62e-01 5.96e-01 6.31e-01 6.68e-01 7.08e-01 7.50e-01 7.94e-01 &               
8.41e-01 8.91e-01 9.44e-01 1.00e+00 1.06e+00 1.12e+00 1.19e+00 &               
1.26e+00 1.33e+00 1.41e+00 1.50e+00 1.58e+00 1.68e+00 1.78e+00 &               
1.88e+00 2.00e+00 2.11e+00 2.24e+00 2.37e+00 2.51e+00 2.66e+00 &               
2.82e+00 2.99e+00 3.16e+00 3.35e+00 3.55e+00 3.76e+00 3.98e+00 &               
4.22e+00 4.47e+00 4.73e+00 5.01e+00 5.31e+00 5.62e+00 5.96e+00 &               
6.31e+00 6.68e+00 7.08e+00 7.50e+00 7.94e+00 8.41e+00 8.91e+00 &               
9.44e+00 1.00e+01 1.06e+01 1.12e+01 1.19e+01 1.26e+01 1.33e+01 &               
1.41e+01 1.50e+01 1.58e+01 1.68e+01 1.78e+01 1.88e+01 2.00e+01
F4:N (14 < (U=11) < (U=1))
SD4 5.35729E+01
""")
    os.system(f'mcnp5 in={filename} out={Output}')
    scan_k = open(f"{Output}", "r")
    for line in scan_k:
        global k_eff, sigma, Thermal, Epithermal, Fast
        line = line.strip()  # strip end-on-line
        if re.search(r'keff = ....... w', line):  # and len(line)==7:
            target_line = line.split(' ')
            #print(f'{R}', target_line[8], f' \u00B1{target_line[15]}')
            k_eff = float(target_line[8])
            sigma = float(target_line[15])
            # print(f'{R:.2f}', target_line[8], f' \u00B1{target_line[15]}')
        # Modification of 28/11/2020: adding the fissionning neutrons
        if re.search(r'(<0.625 ev)', line):
            target_neutrons = line.split(' ')
            Thermal = f'{(target_neutrons[12])}'
            Thermal = float(Thermal.replace('%', ''))
            Epithermal = f'{target_neutrons[27]+target_neutrons[28]}'
            Epithermal = float(Epithermal.replace('%', ''))
            Fast = f'{target_neutrons[40]+target_neutrons[41]+target_neutrons[42]}'
            Fast = float(Fast.replace('%', ''))
            #print(f'{k_eff:.5f} {sigma:.5f}', Thermal, Epithermal, Fast)
    scan_k.close()
    return k_eff, sigma, Thermal, Epithermal, Fast

# header(enr_U5, f_UO2, H_core)
  
# read_Output()

# os.system(f'del runt* srct* "U{int(f_UO2)}e{int(enr_U5)}"')
# os.system(f'del "U{int(f_UO2)}e{int(enr_U5)}O" ')

file = 'Results'

def Crt_info():
    global Rslt
    Rslt = f"{file}"
    with open(f"{Rslt}", "a+") as Crt_info_file:
        Crt_info_file.write(f"""{H_core:.2f} {f_UO2:.5f} {enr_U5:.3f} {k_eff:.5f} {sigma:.5f} {Thermal:.2f} {Epithermal:.2f} {Fast:.2f}
""")
# Crt_info()

#print(f'Pour: enr_U5 = {enr_U5:.5f} f_UO2 = {f_UO2:.5f} H_core = {H_core:.5f}')


def Crt():
    global k_eff, epsilon, delta, enr_U5, f_UO2, filename, Output, H_core
    epsilon = 0.00500
    delta = 0.10
    if (k_eff > 1 + epsilon):
        enr_U5 = enr_U5 - delta
        f_UO2 = f_UO2
        H_core = 100
        os.system(f'del runt* srct* meshta* "U{int(f_UO2)}e{int(enr_U5)}"')
        os.system(f'del "U{int(f_UO2)}e{int(enr_U5)}O" ')
        header(enr_U5, f_UO2, H_core)
        Crt_info()
        Crt()   
    elif (k_eff < 1):
        enr_U5 = enr_U5 + delta
        f_UO2 = f_UO2 
        H_core = 100
        os.system(f'del runt* srct* meshta* "U{int(f_UO2)}e{int(enr_U5)}"')
        os.system(f'del "U{int(f_UO2)}e{int(enr_U5)}O" ')
        header(enr_U5, f_UO2, H_core)
        Crt_info()
        Crt()    
    else:  
        os.system(f'del runt* srct*')
        # Crt_info()
        print(f'keff = {k_eff:.5f} ,On parle d’un réacteur critique.')
        exit()
    return k_eff


header(enr_U5, f_UO2, H_core)
Crt_info()
os.system(f'del runt* srct*')

# Crt()
    

