  # .................. Criticité_Crt ...........................................
# .................. UFAS1 ...................................................
# .................. Python based script .....................................
#................... Author: S.E BENTRIDI, H. BOUZOURDAZ, H. MAKHLOUFI........
# .................. Starting Date: 07/04/2021 ...............................
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
KCODE 20000 1 50 800
C Flux volume tally calculation for real material : mixture of UO2 and ThO2
C Fission rate reaction for fuel material
F104:N (14 < (U=11) < (U=1))
SD104 5.35729E+01
FM104:N ({N_Fuel:.3e} 1 -6)
C Capture rate reaction for fuel material 
F114:N (14 < (U=11) < (U=1))
SD114 5.35729E+01
FM114:N ({N_Fuel:.3e} 1 -2)
C 
C Flux volume tally calculation for virtual material : pure UO2 with 100% of U5
C Fission rate reaction for fuel material
F504:N (14 < (U=11) < (U=1))
SD504 5.35729E+01
FM504:N ({N_U5:.3e} 235 -6)
C Capture rate reaction for fuel material
F514:N (14 < (U=11) < (U=1))
SD514 5.35729E+01
FM514:N ({N_U5:.3e} 235 -2)
C
C Flux volume tally calculation for virtual material : pure UO2 with 100% of U8
C Fission rate reaction for fuel material
F804:N (14 < (U=11) < (U=1))
SD804 5.35729E+01
FM804:N ({N_U8:.3e} 238 -6)
C Capture rate reaction for fuel material
F814:N (14 < (U=11) < (U=1))
SD814 5.35729E+01
FM814:N ({N_U8:.3e} 238 -2)
C
C Flux volume tally calculation for virtual material : pure ThO2
C Fission rate reaction for fuel material
F204:N (14 < (U=11) < (U=1))
SD204 5.35729E+01
FM204:N ({N_Th:.3e} 232 -6)
C Capture rate reaction for fuel material
F214:N (14 < (U=11) < (U=1))
SD214 5.35729E+01
FM214:N ({N_Th:.3e} 232 -2)
C Fmesh card of fission distribution central CS of the core (10cm thick)
FMESH34:n GEOM=xyz ORIGIN= -140 -140 -5
     IMESH= 140 IINTS= 280 &
     JMESH= 140 JINTS= 280 &
     KMESH= 5 KINTS= 1 
FM34:N -7.62E16 0 -2 $ Spatial distrib. of capture cst for 1MW power equivalent !!!!
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
        Crt_info_file.write(f"""{H_core:.2f} {f_UO2:.5f} {enr_U5:.3f} {k_eff:.5f} {sigma:.5f} {Thermal:.3f} {Epithermal:.3f} {Fast:.3f}
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
        os.system(f'del runt* srct* "U{int(f_UO2)}e{int(enr_U5)}"')
        os.system(f'del "U{int(f_UO2)}e{int(enr_U5)}O" ')
        header(enr_U5, f_UO2, H_core)
        Crt_info()
        Crt()   
    elif (k_eff < 1):
        enr_U5 = enr_U5 + delta
        f_UO2 = f_UO2 
        H_core = 100
        os.system(f'del runt* srct* "U{int(f_UO2)}e{int(enr_U5)}"')
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
Crt()
    

