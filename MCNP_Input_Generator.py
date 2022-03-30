# .................. MCNP Material card generator ............................
# .................. Oklo type Reaction Zone .................................
# .................. Python based script .....................................
#................... Author: S.E BENTRIDI ....................................
# .................. Starting Date: 08/09/2020 ...............................
# .................. present release : 0.1 ...................................
# .................. Last modif. date: 08/09/2020 ............................
import math
import os
import re
import numpy as np
from Card_DATA import *
from decimal import *
from Material_Card_MCNP import *

getcontext().prec = 12

# Pressure and Temperature values
#P = 200
#T = 150

# The water density
#def Water_density(P*, T):
#    global Rho_wtr
#    Rho_wtr = -3.51e-6*((T+273.15)**2)+2.01e-3*(T+273.15)+0.7125
#    return Rho_wtr

# def room_tmp(T):
#     KT = 8.617e-11 * (T + 273.15)
#     tmp = (N_Cells + 1) * f'{KT:.3e} '
#     return KT


# Reading of enrichment and UO2 fraction given by the user:
# float is added to convert charaters into digits
#enr_U5 = float(input('Give the value of U5 enrichment in [%]:..... '))
#f_UO2 = float(input('Give the value of UO2 volume fraction [%] in fuel pellet: ...... '))
#H_core = float(input('Give the active hight of the core [cm]:..... '))

# Files naming protocol
filename = f"U{int(f_UO2)}e{int(enr_U5)}"
Output = f"U{int(f_UO2)}e{int(enr_U5)}O"

# Header Definition
def header(enr_U5, f_UO2, H_core):
    global filename, Output
    Rho_Fuel=Fuel(enr_U5, f_UO2)
    d_U8, d_U5, d_Th, d_O16_Fuel=ZAID_Fuel(enr_U5, f_UO2)
    with open(f"{filename}", "w+") as header_card:
        header_card.write(
f"""c -------------------------- Head of the input file ---------------------------
c -------------------------- MCNP Simulations input ---------------------------
c -----------------------------------------------------------------------------
c ------------------ U5 Enrichment = {enr_U5:.5f} % -----------------------------------
c ------------------ UO2 volume fraction = {f_UO2:.5f} % ------------------------------------
c ------------------ Giving a Fuel density = {Rho_Fuel:.6f} [g/cc] -----------------
c -----------------------------------------------------------------------------
c ------------------- Cell definition -----------------------------------------
c start file for MCNP project 2018
C 4 radial fuel zones (with same fuel composition here)
C surrounded by 20-cm (D2O) axial and radial reflectors
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
C les univers 1 et 4 sont definis ici comme des reseaux
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
5 CZ   70
6 PZ  -{H_core/2:.2f}
7 PZ   {H_core/2:.2f}
10 PZ -{H_core/2+28.00:.2f}
20 PZ  {H_core/2+28.00:.2f}
30 CZ  75 $ ~ external radius of zone4 
40 PZ -{H_core/2+38.00:.2f} 
50 PZ  {H_core/2+38.00:.2f}
60 CZ  85 $ for 20-cm radial reflector
C another cylinder
1000 PZ -{H_core/2+122.40:.2f}
2000 PZ  {H_core/2+122.40:.2f}
3000 CZ  122
4000 PZ -{H_core/2+137.40:.2f}
5000 PZ  {H_core/2+137.40:.2f}
6000 CZ  137
C for 7 cells per assembly side of 21.42 cm (factor of 17/7 on standard dimensions) 
100 PX   0
200 PX   1.2626
300 PY   0
400 PY   1.2626
700 C/Z 0.6313 0.6313  0.4761  $
800 C/Z 0.6313 0.6313  0.41895 $
900 C/Z 0.6313 0.6313  0.41295

C  
C General cards   
C  
C PRDMP 2J 1
C  
C The Materials (fuel at 900 K except Pu-240 not available, D2O and Zircalloy at 300 K)
C  
m1 92238.60c {-d_U8:.7f} &
   92235.60c {-d_U5:.7f} &
   90232.60c {-d_Th:.7f} &
   8016.60c  {-d_O16_Fuel:.7f} 
m2 6000.60c 0.08 14000.60c 2.00 24000.42c 19.50 25055.60c 1.50 &
   26000.42c 67.34 28000.42c 9.50
m3 1002.60c 2 8016.60c 1
mt3 hwtr.01t
m4 40000.60c 0.988548 8016.60c 0.00773128 26056.60c 0.00372108
C  
C The Source  
C  
KSRC -33.4 33.4 0  33.4 33.4 0  -33.4 -33.4 0  33.4 -33.4 0 &
     -22 22 0  22 22 0  -22 -22 0  22 -22 0
KCODE 1000 1 20 100
C f4:N 2
FMESH14:n GEOM=xyz ORIGIN= -75 -75 -5 &
IMESH= 75 IINTS= 150 &
JMESH= 75 JINTS= 150 &
KMESH= 5 KINTS= 1 
FM14 6.022E23
FM14:N
""")


# execution of the header function which will create a file
#header(enr_U5, f_UO2, H_core)
#os.system(f'mcnp5 in={filename} out={Output}')
#os.system(f'del runt* srct*')