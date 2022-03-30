# .................. MCNP Material card generator ............................
import math
import os
import re
import numpy as np
from Card_DATA import *
from decimal import *

# digital precision
getcontext().prec = 12

# def room_tmp(T):
#     KT = 8.617e-11 * (T + 273.15)
#     tmp = (N_Cells + 1) * f'{KT:.3e} '
#     return KT

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
    
# Reading of enrichment and UO2 fraction given by the user:
# float is added to convert charaters into digits
# enr_U5 = float(input('Give the value of U5 enrichment in [%]:..... '))
# f_UO2 = float(input('Give the value of UO2 volume fraction [%] in fuel pellet: ...... '))

# Execution of both functions with given parameters enr_U5 and f_UO2
# Fuel(enr_U5, f_UO2)
# ZAID_Fuel(enr_U5, f_UO2)
# print(f"""
# --------------------------------------------------------------------------
# m1 8016.60c  {-d_O16_Fuel:.7f} &
#    92238.60c {-d_U8:.7f} &
#    92235.60c {-d_U5:.7f} &
#    90232.60c {-d_Th:.7f}
# --------------------------------------------------------------------------
# """)

