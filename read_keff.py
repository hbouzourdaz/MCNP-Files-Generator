# .................. MCNP Material card generator ............................
# .................. Oklo type Reaction Zone .................................
# .................. Python based script .....................................
#................... Author: S.E BENTRIDI ....................................
# .................. Starting Date: 08/09/2020 ...............................
# .................. present release : 0.4 ...................................
# .................. Last modif. date: 21/11/2020 ............................
import math
import os
import re
import platform
import numpy as np
from Card_DATA import *
from decimal import *
Output = 'U30e5'

def read_Output():
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
read_Output()
x, y, z, u, v = read_Output()
print(f'{k_eff:.5f} {sigma:.5f}', Thermal, Epithermal, Fast)
print(x, y, z, u, v)
