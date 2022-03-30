# Physics and chemical data needed to generate materials in MCNP input files
#................... Author: S.E BENTRIDI ....................................
# .................. Starting Date: 31/03/2020 ...............................
# .................. Lockdown because of COVID-19 ............................
# .................. present release : 0.1 ...................................
# .................. Last modif. date: 31/03/2020 ............................
import math

# Natural U235 enrichment is given by:
Age = 2.0e+9
enr_U5 = 100 / (1 + 137.84261 * math.exp(-8.29517579948633e-10 * Age))

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

# Uraninite definition:
M_UO2 = (M_O16*2)+0.01*(M_U5*enr_U5+M_U8*(100-enr_U5))
M_ThO2 = (M_O16*2)+M_Th

# Sandstone : Silica and clay Molar masses
M_Si28 = 28.08550 # Natural Silicon
M_Al27 = 26.08600
M_Mg12 = 24.30500 # Natural Magnesium
M_Fe56 = 55.84500 # Natural Iron
M_Na23 = 23.00000
M_K39 = 39.09800
M_Zr = 90.4360 # Natural Zirconium
# Neutron Poisons Molar masses :
# Boron 10
M_B10 = 9.92690
# Samarium
M_Sm144 = 143.91199 
Ab_Sm144 = 3.07 #%
M_Sm147 = 146.91489
Ab_Sm147 = 14.99 #%
M_Sm148 = 147.91482
Ab_Sm148 = 11.24 #%
M_Sm149 = 148.91718
Ab_Sm149 = 13.82 #%
M_Sm150 = 149.91727
Ab_Sm150 = 7.38 #%
M_Sm152 = 151.91973
Ab_Sm152 = 26.75 #%
M_Sm154 = 153.92221
Ab_Sm154 = 22.75 #%
M_Sm_nat = 150.36635
# Gadolinium
M_Gd152 = 151.91979
Ab_Gd152 = 0.20 #%
M_Gd154 = 153.92086
Ab_Gd154 = 2.18 #%
M_Gd155 = 154.92262
Ab_Gd155 = 14.80 #%
M_Gd156 = 155.92212
Ab_Gd156 = 20.47 #%
M_Gd157 = 156.92396
Ab_Gd157 = 15.65 #%
M_Gd158 = 157.92410
Ab_Gd158 = 24.84 #%
M_Gd160 = 159.92705
Ab_Gd160 = 21.86 #%
M_Gd_nat = 155.8891
# basical Compound densities [g/cc]
Rho_UO2 = 10.60000 # Uraninite [g/cc]
Rho_ThO2 = 9.86 # Thorite [g/cc]
Rho_Wtr0 = 1.00 # Water at normal P,T conditions [g/cc]
Rho_SiO2 = 2.65000 # Silica or Quartz [g/cc]
Rho_Chlorite = 2.70000 # Chlorite [g/cc]
Rho_Illite = 2.80000 # Illite [g/cc]
Rho_OM1 = 0.77 # Alkaline [g/cc]
Rho_OM2 = 0.88 # Aromatic [g/cc]
# First degree mixture of Chlorite and Illite elements
# Chlorite
Ch_Si = 0.2237486313
Ch_Al = 0.1688436448
Ch_Mg = 0.0109945982
Ch_Fe = 0.0196489396
Ch_Na = 0.0017340448
Ch_K = 0.0707454107
Ch_O = 0.4965591964
Ch_H = 0.0077255342
# Illite
Ill_Si = 0.1305860504
Ill_Al = 0.1407521405
Ill_Mg = 0.0409500054
Ill_Fe = 0.2605657533
Ill_Na = 0.0007452171
Ill_K = 0.0019002065
Ill_O = 0.4147036135
Ill_H = 0.0097970132

# S(a,b) Treatment
S_treatment = ['nothing','lwtr.61t']
S = S_treatment.index('lwtr.61t')
