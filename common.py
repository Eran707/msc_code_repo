# -*- coding: utf-8 -*-
"""
"""
# constants
# R with different units
# 8.31451                   J K-1 mol-1
# 8.20578 x 10-2            L atm K-1 mol-1
# 8.31451 x 10-2                L bar K-1 mol-1
# 8.31451                           Pa m3 K-1 mol-1
# 62.364                             L Torr K-1 mol-1
# 1.98722                           cal K-1 mol-1
R = 8.31446
F: float = 96485.33  # Faraday's constant        C mol-1
k = 1.38e-23  # Boltzmann constant        J K-1
q = 1.602176620898e-19  # elementary charge         C
Na = 6.022e23  # Avogadro's constant       mol-1

T = 37 + 273.15
RTF = R * T / F
RT = R * T
# permeabilities
gna = (2e-3) / (F)
gk = (7e-3) / (F)
gcl = (2e-3) / (F)
gx = 1e-8  # gna,gk,gcl: conductances in mS/cm^2 conv to S/dm^2 (10^-3/10^-2) - corrected for neuron
p_kcc2 = (2e-3) / (F)
p_atpase = 0.1 / F  # C/(dm2·s)

# concentrations
nao = 145e-3
clo = 119e-3
ko = 3.5e-3
xo = 29.5e-3
# xo = -1.0 * (clo - nao - ko)  # nao,clo,ko,xo: extracellular concentrations (mM converted to M)
xo_z = xo * 0.02
oso = xo + nao + clo + ko

vw = 0.018  # partial molar volume of water, dm3/mol
pw = 0.0015  # osmotic permeability, biological membrane (muscle? unknown), dm s
km = 5 * 10 ** (-14)  # extensional rigidity of RBC at 23 deg, Mohandas and Evans (1994), N/dm

# cm = 2e-4 #default membrane capacitance (F/dm^2)
cm = 2e-4

val = {"na": 1, "k": 1, "cl": -1, "x": -0.85}

diff_constants = {"na": 1.33e-7, "k": 1.96e-7, "cl": 2.03e-7,
                  "x": 0}  # diffusion coefficients for the various ions in dm2/s


# diff_constants = {"na" : 1.33e-8, "k": 1.96e-8, "cl":2.03e-8, "x":0}


hhv0 = 0
hhv1 = 100*1e-3*1e-3  # mV.ms --> v.s
hhv2 = -25*1e-3  # mV -->v
hhv3 = 10*1e-3  # mV --v
hhv4 = 0.125*1e3  # ms^-1 --> s^-1
hhv5 = -72.6*1e-3 # mV --> v
hhv6 = 80*1e-3  # mV --> v
hhv7 = 10*1e-3*1e-3  # mV.ms --> v.s
hhv8 = -25*1e-3 # mV --> v
hhv9 = 10*1e-3 # mV --> v
hhv10 = 4*1e3  # ms^-1 --> s^-1
hhv11 = -72.6*1e-3 # mV --> v
hhv12 = 18*1e-3 # mV --> v
hhv13 = 0.07*1e3  # ms^-1 --> s^-1
hhv14 = -72.6*1e-3 # mV --> v
hhv15 = 20*1e-3 # mV --> v
hhv16 = 1*1e-3  # ms --> s
hhv17 = -35*1e-3 # mV --> v
hhv18 = 10*1e-3 # mV --> v

hh_v = (
    hhv0, hhv1, hhv2, hhv3, hhv4, hhv5, hhv6, hhv7, hhv8, hhv9, hhv10, hhv11, hhv12, hhv13, hhv14, hhv15, hhv16, hhv17,
    hhv18)

gk_bar = (36 * 1e-3 / 1e-2) / F  # mS/cm^2 to S/dm^2
gna_bar = (120 * 1e-3 / 1e-2) / F  # mS/cm^2 to S/dm^2 # mS/cm^2

m = -0.0809555699999999 #for adjustment of x based on z
b = 0.08697060950000013 #for adjustment of x based on z

