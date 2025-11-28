# Chocked flow Calculation for a de Laval nozzle

import numpy as np

Cd = 0.72 # Discharge coefficient
A = 3.0E-3 * 14E-3 # Thoat Area (m2)
Ap = 0.25 * np.pi * (25.4E-3 * 3/4) **2 # pipe area
Pu = 689476 + 101325 # Total upstream pressure
gamma = 1.4 # Specific heat ration of the gas
Tt = 293
Ma = 1 # Mach number
M = 28.96E-3 # Molecular mass of the gas
R = 8.3145

dmdt = (A * Pu / Tt**0.5) * (gamma * M / R)**0.5 * Ma * (0.5 * (gamma + 1))**(-(gamma + 1)/(2 * (gamma - 1)))

print( " Sonic Mass flow rate at the throat: ", dmdt, " kg/s")
print( " Sonic Mass flow rate at the throat: ", dmdt * 1000 * 60, " L/min")



# Subsonic flow through a restriction with the same outlet and inlet area

dP = 750000  # Pressure difference in Pa
Q = 0.5 / 60 # Gas flow rate m3/s

print ("The power consumed by the flow rate is: ", Q * dP, " W")


# Orifice plate pressure drop Calculation
od = ((14E-3 * 2.5E-3) / (np.pi / 4))**0.5 # orifice diameter
beta = A / Ap # Area relation between pipe and orifice

rho1 = 1 # Fluid density in plane of upstream tapping
epsilon = 1 # Expansibility factor
dmdt = (Cd / (1 -  beta**4)**0.5) * epsilon * np.pi * 0.25 * od**2 * (2 * dP * rho1)**0.5

print( " The flow rate in the orifice plate is; ", dmdt, " kg/s")