#libraries
import math
import numpy as np
# Base equations for Taylor-Couette (Concentric cylinder) Flow

#CONSTANTS
visc_37 = 0.0035 # Viscosity, Pa*s 
temp_factor = 1.27 # Conversion factor from 37 to 22

# Conversion from rad/s to RPM and back
def RPM(rad):
    return rad*30/math.pi
def rad(RPM):
    return RPM*math.pi/30

# Calculate shear rate in concentric cylinder flow (laminar, idealized)
# r is optional and represents distance from inner to outer (0 to 1)
def shear_Rate(r_i, r_0, w, r=1):
    r = r_i + r*(r_0-r_i)
    out = 2*(r_0**2 * r_i**2 * rad(w))   /   ((r_0**2 - r_i**2) * r**2)
    return out

# Same as above but for calculating RPM
def ang_vel(r_i, r_0, w, r=1):
    r = r_i + r*(r_0-r_i)
    out = rad(w)*((r_0**2 - r_i**2) * r**2) / 2*(r_0**2 * r_i**2 )  
    return out

# Calculate stability (Taylor Number) of flow
def Ta(r_i, r_0, w):
    return rad(w)**2*r_i*(r_0-r_i)**3/((4.566e-6)**2)

# Concentric cylinder wetted surface area approximation (slight overestimate)
def surf_Area(r_i, r_0):
    H = 0.015 # 1.5 cm
    return 0.01 * (2*H*(r_i+r_0)+(r_0**2+r_i**2) / (H*(r_0**2 - r_i**2)))

# Surface area of a tube for comparison
def surf_Tube(rad):
    return 0.01*2/rad