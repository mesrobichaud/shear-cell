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
def ang_vel(r_i, r_0, shear, r=1):
    r = r_i + r*(r_0-r_i)
    out = shear*((r_0**2 - r_i**2) * r**2) / 2*(r_0**2 * r_i**2 )  
    return RPM(out)

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

# native sizing of outer cylinder for 24 well plate 
# height is height of center of cone over flat surface
# hmax is max height of fluid in well
# alpha in radians
# OUTPUT in microlitres 
def volume_24well(r_i,alpha,hmax):
    from scipy.integrate import quad
    
    beta = math.atan(0.3/17) #eppendorf draft angle
    height = (r_i/16.2/2)*(1-math.tan(alpha))
    
    r_omax = 16.2 + hmax*math.tan(beta)
    
    def int1(r):
        return (height + r*math.tan(alpha))*2*math.pi*r
    
    def int2(r): 
        return hmax*2*math.pi*r
    
    def int3(r):
        return (hmax - (r-16.2/2)/math.tan(beta))*2*math.pi*r
    
    return quad(int1,0,r_i)[0] + quad(int2,r_i,16.2/2)[0] + quad(int3, 16.2/2,r_omax/2)[0]

# accepts array of radii from centre of rotation and returns shear rates
# omega(r)/h
# omega is given in RPM for ease of use
# r and h must be in same units
def disk_Shear(omega,r,h):
    return rad(omega)*r/h
    
# accepts array of radii from centre of rotation and returns shear rates
# omega(r)/h
# omega is given in RPM for ease of use
# r and h must be in same units
# SDM stepped disk model
# h1 is inner
# ratio: of surface areas, zone 1 : zone 2
def SDM_Shear(omega,r,h1,h2,frac):
    
    
    ratio = frac / (1 - frac)
    
    r_star = math.sqrt(4*ratio + max(r)**2) / math.sqrt(ratio+1)
    
    out = []
    
    for rs in r: 
        if (rs < r_star):
            out.append(rad(omega)*rs/h1)
        else: 
            out.append(rad(omega)*rs/h2)
    return np.array(out)