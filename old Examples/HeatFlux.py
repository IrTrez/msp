import math
from math import sqrt, cos, pi
import time
import numpy as np

AU = 149.6e6  # km
muSun = 1.327178e11
muMars = 42828.37
currentTime = time.time()

AU = 149.6e6  # km
r_a=1.6660*AU
r_p=1.3814*AU
e_m=(r_a-r_p)/(r_a+r_p)
a=0.5*(r_a+r_p)
i_m=25*pi/180 #rot. axis inclination
F_0=590 #W/m^2
L=3.9e26 #[W]
T=2*pi*sqrt(a**3/muMars)
n=sqrt(muMars/a**3) #mean motion, useful for the longitude calculations [rad/s]

theta_tab=np.arange(0,2*np.pi)
flux_tab=[]
rtab=[]

def Flux(theta, l): #with clear sky #heat flux at midday
    """Find the heat flux from the sun at any point of the martian year and at any latitude
    
    Arguments:
        theta  -- true anomaly [rad]
        l  -- latitude angle [rad]
    
    Returns:
        flux -- Heat flux [W/m^2]
    """
    r=(a*(1-e_m**2))/(1+e_m*cos(theta))
    flux_str = L / (4*pi*(r*10**3)**2) #Now we need to adjust for latitude and axis inclination
    #tilt_eff=i_m*cos(theta)
    flux=flux_str*cos(l)
    return flux, r

for i, theta in enumerate(theta_tab):
    flux_tab.append(Flux(theta,0)[0])
    rtab.append(Flux(theta,0)[1])

plt.polar(theta_tab, rtab)