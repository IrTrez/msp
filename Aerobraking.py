import main as m
import simtools
import time
from tqdm import tqdm
import math
import numpy as np

DATAFILE = "runs/Aerobraking.csv"
ATMOSPHEREDATA = "densityModels/MarsDensity.csv"
SPEED = 300  # __ times speed
RUNTIME = 10000
print("Total runtime will be:", RUNTIME, "s or:", RUNTIME/3600, "hours or:", RUNTIME/86400, "days")

# USE km AS STANDARD DISTANCE UNIT
# USE s AS STANDARD TIME UNIT
AU = 149.6e6  # km
muSun = 1.327178e11
currentTime = time.time()
limitAltitude = 200 # 260  #[km]. At this altitude density is just below 1*10^-10


Mars_atmosphere=m.Atmosphere(limitAltitude, densityFile=ATMOSPHEREDATA)
Mars = m.Planet(4.282837e4, 3396.2, 1.52367934 * AU, muSun, Mars_atmosphere)

r = np.array([21508.114845629447, 0.0, 992.3450283462487])
v = np.array([-2.968111925169866, 0.0, -1.4808260236254678])

CD = 1.23 * 50
surfaceArea = 3.6**2 * math.pi

spacecraft = m.Body(Mars, 3000, CD, surfaceArea )
spacecraft.initPositionOrbit(r,v)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Sphere:
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = Mars.r * np.outer(np.cos(u), np.sin(v))
y = Mars.r * np.outer(np.sin(u), np.sin(v))
z = Mars.r * np.outer(np.ones(np.size(u)), np.cos(v))
# Plot the surface
ax.plot_surface(x, y, z, color='tab:orange')


# PROPAGATE Here
dt = 1
# spacecraft.AddManoeuverByDirection()
# spacecraft.addManouvreByDirection(spacecraft.start + 1, 2, "r")
spacecraft.addManouvreByDirection(spacecraft.start + 100, -1.332, "t")
spacecraft.addManouvreByDirection(spacecraft.start + 8900, -1.1, "t")

rlist = spacecraft.propagate(RUNTIME, DATAFILE, True, dtAtmospheric = dt, dtNormal = dt)

print(np.sqrt(spacecraft.v.dot(spacecraft.v)))
print(spacecraft.e)
print(spacecraft.periapsis)
print("Periaosis alt", spacecraft.periapsis-Mars.r)
print("Apoapsis alt", spacecraft.apoapsis-Mars.r)

simtools.quickAnimate(SPEED,DATAFILE,Mars)