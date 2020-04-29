from mpl_toolkits.mplot3d import Axes3D
import main as m
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
from tqdm import tqdm
import math
from math import exp
plt.style.use('seaborn-pastel')

DATAFILE = "runs/propagateTest.csv"
ATMOSPHEREDATA = "MarsDensity.csv"
SPEED = 200  # __ times speed

# USE km AS STANDARD DISTANCE UNIT
# USE s AS STANDARD TIME UNIT
AU = 149.6e6  # km
muSun = 1.327178e11
currentTime = time.time()
limitAltitude = 500 # 260  #[km]. At this altitude density is just below 1*10^-10


Mars_atmosphere=m.Atmosphere(limitAltitude, densityFile=ATMOSPHEREDATA)
Earth = m.Planet(398600.441, 6378.136, AU, muSun, Mars_atmosphere)

p = 18067.790
e = 0.58285
i = math.radians(45)
Omega = math.radians(30.89)
omega = math.radians(53.38)
trueAnomaly = math.radians(40)
a = p / (1 - e**2)

CD = 1.2
surfaceArea = 3.6**2 * math.pi

spacecraft = m.Body(Earth, 100, CD, surfaceArea )
spacecraft.initPositionOrbit()

print(spacecraft.apoapsis)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Sphere:
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = Earth.r * np.outer(np.cos(u), np.sin(v))
y = Earth.r * np.outer(np.sin(u), np.sin(v))
z = Earth.r * np.outer(np.ones(np.size(u)), np.cos(v))
# Plot the surface
ax.plot_surface(x, y, z, color='tab:cyan')


# PROPAGATE Here
rlist = spacecraft.propagate(2*spacecraft.orbitalPeriod, DATAFILE, True, dtAtmospheric = -1, dtNormal = -1)

ax.set_ylim(-30000, 30000)
ax.set_xlim(-30000, 30000)
ax.set_zlim(-30000, 30000)

data = pd.read_csv(DATAFILE, names=["x", "y", "z"])
droppedPoints = []
for u in range(len(data)):
    if u % SPEED != 0:
        droppedPoints.append(u)

data = data.drop(droppedPoints, axis=0)
x = data.loc[:,"x"].to_numpy()
y = data.loc[:,"y"].to_numpy()
z = data.loc[:,"z"].to_numpy()

# print(rlist)
plt.pause(1)
for u in tqdm(range(len(x))):
    # print(deltat)
    ax.plot(x[0:u], y[0:u], z[0:u], color="g")
    # ax.quiver(*parentradiuslist[u], *hlist[u], color="red")
    plt.pause(0.0000001)
    plt.clf
plt.pause(2)