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

DATAFILE = "runs/reverseKepler.csv"
ATMOSPHEREDATA = "densityModels/MarsDensity.csv"
SPEED = 200  # __ times speed

# USE km AS STANDARD DISTANCE UNIT
# USE s AS STANDARD TIME UNIT
AU = 149.6e6  # km
muSun = 1.327178e11
currentTime = time.time()
limitAltitude = 500 # 260  #[km]. At this altitude density is just below 1*10^-10


Mars_atmosphere=m.Atmosphere(limitAltitude, densityFile=ATMOSPHEREDATA)
Earth = m.Planet(398600.441, 6378.136, AU, muSun, Mars_atmosphere)

r = np.array([0,0, -6378.136-600])
v = np.array([14,0,0])

CD = 1.2
surfaceArea = 3.6**2 * math.pi

spacecraft = m.Body(Earth, 100, CD, surfaceArea )
spacecraft.initPositionOrbit(r,v)


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
# rlist = spacecraft.propagate(40000, DATAFILE, False, dtAtmospheric = 1, dtNormal = 1)
print(np.sqrt(spacecraft.r.dot(spacecraft.r)))
# print(spacecraft.e)
# rrr, _ = spacecraft.keplerTime(r, v, -50000)
# print(spacecraft.keplerTime(r, v, -50000))
# print(np.sqrt(rrr.dot(rrr)))
# print(Earth.rsoi)



ax.set_ylim(-45000, 45000)
ax.set_xlim(-45000, 45000)
ax.set_zlim(-45000, 45000)

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
