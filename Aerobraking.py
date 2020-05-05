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

DATAFILE = "runs/Aerobraking.csv"
ATMOSPHEREDATA = "densityModels/MarsDensity.csv"
SPEED = 250  # __ times speed

# USE km AS STANDARD DISTANCE UNIT
# USE s AS STANDARD TIME UNIT
AU = 149.6e6  # km
muSun = 1.327178e11
currentTime = time.time()
limitAltitude = 260 # 260  #[km]. At this altitude density is just below 1*10^-10


Mars_atmosphere=m.Atmosphere(limitAltitude, densityFile=ATMOSPHEREDATA)
Mars = m.Planet(4.282837e4, 3396.2, 1.52367934 * AU, muSun, Mars_atmosphere)

r = np.array([19466.487519400343, 0.0, 8689.063496665178])
v = np.array([-2.3959379440133333, 0.0, -2.0502993355905814])

CD = 1.23
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
dt = 8
# spacecraft.AddManoeuverByDirection()
rlist = spacecraft.propagate(14500, DATAFILE, True, dtAtmospheric = dt, dtNormal = dt)
# print("Absolute distance:", np.sqrt(spacecraft.r.dot(spacecraft.r)))
# print("Sphere of influence:", Mars.rsoi)
# print(spacecraft.r)
print(np.sqrt(spacecraft.v.dot(spacecraft.v)))
print(spacecraft.e)



ax.set_ylim(-50000, 50000)
ax.set_xlim(-50000, 50000)
ax.set_zlim(-50000, 50000)

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
plt.pause(5)
