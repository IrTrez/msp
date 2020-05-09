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
SPEED = 2.5  # __ times speed
PERIAPSEALTITUDE = 3000

# USE km AS STANDARD DISTANCE UNIT
# USE s AS STANDARD TIME UNIT
AU = 149.6e6  # km
muSun = 1.327178e11
currentTime = time.time()
limitAltitude = 260 # 260  #[km]. At this altitude density is just below 1*10^-10


Mars_atmosphere=m.Atmosphere(limitAltitude, densityFile=ATMOSPHEREDATA)
Mars = m.Planet(4.282837e4, 3396.2, 1.52367934 * AU, muSun, Mars_atmosphere)

r = np.array([0, 0, -Mars.r-PERIAPSEALTITUDE])
v = np.array([-4.518996301, 0, 0])

CD = 1.2
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
dt = -5
rlist = spacecraft.propagate(6000, dtAtmospheric = dt, dtNormal = dt)
# print("Absolute distance:", np.sqrt(spacecraft.r.dot(spacecraft.r)))
# print("Sphere of influence:", Mars.rsoi)
print("Starting vectors: ")
print("r = np.array([", spacecraft.r[0], ",", spacecraft.r[1], ",", spacecraft.r[2], "])")
print("v = np.array([", spacecraft.v[0], ",", spacecraft.v[1], ",", spacecraft.v[2], "])")




# ax.set_ylim(-50000, 50000)
# ax.set_xlim(-50000, 50000)
# ax.set_zlim(-50000, 50000)

# data = pd.read_csv(DATAFILE, names=["x", "y", "z"])
# droppedPoints = []
# for u in range(len(data)):
#     if u % SPEED != 0:
#         droppedPoints.append(u)

# data = data.drop(droppedPoints, axis=0)
# x = data.loc[:,"x"].to_numpy()
# y = data.loc[:,"y"].to_numpy()
# z = data.loc[:,"z"].to_numpy()

# # print(rlist)
# plt.pause(1)
# for u in tqdm(range(len(x))):
#     # print(deltat)
#     ax.plot(x[0:u], y[0:u], z[0:u], color="g")
#     # ax.quiver(*parentradiuslist[u], *hlist[u], color="red")
#     plt.pause(0.0000001)
#     plt.clf
# plt.pause(5)
