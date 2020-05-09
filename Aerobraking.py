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
SPEED = 500  # __ times speed
RUNTIME = 3000000
print("Total runtime will be:", RUNTIME, "s or:", RUNTIME/3600, "hours or:", RUNTIME/86400, "days")

# USE km AS STANDARD DISTANCE UNIT
# USE s AS STANDARD TIME UNIT
AU = 149.6e6  # km
muSun = 1.327178e11
currentTime = time.time()
limitAltitude = 260 # 260  #[km]. At this altitude density is just below 1*10^-10


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
dt = 4
# spacecraft.AddManoeuverByDirection()
# spacecraft.addManouvreByDirection(spacecraft.start + 1, 2, "r")
spacecraft.addManouvreByDirection(spacecraft.start + 100, -1.332, "t")
spacecraft.addManouvreByDirection(spacecraft.start + 8900, -1.1, "t")

rlist = spacecraft.propagate(RUNTIME, DATAFILE, True, dtAtmospheric = dt, dtNormal = dt)
# print("Absolute distance:", np.sqrt(spacecraft.r.dot(spacecraft.r)))
# print("Sphere of influence:", Mars.rsoi)
# print(spacecraft.r)
print(np.sqrt(spacecraft.v.dot(spacecraft.v)))
print(spacecraft.e)
print(spacecraft.periapsis)
print("Periaosis alt", spacecraft.periapsis-Mars.r)
print("Apoapsis alt", spacecraft.apoapsis-Mars.r)



ax.set_ylim(-50000, 50000)
ax.set_xlim(-50000, 50000)
ax.set_zlim(-50000, 50000)


data = pd.read_csv(DATAFILE, index_col=0)
droppedPoints = []
for u in range(len(data)):
    if u % SPEED != 0:
        droppedPoints.append(u)

data = data.drop(droppedPoints, axis=0)
clock = data.loc[:, "clock"].to_numpy()
x = data.loc[:, "x"].to_numpy()
y = data.loc[:, "y"].to_numpy()
z = data.loc[:, "z"].to_numpy()

manoeuvreData = pd.read_csv(
    (DATAFILE[:-4] + "_man.csv"), index_col="ID").loc[:, ["clock", "manType", "r"]]

# print(rlist)
plt.pause(1)
for u in tqdm(range(len(x))):
    currentClock = clock[u]
    ax.plot(x[0:u], y[0:u], z[0:u], color="g")
    for i, manoeuver in manoeuvreData.iterrows():
        if manoeuver["clock"] > currentClock and manoeuver["clock"] < currentClock+SPEED:
            if manoeuver["manType"] == "t":
                manMarker = "o"
                manColor = "lime"
            if manoeuver["manType"] == "r":
                manMarker = "o"
                manColor = "cyan"
            if manoeuver["manType"] == "n":
                manMarker = "^"
                manColor = "m"
            # general manoeuvre
            if manoeuver["manType"] == "g":
                manMarker = "+"
                manColor = "r"
            manoeuverPos = manoeuver["r"][1:-1].split(" ")
            manoeuverPos = [float(x) for x in manoeuverPos if x]
            ax.scatter(*manoeuverPos, marker=manMarker, color=manColor)
    # ax.quiver(*parentradiuslist[u], *hlist[u], color="red")
    plt.pause(0.000000001)
    plt.clf
plt.pause(5)
