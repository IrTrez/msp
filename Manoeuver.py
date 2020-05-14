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

DATAFILE = "runs/ManouvreTest.csv"
ATMOSPHEREDATA = "densityModels/MarsDensity.csv"
SPEED = 600  # __ times speed

# USE km AS STANDARD DISTANCE UNIT
# USE s AS STANDARD TIME UNIT
AU = 149.6e6  # km
muSun = 1.327178e11
currentTime = time.time()
limitAltitude = 500 # 260  #[km]. At this altitude density is just below 1*10^-10


Mars_atmosphere=m.Atmosphere(limitAltitude, densityFile=ATMOSPHEREDATA)
Earth = m.Planet(398600.441, 6378.136, AU, muSun, Mars_atmosphere)

a = 10000
e = 0
i = math.radians(1)
Omega = math.radians(30.89)
omega = math.radians(80.0)
trueAnomaly = math.radians(40)

CD = 1.2
surfaceArea = 3.6**2 * math.pi

spacecraft = m.Body(Earth, 100, CD, surfaceArea )
spacecraft.initKeplerOrbit(a,e,i,Omega,omega,trueAnomaly)


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

# Add manoeuvers before propagate

spacecraft.addManouvreByDirection(spacecraft.start + 1 * spacecraft.orbitalPeriod, 6, "n")
spacecraft.addManouvreByDirection(spacecraft.start + 1 * spacecraft.orbitalPeriod, -6, "t")
spacecraft.addManouvreByDirection(spacecraft.start + 2 * spacecraft.orbitalPeriod, 1.5, "r")

# PROPAGATE Here
rlist = spacecraft.propagate(4*spacecraft.orbitalPeriod, DATAFILE, False)

ax.set_ylim(-30000, 30000)
ax.set_xlim(-30000, 30000)
ax.set_zlim(-30000, 30000)

data = pd.read_csv(DATAFILE, index_col=0)
droppedPoints = []
for u in range(len(data)):
    if u % SPEED != 0:
        droppedPoints.append(u)

data = data.drop(droppedPoints, axis=0)
clock = data.loc[:,"clock"].to_numpy()
x = data.loc[:,"x"].to_numpy()
y = data.loc[:,"y"].to_numpy()
z = data.loc[:,"z"].to_numpy()

manoeuvreData = pd.read_csv((DATAFILE[:-4] + "_man.csv"), index_col="ID").loc[:,["clock", "manType","r"]]

# print(rlist)
plt.pause(1)
for u in tqdm(range(len(x))):
    currentClock=clock[u]
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
            print(manoeuver["r"])
            manoeuverPos = manoeuver["r"][1:-1].split(" ")
            manoeuverPos = [float(x) for x in manoeuverPos if x]
            ax.scatter(*manoeuverPos, marker = manMarker, color = manColor)
    # ax.quiver(*parentradiuslist[u], *hlist[u], color="red")
    plt.pause(0.0000001)
    plt.clf
plt.pause(5)
