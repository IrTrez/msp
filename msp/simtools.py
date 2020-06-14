from mpl_toolkits.mplot3d import Axes3D
from . import msp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
import os

AU = 149.6e6  # km
muSun = 1.327178e11 #kg
Earth = msp.Planet(398600.441, 6378.136, AU, muSun)
Mars = msp.Planet(4.282837e4, 3396.2, 1.52367934 * AU, muSun, False)


def quickAnimate(speed, dataFile, Body=None, bodyColor="cyan", plotLimits=30000):
    plt.style.use('seaborn-pastel')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    if Body is not None:
        assert (type(Body) == msp.Planet), "Incorrect Body type, should be main.Planet"
        # Sphere:
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = Body.r * np.outer(np.cos(u), np.sin(v))
        y = Body.r * np.outer(np.sin(u), np.sin(v))
        z = Body.r * np.outer(np.ones(np.size(u)), np.cos(v))
        # Plot the surface
        ax.plot_surface(x, y, z, color=str(bodyColor))
    
    ax.set_ylim(-plotLimits, plotLimits)
    ax.set_xlim(-plotLimits, plotLimits)
    ax.set_zlim(-plotLimits, plotLimits)

    data = pd.read_csv(dataFile, index_col=0)
    droppedPoints = []
    for u in range(len(data)):
        if u % speed != 0:
            droppedPoints.append(u)

    data = data.drop(droppedPoints, axis=0)
    clock = data.loc[:, "clock"].to_numpy()
    x = data.loc[:, "x"].to_numpy()
    y = data.loc[:, "y"].to_numpy()
    z = data.loc[:, "z"].to_numpy()

    manoeuvreFile = dataFile[:-4] + "_man.csv"
    showManoeuvres = True if os.path.isfile(manoeuvreFile) else False
    if showManoeuvres:
        manoeuvreData = pd.read_csv((manoeuvreFile), index_col="ID")


    plt.pause(1)
    for u in tqdm(range(len(x))):
        currentClock = clock[u]
        ax.plot(x[0:u], y[0:u], z[0:u], color="g", lw=1)
        if showManoeuvres:
            for i, manoeuver in manoeuvreData.iterrows():
                if manoeuver["clock"] > currentClock and manoeuver["clock"] < currentClock + speed:
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

        plt.pause(0.0000000001)
        plt.clf
    plt.pause(5)
