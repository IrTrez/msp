from mpl_toolkits.mplot3d import Axes3D
import main as m
import numpy as np
import matplotlib.pyplot as plt
import time
from tqdm import tqdm
import math
from matplotlib.animation import FuncAnimation
import pandas as pd
import mpl_toolkits.mplot3d.axes3d as p3
import os.path
plt.style.use('dark_background')

# INPUT
DATAFILE = "runs/Aerobraking.csv"
SPEED = 500 # __ times speed

# USE km AS STANDARD DISTANCE UNIT
# USE s AS STANDARD TIME UNIT
AU = 149.6e6  # km
muSun = 1.327178e11
Mars = m.Planet(4.282837e4, 3396.2, 1.52367934 * AU, muSun, False)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')

# Sphere:
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = Mars.r * np.outer(np.cos(u), np.sin(v))
y = Mars.r * np.outer(np.sin(u), np.sin(v))
z = Mars.r * np.outer(np.ones(np.size(u)), np.cos(v))
# Plot the surface
ax.plot_surface(x, y, z, color='tab:orange')


line, = ax.plot([], [], lw=2)
scatt, = ax.plot([], [], linestyle="", marker="o", color="white")

ax.set_ylim(-40000, 40000)
ax.set_xlim(-40000, 40000)
ax.set_zlim(-40000, 40000)
plt.axis("off")

manoeuvreFileAvailable = os.path.isfile(DATAFILE[:-4] + "_man.csv")

data = pd.read_csv(DATAFILE, index_col=0).loc[:, ["x", "y", "z"]]
clock = pd.read_csv(DATAFILE, index_col=0).loc[:, ["clock"]]
if manoeuvreFileAvailable:
    manoeuvreData = pd.read_csv(
    (DATAFILE[:-4] + "_man.csv"), index_col="ID").loc[:, ["ID", "clock", "r"]]

print(len(data))

droppedPoints = []
for u in range(len(data)):
    if u % SPEED != 0:
        droppedPoints.append(u)

data = data.drop(droppedPoints, axis=0).to_numpy()
print(len(data))

savename = DATAFILE[5:-4]
headerText = savename + " " + str(SPEED) + "x speed"
title = ax.set_title(headerText)
steps = len(data)

def animate(i):
    if manoeuvreFileAvailable:
        # MANOEUVRE
        manoeuverPosList = []
        manoeuverIDList = []
        for d, manoeuver in manoeuvreData.iterrows():
            if int((manoeuver["clock"] - clock.to_numpy()[0][0])/SPEED) < i:
                manoeuverPos = manoeuver["r"][1:-1].split(" ")
                manoeuverPos = [float(x) for x in manoeuverPos if x]
                if manoeuver["ID"] not in manoeuverIDList:
                    manoeuverPosList.append(manoeuverPos)
                    manoeuverIDList.append(manoeuver["ID"])

        if len(manoeuverPosList) > 0:
            mandisplaydata = np.array(manoeuverPosList).T
            scatt.set_data(mandisplaydata[0], mandisplaydata[1])
            scatt.set_3d_properties(mandisplaydata[2])

        if i == 0:
            scatt.set_data(100000, 100000)
            scatt.set_3d_properties(100000)

    # TRAJECTORY
    displaydata = np.array(data[0: i+1]).T
    line.set_data(displaydata[0], displaydata[1])
    line.set_3d_properties(displaydata[2])
    # title.set_text(headerText + ' radius = {}'.format(round(np.sqrt(displaydata[0].dot(displaydata[0])),0)))
    return title, line, scatt,


anim = FuncAnimation(fig, animate, frames= steps, interval=10, blit=True)

fullSaveName = "animations/" + savename + ".mp4"
print("Saving to " + fullSaveName)
anim.save(fullSaveName, fps=30, extra_args=['-vcodec', 'libx264'],dpi=300)
# plot the orbit
# ax.plot(rlist[0], rlist[1], rlist[2], color="g")

plt.show()
