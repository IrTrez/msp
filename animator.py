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
plt.style.use('dark_background')

# INPUT
DATAFILE = "runs/ManouvreTest.csv"
SPEED = 300 # __ times speed

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
ax.set_ylim(-10000, 10000)
ax.set_xlim(-10000, 10000)
ax.set_zlim(-10000, 10000)
# ax.grid(False)
plt.axis("off")

data = pd.read_csv(DATAFILE, names=["x", "y", "z"])
print(len(data))

droppedPoints = []
for u in range(len(data)):
    if u % SPEED != 0:
        droppedPoints.append(u)

data = data.drop(droppedPoints, axis=0).to_numpy()
print(len(data))
savename = DATAFILE[5:-4]
plt.suptitle(savename + " " + str(SPEED) + "x speed")

steps = len(data)
def animate(i):
    displaydata = np.array(data[0: i+1]).T
    line.set_data(displaydata[0], displaydata[1])
    line.set_3d_properties(displaydata[2])
    return line,


anim = FuncAnimation(fig, animate, frames= steps, interval=10, blit=True)

fullSaveName = "animations/" + savename + ".mp4"
print("Saving to " + fullSaveName)
anim.save(fullSaveName, fps=30, extra_args=['-vcodec', 'libx264'])
# plot the orbit
# ax.plot(rlist[0], rlist[1], rlist[2], color="g")

plt.show()
