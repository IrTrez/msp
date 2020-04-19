from mpl_toolkits.mplot3d import Axes3D
import main as m
import numpy as np
import matplotlib.pyplot as plt
import time
from tqdm import tqdm
import math
from matplotlib.animation import FuncAnimation
import mpl_toolkits.mplot3d.axes3d as p3
plt.style.use('seaborn-pastel')

# USE km AS STANDARD DISTANCE UNIT
# USE s AS STANDARD TIME UNIT
AU = 149.6e6  # km
muSun = 1.327178e11
currentTime = time.time()
Mars = m.Planet(42800, 3402, 1.524 * AU, muSun)
Earth = m.Planet(398600.441, 6378.136, AU, muSun)

p = 15067.790
e = 0.53285
i = math.radians(35)
Omega = math.radians(30.89)
omega = math.radians(53.38)
trueAnomaly = math.radians(0)
a = p / (1 - e**2)

initTest = m.Body(Earth, 100)
initTest.initKeplerOrbit(a,e,i,Omega,omega,trueAnomaly)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')


# Sphere:
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = Earth.r * np.outer(np.cos(u), np.sin(v))
y = Earth.r * np.outer(np.sin(u), np.sin(v))
z = Earth.r * np.outer(np.ones(np.size(u)), np.cos(v))
# Plot the surface
ax.plot_wireframe(x, y, z, color='b')

timestep = 200
rlist = []
for deltat in tqdm(range(int((initTest.orbitalPeriod * 1.1)/ timestep) + 1)):
    initTest.refreshByTimestep(timestep)
    rlist.append(initTest.r)
rlist = np.array(rlist)

line, = ax.plot([], [], lw=3)
ax.set_ylim(-20000, 20000)
ax.set_xlim(-20000, 20000)
ax.set_zlim(-20000, 20000)


# currenttime = 0
def animate(i):
    displaydata = np.array(rlist[0: i+1]).T
    # print(displaydata)

    line.set_data(displaydata[0], displaydata[1])
    line.set_3d_properties(displaydata[2])
    return line,


anim = FuncAnimation(fig, animate, frames=int((initTest.orbitalPeriod*1.1) / timestep), interval=20, blit=True)

anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
# plot the orbit
# ax.plot(rlist[0], rlist[1], rlist[2], color="g")

plt.show()
