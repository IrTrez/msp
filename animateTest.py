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

p = 10067.790
e = 0.58285
i = math.radians(35)
Omega = math.radians(30.89)
omega = math.radians(53.38)
trueAnomaly = math.radians(40)
a = p / (1 - e**2)

CD = 1.2
surfaceArea = 3.6**2 * math.pi

initTest = m.Body(Earth, 100, CD, surfaceArea)
initTest.initKeplerOrbit(a, e, i, Omega, omega, trueAnomaly)
print(2 * int(initTest.orbitalPeriod/ 200))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.set_aspect('equal')


# Sphere:
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = Earth.r * np.outer(np.cos(u), np.sin(v))
y = Earth.r * np.outer(np.sin(u), np.sin(v))
z = Earth.r * np.outer(np.ones(np.size(u)), np.cos(v))
# Plot the surface
ax.plot_wireframe(x, y, z, color='b')

timestep = 150
steps = 200
rlist = []
for deltat in tqdm(range(steps)):
    initTest.refreshByTimestep(timestep)
    rlist.append(initTest.r)
rlist = np.array(rlist)

line, = ax.plot([], [], lw=2)
ax.set_ylim(-20000, 20000)
ax.set_xlim(-20000, 20000)
ax.set_zlim(-20000, 20000)


# currenttime = 0
def animate(i):
    r = np.sqrt(rlist[i].dot(rlist[i]))
    # ax.annotate("radius =" + str(r), xy=(0,0))
    # ax.text(1, 2, 2, "hey")
    ax.text(10000, 10000, 10000, "radius: {}".format(r))
    displaydata = np.array(rlist[0: i+1]).T
    line.set_data(displaydata[0], displaydata[1])
    line.set_3d_properties(displaydata[2])
    return line,


anim = FuncAnimation(fig, animate, frames= steps, interval=10, blit=True)

anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
# plot the orbit
# ax.plot(rlist[0], rlist[1], rlist[2], color="g")

plt.show()
