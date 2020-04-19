from mpl_toolkits.mplot3d import Axes3D
import main as m
import numpy as np
import matplotlib.pyplot as plt
import time
from tqdm import tqdm
import math
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
print(initTest.orbitalPeriod)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Sphere:
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = Earth.r * np.outer(np.cos(u), np.sin(v))
y = Earth.r * np.outer(np.sin(u), np.sin(v))
z = Earth.r * np.outer(np.ones(np.size(u)), np.cos(v))
# Plot the surface
ax.plot_wireframe(x, y, z, color='b')

timestep = 1000
rlist = []
hlist = []
parentradiuslist = []
for deltat in tqdm(range(int(initTest.orbitalPeriod / timestep) + 1)):
    initTest.refreshByTimestep(timestep)
    rlist.append(initTest.r)
    parentradius = (initTest.parentRadius * (initTest.r/np.sqrt((initTest.r).dot(initTest.r))))
    h = initTest.r - parentradius
    parentradiuslist.append(parentradius)
    hlist.append(h)
    # print(np.sqrt(parentradius.dot(parentradius)))

hlist = np.array(hlist)
zerolist = np.zeros(hlist.shape)
rlist = np.array(rlist).T
parentradiuslist = np.array(parentradiuslist)

ax.set_ylim(-30000, 30000)
ax.set_xlim(-30000, 30000)
ax.set_zlim(-30000, 30000)

# print(rlist)
plt.pause(1)
for u in tqdm(range(len(rlist[0]))):
    # print(deltat)
    ax.plot(rlist[0][0:u], rlist[1][0:u], rlist[2][0:u], color="g")
    ax.quiver(*parentradiuslist[u], *hlist[u], color="red")
    plt.pause(0.000001)
    plt.clf
plt.pause(5)
# rlist = np.array(rlist).T

# plot the orbit

# ax.set_aspect('equal')
# plt.show()
