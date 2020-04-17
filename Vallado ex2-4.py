from mpl_toolkits.mplot3d import Axes3D
import main as m
import numpy as np
import matplotlib.pyplot as plt
import time
from tqdm import tqdm

# USE km AS STANDARD DISTANCE UNIT
# USE s AS STANDARD TIME UNIT
AU = 149.6e6 # km
muSun = 1.327178e11
currentTime = time.time()
Mars = m.Planet(42800, 3402, 1.524 * AU, muSun)
Earth = m.Planet(398600.441, 6378.136, AU, muSun)

spacecraft = m.Body(Mars, 100)
# spacecraft.initKeplerOrbit(5000,0.1,currentTime-400)
# spacecraft.refreshKeplerOrbit(currentTime)
# spacecraft.findEccentricAnomaly()

ex24 = m.Body(Earth,100)
ex24.initKeplerOrbit(11000,0.5, currentTime)
r0 = np.array([1131.34, -2282.343, 6672.423])
v0 = np.array([-5.64305, 4.30333, 2.42879])
deltat = 40 * 60

# print(ex24.keplerTime(r0, v0, deltat))
rlist = []
for deltat in tqdm(range(20000)):
    r,_ = ex24.keplerTime(r0, v0, deltat)
    rlist.append(r)
rlist = np.array(rlist).T
# print(rlist.T[1])
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Sphere:
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = Earth.r * np.outer(np.cos(u), np.sin(v))
y = Earth.r * np.outer(np.sin(u), np.sin(v))
z = Earth.r * np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the surface
ax.plot_wireframe(x, y, z, color='r')
ax.set_aspect('equal')

# plot the orbit
ax.plot(rlist[0],rlist[1], rlist[2], color="g")
plt.show()
