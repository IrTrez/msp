from mpl_toolkits.mplot3d import Axes3D
import main as m
import numpy as np
import matplotlib.pyplot as plt
import time
from tqdm import tqdm
import math

# USE km AS STANDARD DISTANCE UNIT
# USE s AS STANDARD TIME UNIT
AU = 149.6e6  # km
muSun = 1.327178e11
currentTime = time.time()
Mars = m.Planet(42800, 3402, 1.524 * AU, muSun)
Earth = m.Planet(398600.441, 6378.136, AU, muSun)

spacecraft = m.Body(Mars, 100)
# spacecraft.initKeplerOrbit(5000,0.1,currentTime-400)
# spacecraft.refreshKeplerOrbit(currentTime)
# spacecraft.findEccentricAnomaly()

ex26 = m.Body(Earth, 100)

p = 11067.790
e = 0.83285
i = math.radians(87.87)
Omega = math.radians(227.89)
omega = math.radians(53.38)
trueAnomaly = math.radians(92.335)

r, v = ex26.COEtoRV(p, e, i, Omega, omega, trueAnomaly)
print("r: ", r)
print("v: ", v)
