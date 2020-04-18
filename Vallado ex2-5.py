from mpl_toolkits.mplot3d import Axes3D
import main as m
import numpy as np
import matplotlib.pyplot as plt
import time
from tqdm import tqdm

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

ex25 = m.Body(Earth, 100)
# ex25.initKeplerOrbit(11000,0.5, currentTime)
r0 = np.array([6524.834, 6862.875, 6448.296])
v0 = np.array([4.901327, 5.533756, -1.976341])
deltat = 40 * 60

print(ex25.RVtoCOE(r0, v0))
