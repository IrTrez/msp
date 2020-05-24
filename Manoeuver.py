import main as m
import simtools
import time
from tqdm import tqdm
import math
import numpy as np

DATAFILE = "runs/ManouvreTest.csv"
ATMOSPHEREDATA = "densityModels/MarsDensity.csv"
SPEED = 600  # __ times speed

# USE km AS STANDARD DISTANCE UNIT
# USE s AS STANDARD TIME UNIT

currentTime = time.time()
Earth = simtools.Earth

a = 10000
e = 0.3
i = math.radians(1)
Omega = math.radians(30.0)
omega = math.radians(60.0)
trueAnomaly = math.radians(40)

CD = 1.2
surfaceArea = 3.6**2 * math.pi

spacecraft = m.Body(Earth, 100, CD, surfaceArea )
spacecraft.initKeplerOrbit(a,e,i,Omega,omega,trueAnomaly)


# Add manoeuvers before propagate
# Example of an inclination change
# Note that for best results add Manoeuvres in order.
spacecraft.addManouvreByDirection("p0", 1.3, "t")
spacecraft.addManouvreByDirection("p2", -1.3, "t")
spacecraft.addManouvreByDirection("a1", 0.3, "n")
spacecraft.addManouvreByDirection("a2", 0.3, "n")

# Some more examples of manoeuvres
# spacecraft.addManouvreByDirection(1000, -0.5, "n")
# spacecraft.addManouvreByDirection(spacecraft.start + 1 * spacecraft.orbitalPeriod, 6, "n")
# spacecraft.addManouvreByDirection("p1", 1, "r")\
# spacecraft.addManoeuverByVector(spacecraft.start + 1000, np.array([1,-1,1]))

# PROPAGATE Here
spacecraft.propagate(12*spacecraft.orbitalPeriod, DATAFILE, False)

simtools.quickAnimate(SPEED, DATAFILE, Earth)
