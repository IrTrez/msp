import main as m
import simtools
import time
from tqdm import tqdm
import math

DATAFILE = "runs/ManouvreTest.csv"
ATMOSPHEREDATA = "densityModels/MarsDensity.csv"
SPEED = 600  # __ times speed

# USE km AS STANDARD DISTANCE UNIT
# USE s AS STANDARD TIME UNIT

currentTime = time.time()
limitAltitude = 500 # 260  #[km]. At this altitude density is just below 1*10^-10

Mars_atmosphere=m.Atmosphere(limitAltitude, densityFile=ATMOSPHEREDATA)
Earth = simtools.Earth

a = 10000
e = 0.3
i = math.radians(1)
Omega = math.radians(30.89)
omega = math.radians(60.0)
trueAnomaly = math.radians(40)

CD = 1.2
surfaceArea = 3.6**2 * math.pi

spacecraft = m.Body(Earth, 100, CD, surfaceArea )
spacecraft.initKeplerOrbit(a,e,i,Omega,omega,trueAnomaly)


# Add manoeuvers before propagate

# spacecraft.addManouvreByDirection(spacecraft.start + 1 * spacecraft.orbitalPeriod, 6, "n")
# spacecraft.addManouvreByDirection(spacecraft.start + 1 * spacecraft.orbitalPeriod, -6, "t")
# spacecraft.addManouvreByDirection(spacecraft.start + 2 * spacecraft.orbitalPeriod, 1.5, "r")
# spacecraft.addManouvreByDirection("p0", 2, "n")
# spacecraft.addManouvreByDirection("p1", 1, "r")


# PROPAGATE Here
rlist = spacecraft.propagate(1*spacecraft.orbitalPeriod, DATAFILE, False)

simtools.quickAnimate(SPEED, DATAFILE, Earth)
