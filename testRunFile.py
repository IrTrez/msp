import main as m
import numpy as np
import matplotlib.pyplot as plt
import time

# USE km AS STANDARD DISTANCE UNIT
# USE s AS STANDARD TIME UNIT
AU = 149.6e6 # km
currentTime = time.time()
Mars = m.Planet(42800,3402, 1.524 * AU, 1.327178e11)

spacecraft = m.Body(Mars,100)
spacecraft.initKeplerOrbit(5000,0.1,currentTime-400)
spacecraft.refreshKeplerOrbit(currentTime)
# spacecraft.findEccentricAnomaly()