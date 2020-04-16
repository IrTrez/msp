import numpy as np
from tqdm import tqdm as tqdm
import math
# import matplotlib.pyplot as plt #uncomment when visualisation has been implemented
# import pandas as pd #uncomment when exporting has been implemented.

class Planet:
    def __init__(self, gravitationalParameter, radius, semiMajorAxis, parentGravitationalParameter):
        self.mu = gravitationalParameter
        self.muParent = parentGravitationalParameter
        self.r = radius
        self.a = semiMajorAxis

        self.rsoi = self.a * (self.mu/self.muParent)**0.4
        self.position = [0,0,0]

class Body:
    def __init__(self, parentBody, mass):
        self.mu = parentBody.mu
        self.m = mass

    def initKeplerOrbit(self, semiMajorAxis, eccentricity, timeSincePeriapse):
        self.a = semiMajorAxis
        self.e = eccentricity
        self.tp = timeSincePeriapse

    def refreshKeplerOrbit(self, t):
        self.meanMotion = (self.mu/(self.a**3))
        self.meanAnomaly = self.meanMotion * (t - self.tp)
        if self.e == 1:
            print("Something crazy has happened")
            e = 1 + 10e-10 # ;)
        # True anomaly, normally theta now v
        if self.e < 1:
            E = self.findEccentricAnomaly()
            self.v = math.acos((math.cos(E) - self.e) / (1 - e*math.cos(E)))
            self.fpa = math.asin((self.e*math.sin(E))/math.sqrt(1 - (self.e)**2 * (math.cos(E))**2))
        else:
            H = self.findHyperbolicAnomaly()
            self.v = math.acos((math.cosh(H) - self.e) / (1 - self.e*math.cosh(H)))
            self.fpa = math.sqrt(((self.e)**2 - 1) / ((self.e)**2 * (math.cosh(H))**2) - 1)
        
    def findEccentricAnomaly(self):
        if (self.meanMotion > -1 * math.pi and self.meanMotion < 0 )or self.meanMotion > math.pi:
            E = self.meanAnomaly - self.e
        else:
            E = self.meanAnomaly + self.e

        tolerance = 10e-10
        while True:
            lastE = E
            E += (self.meanAnomaly - E + (self.e*math.sin(E))) / (1 - (self.e * math.cos(E)))
            if abs(lastE - E) < tolerance:
                break
        return E
        
    def findHyperbolicAnomaly(self):
        if self.e < 1.6:
            if (self.meanMotion > -1 * math.pi and self.meanMotion < 0 )or self.meanMotion > math.pi:
                H = self.meanAnomaly - self.e
            else:
                H = self.meanAnomaly + self.e
        else:
            getSign = lambda u: (u > 0) - (u < 0) #lambda function to get sign
            if self.e < 3.6 and abs(self.meanAnomaly) > math.pi:
                H = self.meanAnomaly - math.copysign(e, self.meanAnomaly)  # basically M - sign(M) * e
            else:
                H = self.meanAnomaly/(e-1)

        tolerance = 10e-10
        while True:
            lastE = H
            H += (H + self.meanAnomaly - (self.e*math.sinh(H))) / ((self.e * math.cos(H) - 1))
            if abs(lastE - H) < tolerance:
                break
        return H

