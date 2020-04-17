import numpy as np
from tqdm import tqdm as tqdm
import math
# import matplotlib.pyplot as plt #uncomment when visualisation has been implemented
# import pandas as pd #uncomment when exporting has been implemented.
print("Fix the gdot")

class Planet:
    def __init__(self, gravitationalParameter, radius, semiMajorAxis, parentGravitationalParameter, atmosphere="false"):
        self.mu = gravitationalParameter
        self.muParent = parentGravitationalParameter
        self.r = radius
        self.a = semiMajorAxis

        self.rsoi = self.a * (self.mu/self.muParent)**0.4
        self.position = [0,0,0]

        # Add atmospheric density model here
        if atmosphere:
            pass

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
            self.trueAnomaly = math.acos((math.cos(E) - self.e) / (1 - e*math.cos(E)))
            self.fpa = math.asin((self.e*math.sin(E))/math.sqrt(1 - (self.e)**2 * (math.cos(E))**2))
        else:
            H = self.findHyperbolicAnomaly()
            self.trueAnomaly = math.acos((math.cosh(H) - self.e) / (1 - self.e*math.cosh(H)))
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


    # definition of psi Vallado 2-48
    # def psi(self, UniversalVariable):
    #     return UniversalVariable**2 / self.a
    # Vallado Algorithm 1 page 63
    def c2_c3(self, psi):
        if psi > 10e-6:
            c2 = (1 - math.cos(np.sqrt(psi))) / psi
            c3 =  (math.sqrt(psi) - math.sin(np.sqrt(psi))) / np.sqrt(psi**3)
        elif psi < -10e-6:
            c2 = (1 - math.cosh(np.sqrt(-psi)))/psi
            c3 = (math.sinh(np.sqrt(-psi) - np.sqrt(-psi))) / np.sqrt((-psi)**3)
        else:
            c2 = 0.5
            c3 = 1./6.
        return c2, c3

    # Vallado algorithm 8 page 93 change of positions with time
    def keplerTime(self, r, v, dt):
        # r: numpy array carthesian position
        # v: numpy array velocities
        # We call the universal parameter chi (weird greek capital X)
        vnorm = np.sqrt(v.dot(v))
        rnorm = np.sqrt(r.dot(r))
        def getSign(u): return (u > 0) - (u < 0)

        eta = ((vnorm**2) / 2) - (self.mu/rnorm)
        alpha = ( (-(vnorm**2)) / self.mu) + (2 / rnorm)
        alphainv = 1/alpha #this is actually semi major axis

        if self.e < 1:
            assert (alpha != 1), "Alpha in Kepler Time is 1"
            chi0 = np.sqrt(self.mu) * dt * alpha
        else: # e>1
            chi0 = getSign(dt) * np.sqrt(-alphainv) * np.log((-2 * self.mu * alpha * dt) / (np.dot(r,v) + getSign(dt) * math.sqrt(-self.mu * alphainv) * (1-(rnorm * alpha)) ))

        
        chi = chi0
        while True:
            chiLast = chi
            psi = chi**2 * alpha
            c2, c3 = self.c2_c3(psi)

            rr = (chi**2 * c2) + ((np.dot(r,v) / np.sqrt(self.mu)) * chi * (1 - (psi * c3)) ) + (rnorm * (1 - (psi * c2) ))

            chi += ( (np.sqrt(self.mu) * dt) - (c3 * chi**3) - ((np.dot(r,v) / np.sqrt(self.mu)) * c2 * chi**2) - (rnorm * chi *(1-(psi*c3)))) / rr

            if abs(chiLast-chi) < 10e-6:
                break
        
        f = 1 - ((chi**2 / rnorm) * c2)
        fdot = (np.sqrt(self.mu) / (rr * rnorm) ) * chi * ( (psi*c3) - 1)
        g = dt - ( (chi**3 / np.sqrt(self.mu) ) * c3)
        gdot = 1 - ((chi**2 / rr) * c2)

        # print("test", gdot)
        
        rnew = (f * r) + (g * v)
        vnew = (fdot * r) + (gdot * v)

        succes = (f*gdot) - (fdot*g)
        # print("1 =", succes)
        return rnew, vnew

    # Vallado algorithm 9
    def RVtoCOE(self, r,v):
        I=[1,0,0]
        J=[0,1,0]
        K=[0,0,1]
        vnorm = np.sqrt(v.dot(v)) #is this the same as vnorm=np.sqrt(v[0]**2+v[1]**2)
        rnorm = np.sqrt(r.dot(r))

        h=np.cross(r,v)
        hnorm=np.sqrt(h.doth(h))

        n=np.cross(K,h)
        nnorm=np.sqrt(n.dot(n))

        e=((vnorm2-self.mu/rnorm)*r-np.dot(r,v)*v)/(self.mu) #is a vector
        enorm=np.sqrt(e.dot(e))
        
        if enorm == 0:
            p=hnorm**2/self.mu
            a=math.inf
        else:
            a=-self.mu/(2*(vnorm**2/2-self.mu/rnorm))

        i=arccos(h[2]/hnorm)

        Omega=acos(n[0]/nnorm) #acos should work here instead of np.arccos since Omega should not be an array
        if n[1] < 0:
            Omega=2*Pi()-Omega

        omega=acos(np.dot(n,e)/(nnorm*enorm))
        if e[2] < 0:
            omega=2*Pi()-omega

       nu=acos(np.dot(e,r)/(enorm*rnorm)) #nu=greek "v" used for poisson ratio
        if np.dot(r,v) < 0:
            nu=2*Pi()-nu
        
        return(p,a,e,i,Omega, omega,nu) 
        
        