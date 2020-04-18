import numpy as np
from tqdm import tqdm as tqdm
import math
from math import acos, cos, cosh, asin, sin,sinh
# import matplotlib.pyplot as plt #uncomment when visualisation has been implemented
# import pandas as pd #uncomment when exporting has been implemented.
print("Fix the gdot")


class Atmosphere:
    def __init__(self, densityFunction = None, densityFile = None):
        assert (densityFunction is None and densityFile is None), "Only one input method of atmosphere can be given"
        assert (densityFunction is not None and densityFile is not None), "Multiple methods of atmosphere are given"
        if densityFunction is not None:
            #create dictionary with density per meter altitude. density in kg/m^3 altitude in km
            pass
        if densityFile is not None:
            # take dictionary from density file
            pass


class Planet:
    def __init__(self, gravitationalParameter, radius, semiMajorAxis, parentGravitationalParameter, atmosphere="false"):
        self.mu = gravitationalParameter
        self.muParent = parentGravitationalParameter
        self.r = radius
        self.a = semiMajorAxis

        self.rsoi = self.a * (self.mu/self.muParent)**0.4
        self.position = [0, 0, 0]

        # Add atmospheric density model here
        if atmosphere != False:
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
            e = 1 + 10e-10  # ;)
        # True anomaly, normally theta now v
        if self.e < 1:
            E = self.findEccentricAnomaly()
            self.trueAnomaly = acos(
                (cos(E) - self.e) / (1 - e*cos(E)))
            self.fpa = asin((self.e*sin(E)) /
                                 math.sqrt(1 - (self.e)**2 * (cos(E))**2))
        else:
            H = self.findHyperbolicAnomaly()
            self.trueAnomaly = acos(
                (cosh(H) - self.e) / (1 - self.e*cosh(H)))
            self.fpa = math.sqrt(((self.e)**2 - 1) /
                                 ((self.e)**2 * (cosh(H))**2) - 1)

    def findEccentricAnomaly(self):
        """Find eccentric anomaly taking the body's current keplerian values
        
        Returns:
            float -- eccentric anomaly [rad]
        """
        if (self.meanMotion > -1 * math.pi and self.meanMotion < 0)or self.meanMotion > math.pi:
            E = self.meanAnomaly - self.e
        else:
            E = self.meanAnomaly + self.e

        tolerance = 10e-10
        while True:
            lastE = E
            E += (self.meanAnomaly - E + (self.e*sin(E))) / \
                (1 - (self.e * cos(E)))
            if abs(lastE - E) < tolerance:
                break
        return E

    def findHyperbolicAnomaly(self):
        """Find hyperbolic anomaly taking the body's current keplerian values
        
        Returns:
            float -- hyperbolic anomaly [rad]
        """
        if self.e < 1.6:
            if (self.meanMotion > -1 * math.pi and self.meanMotion < 0)or self.meanMotion > math.pi:
                H = self.meanAnomaly - self.e
            else:
                H = self.meanAnomaly + self.e
        else:
            # lambda function to get sign
            def getSign(u): return (u > 0) - (u < 0)
            if self.e < 3.6 and abs(self.meanAnomaly) > math.pi:
                # basically M - sign(M) * e
                H = self.meanAnomaly - math.copysign(e, self.meanAnomaly)
            else:
                H = self.meanAnomaly/(e-1)

        tolerance = 10e-10
        while True:
            lastE = H
            H += (H + self.meanAnomaly - (self.e*sinh(H))) / \
                ((self.e * cos(H) - 1))
            if abs(lastE - H) < tolerance:
                break
        return H

    # Vallado Algorithm 1 page 63
    def c2_c3(self, psi):
        """Obtain c2 and c3 values used for the universal solution method of kepler's problem
        
        Arguments:
            psi {float} -- psi value for universal solution method [rad]
        
        Returns:
            float, float -- c2, c3
        """
        if psi > 10e-6:
            c2 = (1 - cos(np.sqrt(psi))) / psi
            c3 = (math.sqrt(psi) - sin(np.sqrt(psi))) / np.sqrt(psi**3)
        elif psi < -10e-6:
            c2 = (1 - cosh(np.sqrt(-psi)))/psi
            c3 = (sinh(np.sqrt(-psi) - np.sqrt(-psi))) / np.sqrt((-psi)**3)
        else:
            c2 = 0.5
            c3 = 1./6.
        return c2, c3

    # Vallado algorithm 8 page 93 change of positions with time
    def keplerTime(self, r, v, dt):
        """Returns new velocity and position vectors with change over time
        
        Arguments:
            r {ndarray} -- Carthesian position vector [km]
            v {ndarray} -- Carthesian velocity vector [km/s]
            dt {float} -- Change in time
        
        Returns:
            ndarray, ndarray -- new velocity and position vectors
        """

        # We call the universal parameter chi (weird greek capital X)
        vnorm = np.sqrt(v.dot(v))
        rnorm = np.sqrt(r.dot(r))

        def getSign(u): return (u > 0) - (u < 0)

        eta = ((vnorm**2) / 2) - (self.mu/rnorm)
        alpha = ((-(vnorm**2)) / self.mu) + (2 / rnorm)
        alphainv = 1/alpha  # this is actually semi major axis

        if self.e < 1:
            assert (alpha != 1), "Alpha in Kepler Time is 1"
            chi0 = np.sqrt(self.mu) * dt * alpha
        else:  # e>1
            chi0 = getSign(dt) * np.sqrt(-alphainv) * np.log((-2 * self.mu * alpha * dt) / (
                np.dot(r, v) + getSign(dt) * math.sqrt(-self.mu * alphainv) * (1-(rnorm * alpha))))

        chi = chi0
        while True:
            chiLast = chi
            psi = chi**2 * alpha
            c2, c3 = self.c2_c3(psi)

            rr = (chi**2 * c2) + ((np.dot(r, v) / np.sqrt(self.mu)) *
                                  chi * (1 - (psi * c3))) + (rnorm * (1 - (psi * c2)))

            chi += ((np.sqrt(self.mu) * dt) - (c3 * chi**3) - ((np.dot(r, v) /
                                                                np.sqrt(self.mu)) * c2 * chi**2) - (rnorm * chi * (1-(psi*c3)))) / rr

            if abs(chiLast-chi) < 10e-6:
                break

        f = 1 - ((chi**2 / rnorm) * c2)
        fdot = (np.sqrt(self.mu) / (rr * rnorm)) * chi * ((psi*c3) - 1)
        g = dt - ((chi**3 / np.sqrt(self.mu)) * c3)
        gdot = 1 - ((chi**2 / rr) * c2)

        # print("test", gdot)

        rnew = (f * r) + (g * v)
        vnew = (fdot * r) + (gdot * v)
        succes = (f*gdot) - (fdot*g)

        return rnew, vnew

    # Vallado algorithm 9
    def RVtoCOE(self, r, v):
        """Converts carthesian position and velocity vectors to keplerian elements
        
        Arguments:
            r {ndarray} -- Carthesian position vector
            v {ndarray} -- Carthesian velocity vector
        
        Returns:
            floats -- semi-latus rectum, semi-major axis, eccentricity, inclination, Omega, omega, True Anomaly, (True omega, u, true lambda) [km] [rad]
        """
        I = [1, 0, 0]
        J = [0, 1, 0]
        K = [0, 0, 1]
        # is this the same as vnorm=np.sqrt(v[0]**2+v[1]**2)
        vnorm = np.sqrt(v.dot(v))
        rnorm = np.sqrt(r.dot(r))

        h = np.cross(r, v)
        hnorm = np.sqrt(h.dot(h))

        n = np.cross(K, h)
        nnorm = np.sqrt(n.dot(n))

        e = (( (vnorm**2-self.mu/rnorm) * r) - (np.dot(r, v) * v)) / \
            (self.mu)  # is a vector
        enorm = np.sqrt(e.dot(e))

        eta = ((vnorm**2/2) - (self.mu/rnorm))

        if enorm != 1.0:
            a = -self.mu / (2 * eta)
            p = a * (1 - enorm**2)
        else:
            p = hnorm**2 / self.mu
            a = math.inf

        i = acos(h[2]/hnorm)

        # acos should work here instead of np.arccos since Omega should not be an array
        Omega = acos(n[0]/nnorm)
        if n[1] < 0:
            Omega = 2*math.pi-Omega

        omega = acos(np.dot(n, e)/(nnorm*enorm))
        if e[2] < 0:
            omega = 2*math.pi-omega

        # nu=greek "v" used for poisson ratio
        trueAnomaly = acos(np.dot(e, r)/(enorm*rnorm))
        if np.dot(r, v) < 0:
            trueAnomaly = 2*math.pi-trueAnomaly
        
        omega_true, lambda_true, u = None, None, None
        if enorm < 1 and i == 0:
            omega_true = acos(e[1]/enorm)
            if e[1] < 0:
                omega_true = 2*math.pi-omega_true
        
        elif enorm == 0 and i != 0:
            u = acos(np.dot(n, r)/(nnorm*rnorm))
            if r[2] < 0:
                u = 2*math.pi-u

        elif enorm == 0 and i == 0:
            lambda_true = acos(r[0]/rnorm)
            if r[1] < 0:
                lambda_true = 2*math.pi-lambda_true

        return p, a, enorm, i, Omega, omega, trueAnomaly, omega_true, u, lambda_true

    # Vallado algorithm 10
    def COEtoRV(self, a, e, i, Omega, omega, trueAnomaly, omega_true=None, u=None, lambda_true=None, p = None):
        """Takes in keplerian elements and outputs carthesian position and velocity vectors
        
        Arguments:
            a {float} -- semi-major axis
            e {float} -- eccentricity
            i {float} -- inclination [rad]
            Omega {float} -- Big Omega [rad]
            omega {float} -- argument of periapse [rad]
            trueAnomaly {float} -- True Anomaly [rad]
        
        Keyword Arguments:
            omega_true {float} -- take omega true as input [rad] (default: {None})
            u {float} -- take omega true as input [rad] (default: {None})
            lambda_true {float} -- take omega true as input [rad] (default: {None})
            p {float} -- Take semi latus rectum as an input instead of semi major axis
        
        Returns:
            r, v -- carthesian position and velocity vectors
        """

        # Semi latus rectum this is different from source material where p is taken as input:
        # when a p is also given it uses p instead of a for calculations
        if p is None:
            p = a (1 - e**2)

        # circular equatorial
        if e < 1 and i == 0:
            assert (lambda_true is not None), "Lambda_True is not given"
            omega, Omega = 0.0, 0.0
        
        # circular inclined
        elif e == 0 and i != 0:
            assert (u is not None), "u is not given"
            omega = 0.0
            trueAnomaly = u

        # elliptical equatorial
        elif e == 0 and i == 0:
            assert (u is not None), "omega_true is not given"
            Omega = 0.0
            omega = omega_true

        rpqw = np.array([((p * cos(trueAnomaly)) / (1 + (e * cos(trueAnomaly)))),
                         ((p * sin(trueAnomaly)) / (1 + (e * cos(trueAnomaly)))),
                         0])

        vpqw = np.array([-1 * np.sqrt(self.mu/p) * sin(trueAnomaly),
                         np.sqrt(self.mu/p) * (e + cos(trueAnomaly)),
                         0])

        PQWtoIJKtransform = np.array([[cos(Omega) * cos(omega) - sin(Omega) * sin(omega) * cos(i), -cos(Omega) * sin(omega) - sin(Omega) * cos(omega) * cos(i), sin(Omega) * sin(i)],
                                        [sin(Omega) * cos(omega) + cos(Omega) * sin(omega) * cos(i), -sin(Omega) * sin(omega) + cos(Omega) * cos(omega) * cos(i), -cos(Omega) * sin(i)],
                                        [sin(omega) * sin(i), cos(omega) * sin(i), cos(i)]])

        rijk = np.matmul(PQWtoIJKtransform, rpqw)
        vijk = np.matmul(PQWtoIJKtransform, vpqw)
        return rijk, vijk