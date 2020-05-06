import numpy as np
from tqdm import tqdm as tqdm
import math
from math import acos, cos, cosh, asin, sin,sinh
import time
from math import acos, cos, cosh, asin, sin,sinh, exp
# import matplotlib.pyplot as plt #uncomment when visualisation has been implemented
import pandas as pd


class Atmosphere:
    def __init__(self, limitAltitude, densityFunction = None, densityFile = None):
        self.limitAltitude = limitAltitude
        assert not (densityFunction is not None and densityFile is not None), "Only one input method of atmosphere can be given"
        assert not (densityFunction is None and densityFile is None), "Multiple methods of atmosphere are given"
        if densityFunction is not None:
            # return the function
            # or
            # create dictionary with density per meter altitude. density in kg/m^3 altitude in km
            print("DensityFunction is temporarily deprecated, use densityFile instead")
            self.densityFunction = densityFunction
        if densityFile is not None:
            density = pd.read_csv(densityFile)
            density.iloc[:, 0] = density.iloc[:, 0].apply(round,args=[3])
            density.set_index("Altitude", inplace=True)
            self.densityDict = density.to_dict()["Density"]



class Planet:
    def __init__(self, gravitationalParameter, radius, semiMajorAxis, parentGravitationalParameter, atmosphere:Atmosphere="false"):
        self.mu = gravitationalParameter
        self.muParent = parentGravitationalParameter
        self.r = radius
        self.a = semiMajorAxis

        self.rsoi = self.a * (self.mu/self.muParent)**0.4
        self.position = [0, 0, 0]

        # Add atmospheric density model here
        if atmosphere != False:
            self.atmosphere = atmosphere
            self.atmosphericLimitAltitude = atmosphere.limitAltitude
            self.densityDict = atmosphere.densityDict


class Body:
    def __init__(self, parentBody, mass, DragCoeff=0, surfaceArea=0):
        self.mu = parentBody.mu
        self.parentRadius = parentBody.r
        self.m = mass
        self.CD = DragCoeff
        self.surfaceArea = surfaceArea
        self.clock = time.time()
        self.start = self.clock
        self.parentRSOI = parentBody.rsoi
        self.manoeuvers = []
        if parentBody.atmosphere != False:
            self.atmosphericLimitAltitude = parentBody.atmosphericLimitAltitude
            self.densityDict = parentBody.densityDict

    def getDensity(self, altitude):
        return self.densityDict[round(altitude,3)]


    def initKeplerOrbit(self, semiMajorAxis, eccentricity, inclination, Omega, omega, trueAnomaly = 0.0, useDegrees = False):
        
        """
        if useDegrees:
            i = math.radians(i)
            Omega = math.radians(Omega)
            omega = math.radians(omega)
            trueAnomaly = math.radians(trueAnomaly)
        """

        self.a = semiMajorAxis
        self.e = eccentricity
        self.i = inclination
        self.Omega = Omega
        self.omega = omega
        self.trueAnomaly = trueAnomaly

        self.r, self.v = self.COEtoRV(self.a, self.e, self.i, self.Omega, self.omega, self.trueAnomaly)
        if self.e < 1:
            self.orbitalPeriod = 2 * math.pi * math.sqrt(self.a**3/self.mu)
        self.altitude = self.r - (self.parentRadius * (self.r/np.sqrt(self.r.dot(self.r))))
        self.apoapsis = self.a * (1 + self.e)
        self.periapsis = self.a * (1 - self.e)
        self.counter = 0        #only used to force add a manoeuver
        

        self.dt = 1 # default global timestep


        # self.tp = timeSincePeriapse
    def initPositionOrbit(self, r, v):
        """Set up an orbit using carthesian position and velocity vectors
        
        Arguments:
            r {numpy.ndarray} -- carthesian position vector
            v {numpy.ndarray} -- carthesian velocity vector
        """
        self.r = r
        self.v = v
        self.altitude = self.r - (self.parentRadius * (self.r/np.sqrt(r.dot(r))))

        self.p, self.a, self.e, self.i, self.Omega, self.omega, self.trueAnomaly, _, _, _ = self.RVtoCOE(self.r, self.v)
        if self.e < 1:
            self.orbitalPeriod = 2 * math.pi * math.sqrt(self.a**3/self.mu)
        self.altitude = self.r - (self.parentRadius * (self.r/np.sqrt(self.r.dot(self.r))))
        self.apoapsis = self.a * (1 + self.e)
        self.periapsis = self.a * (1 - self.e)
        self.counter = 0        #only used to force add a manoeuver
        

        self.dt = 1 # default global timestep


    def refreshByTimestep(self, dt, atmospheric):
        self.clock += dt
        rnew, vnew = self.keplerTime(self.r, self.v, dt)
        self.altitude = rnew - (self.parentRadius * (rnew/np.sqrt(rnew.dot(rnew))))
        altitudenorm = np.sqrt(self.altitude.dot(self.altitude))
        if (altitudenorm < self.atmosphericLimitAltitude and atmospheric):
            vnorm = np.sqrt(vnew.dot(vnew))

            adrag = -0.5 * self.getDensity(altitudenorm) * ((self.CD * self.surfaceArea)/self.m) * vnew**2 * (vnew/vnorm)
            vnew += adrag * dt
        
        
        self.initPositionOrbit(rnew, self.Manoeuvers(vnew))

    def propagate(self, timeJump, saveFile = None, atmospheric = False, dtAtmospheric = 1, dtNormal = 1):
        rlist = []
        clocklist = []
        for deltat in tqdm(range((int(timeJump / abs(dtAtmospheric))) + 1)):
            if np.sqrt(self.r.dot(self.r)) + 1000 > self.parentRSOI:
                # sphere of influence check
                print("Body has left sphere of influence")
                print("Radius: ", self.r)
                print("Velocity: ", self.v)
                break
            if self.parentRadius > np.sqrt(self.r.dot(self.r)):
                # Check if body has crashed
                print("Body crashed into surface")
                print("Ending propagation")
                # print("Simulation ran for: " + str(time.time() - self.start))
                print("Radius:", np.sqrt(self.r.dot(self.r)), ":", self.r)
                print("Altitude:", self.altitude.dot(self.altitude))
                break
            #plus 500 is a safety margin as it's takeing the previous altitude
            if np.sqrt(self.altitude.dot(self.altitude)) < self.atmosphericLimitAltitude + 500: 
                self.refreshByTimestep(dtAtmospheric , atmospheric)
                rlist.append(self.r)
                clocklist.append(self.clock)
            else:
                self.refreshByTimestep(dtNormal, atmospheric)
                rlist.append(self.r)
                clocklist.append(self.clock)
        if saveFile is not None:
            rlist = np.array(rlist).T
            data = pd.DataFrame({"clock":clocklist, "x":rlist[0], "y":rlist[1], "z":rlist[2]})
            data.to_csv(saveFile, sep=",")

            manoeuvreSavefile = saveFile[:-4] + "_man" + ".csv"
            manoeuvreData = pd.DataFrame(self.manoeuvers)
            manoeuvreData.to_csv(manoeuvreSavefile)
        return rlist


    # def refreshKeplerOrbit(self, t):
    #     # might be deprecated
    #     print("refreshKeplerOrbit is deprecated")
    #     self.meanMotion = (self.mu/(self.a**3))
    #     self.meanAnomaly = self.meanMotion * (t - self.tp)
    #     if self.e == 1:
    #         print("Something crazy has happened")
    #         e = 1 + 10e-10  # ;)
    #     # True anomaly, normally theta now v
    #     if self.e < 1:
    #         E = self.findEccentricAnomaly()
    #         self.trueAnomaly = acos(
    #             (cos(E) - self.e) / (1 - e*cos(E)))
    #         self.fpa = asin((self.e*sin(E)) /
    #                              math.sqrt(1 - (self.e)**2 * (cos(E))**2))
    #     else:
    #         H = self.findHyperbolicAnomaly()
    #         self.trueAnomaly = acos(
    #             (cosh(H) - self.e) / (1 - self.e*cosh(H)))
    #         self.fpa = math.sqrt(((self.e)**2 - 1) /
    #                              ((self.e)**2 * (cosh(H))**2) - 1)


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
        sqrtmu = self.mu**0.5
        rdotv = r.dot(v)

        eta = ((vnorm**2) / 2) - (self.mu/rnorm)
        alpha = ((-(vnorm**2)) / self.mu) + (2 / rnorm)

        if alpha > 0:
            assert (alpha != 1), "Alpha in Kepler Time is 1"
            chi = sqrtmu * dt * alpha
        elif alpha < 0:  # e>1
            alphainv = 1/alpha  # this is actually semi major axis
            chi = (np.sign(dt) * (-alphainv)**.5 *
                      np.log((-2 * self.mu * alpha * dt) /
                             (rdotv + np.sign(dt) *
                              np.sqrt(-self.mu / alpha) * (1 - rnorm * alpha))))
        else:
            print("parabola in keplerTime")
            chi = sqrtmu * dt / rnorm


        while True:
            chiLast = chi
            psi = chi * chi * alpha
            c2, c3 = self.c2_c3(psi)

            rr = (chi * chi * c2) + ((rdotv/sqrtmu) * chi * (1-psi*c3)) + (rnorm * (1-psi*c2))

            chi += ((np.sqrt(self.mu) * dt) - (chi*chi*chi*c3) - ((np.dot(r, v)/self.mu)*chi*chi*c2) - (rnorm*chi*(1-psi*c3)))/rr

            if abs(chiLast-chi) < 10e-6:
                break


        f = 1 - ((chi*chi / rnorm) * c2)
        fdot = (np.sqrt(self.mu) / (rr * rnorm)) * chi * ((psi*c3) - 1)
        g = dt - ((chi*chi*chi / np.sqrt(self.mu)) * c3)
        gdot = 1 - ((chi*chi / rr) * c2)

        rnew = (f * r) + (g * v)
        vnew = (fdot * r) + (gdot * v)
        succes = (f * gdot) - (fdot * g)
        assert abs(succes-1) < 10e-6, "Succes is not 1 but: " + str(succes)
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

        e = (( (vnorm**2-self.mu/rnorm) * r) - (np.dot(r, v) * v)) / (self.mu)  # is a vector
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
            p = a * (1 - e**2)

        # # circular equatorial
        # if e < 1 and i == 0:
        #     assert (lambda_true is not None), "Lambda_True is not given"
        #     omega, Omega = 0.0, 0.0
        
        # # circular inclined
        # elif e == 0 and i != 0:
        #     assert (u is not None), "u is not given"
        #     omega = 0.0
        #     trueAnomaly = u

        # # elliptical equatorial
        # elif e == 0 and i == 0:
        #     assert (u is not None), "omega_true is not given"
        #     Omega = 0.0
        #     omega = omega_true

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


    def ApoapsisCheck(self):
        """
        Returns:
            True/False statment depending on position of space-craft
        """

        self.checkA = False

        if math.isclose(np.sqrt(self.r.dot(self.r)), self.apoapsis, rel_tol=1e-4):
            self.checkA = True

        return self.checkA


    def PeriapsisCheck(self):
        """
        Returns:
            True/False statment depending on position of space-craft
        """

        self.checkP = False

        if math.isclose(np.sqrt(self.r.dot(self.r)), self.apoapsis, rel_tol=1e-4):
            self.checkP = True

        return self.checkP


    def ParkingOrbit(self):
        '''

            ap {float} -- semi major axis for parking orbit [km]

        Returns:
            POdv {ndarray} -- Carthesian velocity vector for transfer at apoapsis to parking orbit
        '''

        self.ap = 0.5 * (self.apoapsis+self.periapsis+self.atmosphericLimitAltitude)

        self.POdv = abs(math.sqrt(self.mu*(2/self.apoapsis-1/self.a))-math.sqrt(self.mu*(2/self.apoapsis-1/self.ap)))
        self.POdv = (self.v/np.sqrt(self.v.dot(self.v)))*self.POdv

        return self.POdv


    def Manoeuvers(self,v):
        '''

        parameter:
            v (ndarray) -- vnew, the most recent calculations for current velocity in cartesian coordinates

        '''
        if len(self.manoeuvers) != 0:
            for manoeuver in self.manoeuvers:
                if self.clock > manoeuver["clock"] and manoeuver["expired"] == False:
                    v += manoeuver["dv"]
                    self.manoeuvers[manoeuver["ID"]]["r"] = self.r
                    manoeuver["expired"] = True
        return v


    def addManoeuverByVector(self, clock, dv):
        '''

            clock (float) -- clock of manoeuvers
            dv (ndarray) -- delta v for maneuvers in cartesian coordinates

        '''
        # gets latest ID
        if len(self.manoeuvers) == 0:
            latestID = 0
        else:
            latestID = len(self.manoeuvers)

        self.manoeuvers.append({"ID": latestID, "clock": clock, "dv": dv, "expired": False, "direction": (
            dv/np.sqrt(dv.dot(dv))), "manType": "g"})

    def addManouvreByDirection(self, clock, dvMagnitude, manoeuvreType):
        """[summary]

        Arguments:
            time {float} -- time
            dv {float} -- normal value of deltav
            type {string} -- type(tangential, radial and normal) or (t,r,n)
        """
        manoeuvreTypes = ["tangential", "t", "radial", "r", "normal", "n"]
        assert (manoeuvreType in manoeuvreTypes), "Incorrect manouvreType, use tangential, t, radial, r, normal or n"

        # gets latest ID
        if len(self.manoeuvers) == 0:
            latestID = 0
        else:
            latestID = len(self.manoeuvers)

        # Create dv vector based on manoeuvretype and dvMagnitude
        if (manoeuvreType == "tangential" or manoeuvreType == "t"):
            manType = "t"
            direction = self.v/np.sqrt(self.v.dot(self.v))
            dv = dvMagnitude * direction

        if (manoeuvreType == "radial" or manoeuvreType == "r"):
            manType = "r"
            direction = -self.r/np.sqrt(self.r.dot(self.r))
            dv = dvMagnitude * direction

        if (manoeuvreType == "normal" or manoeuvreType == "n"):
            manType = "n"
            directionTangential = self.v/np.sqrt(self.v.dot(self.v))
            directionRadial = -self.r/np.sqrt(self.r.dot(self.r))
            directionNormal = np.cross(directionTangential,directionRadial)
            direction = directionNormal/np.sqrt(directionNormal.dot(directionNormal))
            dv = dvMagnitude * direction

        self.manoeuvers.append({"ID": latestID, "clock": clock, "dv": dv, "expired": False, "direction": (
            dv/np.sqrt(dv.dot(dv))), "manType": manType})
