import time
import numpy as np
import pandas as pd
from tqdm import tqdm as tqdm
from math import acos, cos, cosh, asin, sin, sinh, exp, acosh, asinh, copysign,\
     sqrt, pi, inf, radians


class Atmosphere:
    def __init__(self, limitAltitude, densityFunction = None, densityFile = None):
        self.limitAltitude = limitAltitude
        assert not (densityFunction is not None and densityFile is not None), "Only one input method of atmosphere can be given"
        assert not (densityFunction is None and densityFile is None), "Multiple methods of atmosphere are given"
        if densityFunction is not None:
            # create dictionary with density per meter altitude. density in kg/m^3 altitude in km
            raise DeprecationWarning("DensityFunction is temporarily deprecated, use densityFile instead")
        if densityFile is not None:
            density = pd.read_csv(densityFile)
            density.iloc[:, 0] = density.iloc[:, 0].apply(round,args=[3])
            density.set_index("Altitude", inplace=True)
            self.densityDict = density.to_dict()["Density"]



class Planet:
    def __init__(self, gravitationalParameter, radius, semiMajorAxis, parentGravitationalParameter, atmosphere:Atmosphere=False):
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
        # Parentbody variables
        self.mu = parentBody.mu
        self.parentRadius = parentBody.r
        self.parentRSOI = parentBody.rsoi
        if parentBody.atmosphere != False:
            self.atmosphericLimitAltitude = parentBody.atmosphericLimitAltitude
            self.densityDict = parentBody.densityDict

        # Body variables
        self.m = mass
        self.CD = DragCoeff
        self.surfaceArea = surfaceArea

        # System variables
        self.clock = time.time()
        self.start = self.clock
        self.manoeuvers = []


    def getDensity(self, altitude):
        """Pulls density from self.densityDict

        Arguments:
            altitude {float} -- altitude

        Returns:
            float -- Atmospheric density in kg/m^3
        """
        return self.densityDict[round(altitude,3)]


    def initKeplerOrbit(self, semiMajorAxis, eccentricity, inclination, Omega, omega, trueAnomaly = 0.0, useDegrees = False):
        """Set up or refresh an orbit's parameters using keplerian elements
        
        Arguments:
            semiMajorAxis {float} -- Semi-Major axis [km]
            eccentricity {float} -- eccentricity
            inclination {float} -- inclination [rad]
            Omega {float} -- Big omega [rad]
            omega {float} -- argument of periapse [rad]
        
        Keyword Arguments:
            trueAnomaly {float} -- True Anomaly [rad] (default: {0})
            useDegrees {boolean} -- Option to init angles with degrees  (default: {False})
        """
        if useDegrees:
            i = radians(i)
            Omega = radians(Omega)
            omega = radians(omega)
            trueAnomaly = radians(trueAnomaly)

        self.a = semiMajorAxis
        self.e = eccentricity
        self.i = inclination
        self.Omega = Omega
        self.omega = omega
        self.trueAnomaly = trueAnomaly

     # anomalies and time to Periapse
        if self.e < 1:
            self.orbitalPeriod = 2 * pi * sqrt((self.a * self.a * self.a)/self.mu)

            eccentricAnomaly = asin((sin(self.trueAnomaly) * sqrt(1 - self.e*self.e)) / (1 + (self.e * cos(self.trueAnomaly))))
            meanAnomaly = eccentricAnomaly - (self.e * sin(eccentricAnomaly))
            n = sqrt(self.mu / (self.a * self.a * self.a))
            self.timeToPeriapsis = meanAnomaly / n
            self.timeToApoapsis = pi - meanAnomaly / n

        elif self.e > 1:
            # hyperbolicAnomaly = asinh((sin(self.trueAnomaly) * sqrt(self.e*self.e - 1)) / (1 + (self.e * cos(self.trueAnomaly))))
            hyperbolicAnomaly = acosh((self.e + cos(self.trueAnomaly) / (1 + self.e * cos(self.trueAnomaly))))
            meanAnomaly = self.e * sinh(hyperbolicAnomaly) - hyperbolicAnomaly
            n = sqrt(self.mu/(- self.a * self.a * self.a))
            self.timeToPeriapsis = -meanAnomaly / n

        self.r, self.v = self.COEtoRV(self.a, self.e, self.i, self.Omega, self.omega, self.trueAnomaly)
   
        self.altitude = self.r - (self.parentRadius * (self.r/np.sqrt(self.r.dot(self.r))))
        self.apoapsis = self.a * (1 + self.e)
        self.periapsis = self.a * (1 - self.e)

        self.dt = 1 # default global timestep


    def initPositionOrbit(self, r, v):
        """Set up or refresh an orbit using carthesian position and velocity vectors
        
        Arguments:
            r {numpy.ndarray} -- carthesian position vector
            v {numpy.ndarray} -- carthesian velocity vector
        """
        self.r = r
        self.v = v
        self.altitude = self.r - (self.parentRadius * (self.r/np.sqrt(r.dot(r))))

        self.p, self.a, self.e, self.i, self.Omega, self.omega, self.trueAnomaly, _, _, _ = self.RVtoCOE(self.r, self.v)

                # anomalies and time to Periapse
        if self.e < 1:
            self.orbitalPeriod = 2 * pi * sqrt((self.a * self.a * self.a)/self.mu)

            eccentricAnomaly = asin((sin(self.trueAnomaly) * sqrt(1 - self.e*self.e)) / (1 + (self.e * cos(self.trueAnomaly))))
            meanAnomaly = eccentricAnomaly - (self.e * sin(eccentricAnomaly))
            n = sqrt(self.mu / (self.a * self.a * self.a))
            self.timeToPeriapsis = meanAnomaly / n
            self.timeToApoapsis = pi - meanAnomaly / n

        elif self.e > 1:
            # hyperbolicAnomaly = asinh((sin(self.trueAnomaly) * sqrt(self.e*self.e - 1)) / (1 + (self.e * cos(self.trueAnomaly))))
            hyperbolicAnomaly = acosh((self.e + cos(self.trueAnomaly) / (1 + self.e * cos(self.trueAnomaly))))
            meanAnomaly = self.e * sinh(hyperbolicAnomaly) - hyperbolicAnomaly
            n = sqrt(self.mu/(- self.a * self.a * self.a))
            self.timeToPeriapsis = -meanAnomaly / n

        self.altitude = self.r - (self.parentRadius * (self.r/np.sqrt(self.r.dot(self.r))))
        self.apoapsis = self.a * (1 + self.e)
        self.periapsis = self.a * (1 - self.e)

        self.dt = 1 # default global timestep


    def refreshByTimestep(self, dt, atmospheric):
        """Medium level function to refresh the elements of an orbit. \n
        Adds any acceleration and pertubation effects

        Arguments:\n
            dt {float} -- Time step
            atmospheric {bool} -- If the prerequisite parentbody parameters are available perform atmospheric deceleration
        """
        self.clock += dt
        rnew, vnew = self.keplerTime(self.r, self.v, dt)
        self.altitude = rnew - (self.parentRadius * (rnew/np.sqrt(rnew.dot(rnew))))
        altitudenorm = np.sqrt(self.altitude.dot(self.altitude))
        if (altitudenorm < self.atmosphericLimitAltitude and atmospheric):
            vnorm = np.sqrt(vnew.dot(vnew))
            adrag = -0.5 * self.getDensity(altitudenorm) * ((self.CD * self.surfaceArea)/self.m) * vnew * vnew * (vnew/vnorm)
            vnew += adrag * dt
        
        # Add manoeuvre delta V
        if len(self.manoeuvers) != 0 and any(not x["expired"] for x in self.manoeuvers):
            Adistancer = 0
            Pdistancer = 0
            for manoeuver in self.manoeuvers:
                # this is so that it doesnt propagate multiple manouevres close to ap/peri
                manoeuver["iterSinceCreation"] += 1
                
                rnorm = np.sqrt(rnew.dot(rnew))
                atPeriapsis = (abs(self.periapsis-rnorm) < 0.001)
                atApoapsis = (abs(self.apoapsis-rnorm) < 0.001)

                if (isinstance(manoeuver["clock"], str)) and (self.clock-self.start > 2) and (atApoapsis or atPeriapsis):
                    if (manoeuver["clock"][0] in ["p", "a"]) and (manoeuver["expired"] == False) and (manoeuver["iterSinceCreation"] > 10):
                        location = manoeuver["clock"][0]
                        orbit = int(manoeuver["clock"][1:])
                        
                        # This is to make sure it only propagates the manoeuver once at peri/apoapsis
                        if Adistancer > 0:
                            Adistancer -= 1
                            continue
                        if Pdistancer > 0:
                            Pdistancer -= 1
                            continue

                        if atPeriapsis:
                            Pdistancer = 10
                        if atPeriapsis:
                            Adistancer = 10

                        if location == "p" and atPeriapsis:
                            # print(manoeuver)
                            if orbit == 0:
                                self.manoeuvers.remove(manoeuver)
                                self.addManouvreByDirection(self.clock, manoeuver["dv"], manoeuver["manType"], 0, False)
                            else:
                                orbit -= 1
                                manoeuver["clock"] = location + str(orbit)
                                manoeuver["iterSinceCreation"] = 0

                        elif location == "a" and atApoapsis:
                            # print(manoeuver)
                            if orbit == 0:
                                self.manoeuvers.remove(manoeuver)
                                self.addManouvreByDirection(self.clock, manoeuver["dv"], manoeuver["manType"], 0, False)
                            else:
                                orbit -= 1
                                manoeuver["clock"] = location + str(orbit)
                                manoeuver["iterSinceCreation"] = 0 
                    else:
                        assert (manoeuver["clock"][0] in ["p", "a"]), "Only p and a supported."

                elif isinstance(manoeuver["clock"], int) or isinstance(manoeuver["clock"], float):
                    if self.clock > manoeuver["clock"] and manoeuver["expired"] == False:
                        if isinstance(manoeuver["direction"], str):
                            self.manoeuvers.remove(manoeuver)
                            self.addManoeuverByDirection(self.clock, manoeuver["dv"], manoeuver["manType"], 0, False)
                        else:
                            vnew += manoeuver["dv"]
                            manoeuver["r"] = self.r
                            manoeuver["expired"] = True

        self.initPositionOrbit(rnew, vnew)

    def propagate(self, timeJump, saveFile = None, atmospheric = False, dtAtmospheric = 1, dtNormal = 1):
        """Highest level function to propagate a body's orbit with time.

        Arguments:
            timeJump {int} -- The total time in seconds to propage a spacecraft's orbit

        Keyword Arguments:
            saveFile {string} -- path for CSV savefile, if None: no save is made (default: {None})
            atmospheric {bool} -- Indicate if atmospheric effects should be accounted for (default: {False})
            dtAtmospheric {int} -- atmospheric dt (default: {1})
            dtNormal {int} -- non atmospheric dt (default: {1})

        """
        rlist = []
        clocklist = []
        for deltat in tqdm(range((int(timeJump / abs(dtAtmospheric))) + 1)):
            # sphere of influence check
            if np.sqrt(self.r.dot(self.r)) + 1000 > self.parentRSOI:
                print("Body has left sphere of influence")
                print("Radius: ", self.r)
                print("Velocity: ", self.v)
                break
            # Check if body has crashed
            if self.parentRadius > np.sqrt(self.r.dot(self.r)):
                print("Body crashed into surface")
                print("Ending propagation")
                # print("Simulation ran for: " + str(time.time() - self.start))
                print("Radius:", np.sqrt(self.r.dot(self.r)), ":", self.r)
                print("Altitude:", self.altitude.dot(self.altitude))
                break
            
            # Use refreshByTimestep depending on wheter the body is in atmosphere
            # plus 200 is a safety margin as it's takeing the previous altitude
            if np.sqrt(self.altitude.dot(self.altitude)) < self.atmosphericLimitAltitude + 200: 
                self.refreshByTimestep(dtAtmospheric , atmospheric)
                rlist.append(self.r)
                clocklist.append(self.clock)
            else:
                self.refreshByTimestep(dtNormal, atmospheric)
                rlist.append(self.r)
                clocklist.append(self.clock)

        # Save File in CSV form if path is given
        if saveFile is not None:
            # Save position data
            rlist = np.array(rlist).T
            data = pd.DataFrame({"clock":clocklist, "x":rlist[0], "y":rlist[1], "z":rlist[2]})
            data.to_csv(saveFile, sep=",")

            # Save Manoeuvre data
            manoeuvreSavefile = saveFile[:-4] + "_man" + ".csv"
            manoeuvreData = pd.DataFrame(self.manoeuvers)
            manoeuvreData.to_csv(manoeuvreSavefile)


    # Vallado ed. 4 Algorithmn 1 page 63
    def c2_c3(self, psi):
        """Obtain c2 and c3 values used for the universal solution method of kepler's problem
        
        Arguments:
            psi {float} -- psi value for universal solution method [rad]
        
        Returns:
            float, float -- c2, c3
        """
        if psi > 10e-6:
            c2 = (1 - cos(sqrt(psi))) / psi
            c3 = (sqrt(psi) - sin(sqrt(psi))) / sqrt(psi * psi * psi)
        elif psi < -10e-6:
            c2 = (1 - cosh(sqrt(-psi)))/psi
            c3 = (sinh(sqrt(-psi) - sqrt(-psi))) / sqrt((-psi * psi * psi))
        else:
            c2 = 0.5
            c3 = 1./6.
        return c2, c3


    # Vallado ed. 4 Algorithmn 8 page 93: change of positions with time
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
        sqrtmu = sqrt(self.mu)
        rdotv = r.dot(v)

        eta = ((vnorm * vnorm) / 2) - (self.mu/rnorm)
        alpha = ((-(vnorm * vnorm)) / self.mu) + (2 / rnorm)

        if alpha > 0:
            assert (alpha != 1), "Alpha in Kepler Time is 1"
            chi = sqrtmu * dt * alpha
        elif alpha < 0:  # e>1
            alphainv = 1/alpha  # this is actually semi major axis
            chi = (np.sign(dt) * sqrt(-alphainv) *
                      np.log((-2 * self.mu * alpha * dt) /
                             (rdotv + np.sign(dt) *
                              sqrt(-self.mu / alpha) * (1 - rnorm * alpha))))
        else:
            print("parabola in keplerTime")
            chi = sqrtmu * dt / rnorm


        while True:
            chiLast = chi
            psi = chi * chi * alpha
            c2, c3 = self.c2_c3(psi)

            rr = (chi * chi * c2) + ((rdotv/sqrtmu) * chi * (1-psi*c3)) + (rnorm * (1-psi*c2))

            chi += ((sqrtmu * dt) - (chi*chi*chi*c3) - ((np.dot(r, v)/self.mu)*chi*chi*c2) - (rnorm*chi*(1-psi*c3)))/rr

            if abs(chiLast-chi) < 10e-6:
                break


        f = 1 - ((chi*chi / rnorm) * c2)
        fdot = (sqrtmu / (rr * rnorm)) * chi * ((psi*c3) - 1)
        g = dt - ((chi*chi*chi / sqrtmu) * c3)
        gdot = 1 - ((chi*chi / rr) * c2)

        rnew = (f * r) + (g * v)
        vnew = (fdot * r) + (gdot * v)
        succes = (f * gdot) - (fdot * g)
        assert abs(succes-1) < 10e-6, "Succes is not 1 but: " + str(succes)
        return rnew, vnew


    # Vallado ed. 4 Algorithm 9 page 113
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

        vnorm = np.sqrt(v.dot(v))
        rnorm = np.sqrt(r.dot(r))

        h = np.cross(r, v)
        hnorm = np.sqrt(h.dot(h))

        n = np.cross(K, h)
        nnorm = np.sqrt(n.dot(n))

        e = (( ((vnorm * vnorm)-self.mu/rnorm) * r) - (np.dot(r, v) * v)) / (self.mu)  # is a vector
        enorm = np.sqrt(e.dot(e))

        eta = (((vnorm * vnorm)/2) - (self.mu/rnorm))

        if enorm != 1.0:
            a = -self.mu / (2 * eta)
            p = a * (1 - enorm*enorm)
        else:
            p = hnorm*hnorm / self.mu
            a = inf

        i = acos(h[2]/hnorm)

        # acos should work here instead of np.arccos since Omega should not be an array
        Omega = acos(n[0]/nnorm)
        if n[1] < 0:
            Omega = 2*pi-Omega

        omega = acos(np.dot(n, e)/(nnorm*enorm))
        if e[2] < 0:
            omega = 2*pi-omega

        # nu=greek "v" used for poisson ratio
        trueAnomaly = acos(np.dot(e, r)/(enorm*rnorm))
        if np.dot(r, v) < 0:
            trueAnomaly = 2*pi-trueAnomaly
        
        omega_true, lambda_true, u = None, None, None
        if enorm < 1 and i == 0:
            omega_true = acos(e[1]/enorm)
            if e[1] < 0:
                omega_true = 2*pi-omega_true
        
        elif enorm == 0 and i != 0:
            u = acos(np.dot(n, r)/(nnorm*rnorm))
            if r[2] < 0:
                u = 2*pi-u

        elif enorm == 0 and i == 0:
            lambda_true = acos(r[0]/rnorm)
            if r[1] < 0:
                lambda_true = 2*pi-lambda_true

        return p, a, enorm, i, Omega, omega, trueAnomaly, omega_true, u, lambda_true


    # Vallado ed. 4 Algorithmn 10 page 118
    def COEtoRV(self, a, e, i, Omega, omega, trueAnomaly, omega_true=None, u=None, lambda_true=None, p=None):
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
            p = a * (1 - e*e)


        rpqw = np.array([((p * cos(trueAnomaly)) / (1 + (e * cos(trueAnomaly)))),
                         ((p * sin(trueAnomaly)) / (1 + (e * cos(trueAnomaly)))),
                         0])

        vpqw = np.array([-1 * sqrt(self.mu/p) * sin(trueAnomaly),
                         sqrt(self.mu/p) * (e + cos(trueAnomaly)),
                         0])

        PQWtoIJKtransform = np.array([[cos(Omega) * cos(omega) - sin(Omega) * sin(omega) * cos(i), -cos(Omega) * sin(omega) - sin(Omega) * cos(omega) * cos(i), sin(Omega) * sin(i)],
                                      [sin(Omega) * cos(omega) + cos(Omega) * sin(omega) * cos(i), -sin(Omega) * sin(omega) + cos(Omega) * cos(omega) * cos(i), -cos(Omega) * sin(i)],
                                      [sin(omega) * sin(i), cos(omega) * sin(i), cos(i)]])

        rijk = np.matmul(PQWtoIJKtransform, rpqw)
        vijk = np.matmul(PQWtoIJKtransform, vpqw)
        return rijk, vijk


    def addManoeuverByVector(self, clock, dv, iterSinceCreation=11):
        '''Adds a direct manoeuvre to the manoeuvre list
        Arguments:
            clock {float} -- clock of manoeuvers
            dv {ndarray} -- delta v for maneuvers in cartesian coordinates
            iterSinceCreation {int} -- (DO NOT USE) Background parameter to ensure that manoeuvres run properly


        '''
        # gets latest ID
        if len(self.manoeuvers) == 0:
            latestID = 0
        else:
            latestID = len(self.manoeuvers)

        self.manoeuvers.append({"ID": latestID,
                                "clock": clock,
                                "dv": dv,
                                "expired": False,
                                "direction": (dv/np.sqrt(dv.dot(dv))),
                                "manType": "g",
                                "iterSinceCreation": iterSinceCreation})


    def addManouvreByDirection(self, clock, dvMagnitude, manoeuvreType, iterSinceCreation=11, temp=True):
        """Adds a magnitudal manoeuvre to the manoeuvrelist.

        Arguments:
            clock {float} -- clock to perform manoeuvre
            dvMagnitude {float} -- normal value of deltav
            manoeuvreType {string} -- manoeuvre type (tangential, radial and normal) or (t,r,n)
            iterSinceCreation {int} -- (DO NOT USE) Background parameter to ensure that manoeuvres run properly
            temp {Boolean} -- (DO NOT USE) Background parameter to ensure direction of manoeuvre is correct at burn
        """
        manoeuvreTypes = ["tangential", "t", "radial", "r", "normal", "n"]
        assert (manoeuvreType in manoeuvreTypes), "Incorrect manouvreType, use tangential, t, radial, r, normal or n"

        # gets latest ID
        if len(self.manoeuvers) == 0:
            latestID = 0
        else:
            latestID = len(self.manoeuvers)

        if temp:
            if (manoeuvreType == "tangential" or manoeuvreType == "t"):
                manType = "t"
                dv = dvMagnitude
                direction = "-"
            if (manoeuvreType == "radial" or manoeuvreType == "r"):
                manType = "r"
                dv = dvMagnitude
                direction = "-"
            if (manoeuvreType == "normal" or manoeuvreType == "n"):
                manType = "n"
                dv = dvMagnitude
                direction = "-"
        
        else:
            # Create dv vector based on manoeuvretype and dvMagnitude
            if (manoeuvreType == "tangential" or manoeuvreType == "t"):
                manType = "t"
                direction = self.v / np.sqrt(self.v.dot(self.v))
                dv = dvMagnitude * direction

            if (manoeuvreType == "radial" or manoeuvreType == "r"):
                manType = "r"
                direction = -self.r / np.sqrt(self.r.dot(self.r))
                dv = dvMagnitude * direction

            if (manoeuvreType == "normal" or manoeuvreType == "n"):
                manType = "n"
                directionTangential = self.v / np.sqrt(self.v.dot(self.v))
                directionRadial = -self.r / np.sqrt(self.r.dot(self.r))
                directionNormal = np.cross(directionTangential, directionRadial)
                direction = directionNormal / \
                    np.sqrt(directionNormal.dot(directionNormal))
                dv = dvMagnitude * direction

        self.manoeuvers.append({"ID": latestID,
                                "clock": clock,
                                "dv": dv,
                                "expired": False,
                                "direction": direction,
                                "manType": manType,
                                "iterSinceCreation": iterSinceCreation})
