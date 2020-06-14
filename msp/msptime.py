import numpy as np
import datetime
from math import radians, degrees, sin, asin, tan, atan, cos, acos, pi, floor

def JulianDate(yr,mo,d,h,m,s):
    return 367*yr - int((7 * (yr + int((mo+9)/12) ))/4) + int((275*mo)/9) + d + 1721013.5 + (((s/60 + m)/60 + h) / 24)

def JulianDateFromdatetime(dtobj):
    date = (dtobj.year, dtobj.month, dtobj.day, dtobj.hour,
            dtobj.minute, dtobj.second + 10**-6 * dtobj.microsecond)
    return JulianDate(*date)
    
def LSTTime(JDut1, lmda):
    Tut1 = (JDut1 - 2451545.0) / 36525
    thetaGMST = ((67310.54841 + ((876600 * 3600) + 8640184.812866) * Tut1 +
                  0.093104 * Tut1*Tut1 - 6.2e-6 * Tut1*Tut1*Tut1) % 86400) / 240
    thetaLST = thetaGMST + lmda
    return thetaGMST, thetaLST

def DEGtoHMS(angle):
    temp = radians(angle) * 12/(pi)
    h = floor(temp)
    m = floor(60 * (temp-h))
    s = (temp-h-(m/60)) * 3600
    return h, m, s

def HMStoSeconds(hms):
    return hms[0] * 3600 + hms[1] * 60 + hms[2]

def SecondsToHMS(s):
    temp = s / 3600
    h = floor(temp)
    m = floor(60 * (temp-h))
    s = (temp-h-(m/60)) * 3600
    return h, m, s

def SunriseSunset(date,phigc,lmda, sigma):
    JDut1 = JulianDateFromdatetime(date)
    JDsunrise = JDut1 + (6/24) - (lmda/360)

    def defaultFunc(JDspec):
        Tut1 = (JDspec - 2451545.0) / 36525
        Ttdb = Tut1

        lmdMeanSun = (280.4606184 + 36000.77005361 * Tut1) % 360
        MeanAnomSun = (357.5291092 + 35999.05034 * Ttdb) % 360
        lmdEcliptic = lmdMeanSun + 1.914666471 * sin(radians(MeanAnomSun)) + 0.019994643 * sin(2 * radians(MeanAnomSun))
        eta = 23.439291 - 0.0130042 * Ttdb

        alphaSun = degrees(atan(cos(radians(eta)) * tan(radians(lmdEcliptic))))
        deltaSun = degrees(asin(sin(radians(eta)) * sin(radians(lmdEcliptic))))

        LHAsunset = degrees(acos((cos(radians(sigma)) - sin(radians(deltaSun)) * sin(radians(phigc)))/(cos(radians(deltaSun)) * cos(radians(phigc)))))
        LHAsunrise = 360 - LHAsunset 

        thetaGMST = (100.4606184 + 36000.7700536 * Tut1 + 0.00038793 * Tut1 * Tut1 - 2.6e-8 * Tut1 * Tut1 * Tut1) % 360

        UTsunrise = (LHAsunrise + alphaSun - thetaGMST) % 360
        UTsunset = (LHAsunset + alphaSun - thetaGMST) % 360
        return UTsunrise,UTsunset
    
    return defaultFunc(JDsunrise)

def getSunlight(date, latitude, longitude):
    sunrise, sunset = SunriseSunset(date, latitude, longitude, 90+5/6)
    return DEGtoHMS(sunset-sunrise)

def TDB(yr, mo, d, h, m, s, UTC, dUT1=0, dAT=0):
    # dAT and dUT1 are to be found online or in an almanac for initial estimates they are sufficient left at 0
    UTC = (3600 * (h + UTC)) + (60 * m) + s
    UT1 = UTC + dUT1
    TAI = UTC + dAT 
    GPS = UTC + dAT - 19
    TT = TAI + 32.184
    # I only included the first term. Accurate to about 6 digits
    TDB = TT + 0.001657*sin(628.3076 * TT + 6.24) 
    hmsTDB = SecondsToHMS(TDB)
    JDTDB = JulianDate(yr, mo, d, *hmsTDB)
    return JDTDB
# Vallado algorithm 16


# # test = datetime.datetime(1992, 8, 20, 12, 14, 0)
# # test = datetime.datetime(1992, 8, 20, 12, 14, 0)
# # test = datetime.datetime(1996, 3, 23, 0, 0 ,0 ,0)
# # JDtest = (JulianDateFromdatetime(test))
# # print(LSTTime(JDtest, 0))

# date = datetime.datetime(1996, 3, 23, 0, 0 ,0 ,0)
date = datetime.datetime(2020, 6, 6, 0, 0 ,0 ,0)
longitude = 4.3555600
latitude = 52.0066700

result = (SunriseSunset(date, latitude, longitude, 90+5/6))
print(f"Sunrise: {DEGtoHMS(result[0] + 25.1042)}, Sunset: {DEGtoHMS(result[1] + 25.1042)}")
print(f"Hours: {getSunlight(date, latitude, longitude)} total seconds: {HMStoSeconds(getSunlight(date, latitude, longitude))}")
# print(DEGtoHMS(result))