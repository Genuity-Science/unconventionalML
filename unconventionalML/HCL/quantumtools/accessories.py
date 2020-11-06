# Accessory functions helpful in many contexts

from numpy import sin, cos, tan, exp, pi, round, sign, abs, sqrt, real
import numpy as np
from random import uniform
import time

decimalnum = 16  # small number introduces error ...
tolerance = 10 ** -7

############################
#### Accessory functions ####

def r(inp):
    """Round to preset number of decimals"""
    return round(inp, decimals=decimalnum)

def ruv():
    """Random unit vector"""
    # properties of distribution irrelevant, it is not necessarily uniformly
    # distributed on unit sphere
    a = [uniform(-1, 1) for i in range(0, 3)];
    a = a / np.linalg.norm(a)
    return a

def randompar():
    """return random set of u4 parameters"""
    theta = uniform(0,pi)
    halfsum = uniform(-pi,pi)
    halfdiff = uniform(-pi,pi)
    phi = halfsum + halfdiff
    lam = halfsum - halfdiff
    beta = uniform(0,pi)
    return (theta,phi,lam,beta)

def randomparlist(trials=20):
    return [randompar() for i in range(0,trials)]

def sgnp(num):
    """the sign function with sgn(0)=1"""
    return (2 * (num >= 0) - 1)

def setanglep(angle):
    """make input angle in [0,2pi] """
    angle = angle  % (2 * pi)
    return angle

def setanglepm(angle):
    """make input angle in [-pi,pi] """
    angle = angle % (2 * pi)
    angle -= (angle > pi) * 2 * pi
    return angle

def movebypi(angle):
    """if non-negative, subtract pi, otherwise add it"""
    # print('inverted ptheta ')
    return angle-sgnp(angle) * pi

def equaltol(a,b,tol=tolerance):
    """return true if two inputs arrays are equal entrywise within the tolerance"""
    # incorporate math.is_close (relative tolerance better than absolute)
    return (abs(a-b) < tolerance).all()

def timetounits(seconds):
    """ convert seconds to (days,hours,minutes,seconds) """
    seconds = round(seconds,2)
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    d, h = divmod(h, 24)
    return (d, h, m, s)

def timestring(seconds):
    """ output execution in human readable string format,
    e.g. 1 day 3 hours 20 minuts and 24.34 seconds
    skips unit if 0 """
    (d, h, m, s) = timetounits(seconds)
    out=""
    if d != 0:
        out += "%d d " % d
    if h != 0:
        out += "%d h " % h
    if m != 0:
        out += "%d m " % m
    if s != 0:
        out += "%.2f s." % s
    return out
