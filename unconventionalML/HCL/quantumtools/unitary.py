# script about manipulating custom single qubit unitary gates, by converting
# between different parameterizations in QISKIT, using powers and conjugates.
# Author: Omar Gamel, 2018

# native libraries
import numpy as np
from numpy import sin, cos, tan, exp, pi, round, sign, abs, sqrt, real
from numpy import conjugate as cnj
from cmath import polar, phase

# our libraries
import quantumtools.accessories as ac
#from .accessories import *

#### Read the guide to understand the various unitary parameterizations

# important gate parameterset par=(theta,phi,lam,beta)
parx=(  pi, 0,  pi, pi/2) # Pauli X
pary=(  pi, 0,   0, pi/2) # Pauli Y
parz=(   0, 0,  pi, pi/2) # Pauli Z
parh=(pi/2, 0,  pi, pi/2) # Hadamard
pars=(   0, 0, pi/2,pi/4) # Phase gate
part=(   0, 0, pi/4,pi/2) # Phase gate
parsqrtx=(  pi/2, -pi/2,  pi/2, pi/4) # sqrt(X) gate

####################################################################
#### Functions creating unitary matrices from parameterizations ####

def urotm(alpha,a):
    """returns unitary matrix with in the input parameters angle alpha and unit vector a"""
    c,s=cos(alpha/2),sin(alpha/2)

    return np.array([[c-1j*s*a[2],      -s*(a[1]+1j*a[0])],
                     [s*(a[1]-1j*a[0]),      c+1j*s*a[2]]])

def urotphasem(alpha,a, beta):
    """returns unitary matrix with in the input parameters angle alpha, unit vector a,
    and global phase beta """
    return r(exp(1j*beta)*urotm(alpha,a))

def ubasem(theta,phi,lam):
    """returns U matrix with in the input parameters """
    c,s=cos(theta/2),sin(theta/2)
    psum=exp((phi+lam)/2 * 1j)
    pdiff = exp((phi - lam) / 2 * 1j)

    return np.array([[cnj(psum)*c, -cnj(pdiff)*s],
                     [pdiff*s,            psum*c]])

def u1m(lam):
    """returns u1 matrix with the input parameters """
    return u3m(0,0,lam)

def u2m(phi,lam):
    """returns u2 matrix with the input parameters """
    return u3m(pi/2,phi,lam)

def u3m(theta,phi,lam):
    """returns u3 matrix with the input parameters """
    return r(exp((phi+lam)/2 * 1j)*ubasem(theta,phi,lam))

def u4m(theta,phi,lam,beta):
    """returns newly defined u4 matrix with the input parameters """
    return r(exp(1j*beta) * ubasem(theta,phi,lam))


####################################################################
#### Functions finding parameterizations from unitary matrices  ####

def getrotparams(U):
    """find paramaters such that urotphasem(alpha,a, beta) = U """
    assert U.shape == (2,2), "use single qubit gate"
    assert ac.equaltol(np.dot(U,np.matrix(U).getH()),np.eye(2)), "input must be a unitary matrix"

    # alpha in [0, pi], both c=cosine(alpha/2) and s=sine(alpha/2) postive
    (c,beta)=polar((U[0,0]+U[1,1])/2)
    alpha = 2 * np.arccos(c)
    s=sqrt(1-c**2)

    # remove global phase
    U = exp(-beta*1j) * U
    a = np.array([0.0,0.0,0.0])
    # RHS below is real if U is unitary, but use real() to suppress warning
    a[0] = real((U[1,0]+U[0,1])*1j/(2 * s))
    a[1] = real((U[1,0]-U[0,1])/(2 * s))
    a[2] = real((U[0,0]-U[1,1])*1j/(2 * s))

    return alpha,a, beta

def getu4params(U):
    """find paramaters such that urotphasem(alpha,a, beta) = U """
    assert U.shape == (2,2), "use single qubit gate"
    assert ac.equaltol(np.dot(U,np.matrix(U).getH()),np.eye(2)), "input must be a unitary matrix"

    # theta in [0, pi], both c=cosine(theta/2) and s=sine(theta/2) postive
    # beta in [0, pi] and any ovrall -ve sign absorbed into halfsum = (phi+lam)/2
    # or halfdiff = (phi-lam)/2 both shifting by +/- pi
    (c2,beta2)=polar(U[0,0]*U[1,1])
    theta = np.arccos(2*c2-1) #c2=cos^2(theta/2), then use double angle formula
    beta = ac.setanglep(beta2)/2

    # remove global phase
    U = exp(-beta*1j) * U

    halfsum = phase(U[1,1])
    halfdiff = phase(U[1,0])

    phi = halfsum + halfdiff
    lam = halfsum - halfdiff

    return theta,phi,lam,beta


#############################################################################
#### Main functions returning the parameterizations of operations on u4  ####

def u4conjparams(theta,phi,lam,beta):
    """ return the parameters of the conjugate of u4, such that
    u4m(theta,phi,lam,beta)^dagger = u4m(ctheta,cphi,clam,cbeta)
    """
    ctheta = theta
    cbeta = -beta

    chalfsum = -(phi + lam)/2
    chalfdiff = (phi - lam)/2
    #chalfdiff = movebypi(chalfdiff) # if non-negative, subtract pi, otherwise add it

    if cbeta < 0: # ensure unique parameterization that cbeta in [0, pi]
        cbeta = ac.movebypi(cbeta)
        chalfsum = ac.movebypi(chalfsum)
    else:
        chalfdiff = ac.movebypi(chalfdiff)

    cphi = chalfsum + chalfdiff
    clam = chalfsum - chalfdiff

    return (ctheta,cphi,clam,cbeta)

def u4powerparams(theta,phi,lam,beta,p):
    """ return the parameters of the powered u4
    u4m(theta,phi,lam,beta)^p = u4m(ptheta,pphi,plam,pbeta) ... or more precisely
    u4m(theta,phi,lam,beta)^n = u4m(ptheta,pphi,plam,pbeta)^m ... where p=n/m, rational fraction
    """

    alpha = 2 * np.arccos(cos((phi + lam) / 2) * cos(theta / 2)) # in [0, 2pi]
    # handle case alpha = 0
    sinratio = sin(p * alpha / 2)/ sin(alpha / 2) if not ac.equaltol(alpha,0) else p
    phalfdiff = (phi - lam) / 2  # unchanged

    # get ptheta. 2*arcsin will result in [-pi, pi], then fix it to [0,pi] unique parameter range
    ptheta=2*np.arcsin(sinratio * sin(theta/2))
    if ptheta < 0:
        ptheta *= -1
        phalfdiff = ac.movebypi(phalfdiff)

    # get halfsum
    if ac.equaltol(cos(ptheta/2),0):
        phalfsum = 0 # in this case phalfsum doesn't matter, so set it to zero to avoid errors
    else:
        cosphalfsum = cos(p * alpha /2)/cos(ptheta/2)
        if abs(cosphalfsum) > 1:
            # python floating point inaccuracy can lead to this condition slightly being true,
            # obstructing the ensuing arccos. Fix it with a simple division.
            cosphalfsum/=abs(cosphalfsum)

        # use arccos to get phalfsum, result in [0, pi], then fix it to [-pi, pi] using additional
        # information of sign of sin(phalfsum) =  cos(p*alpha/2)/cos(ptheta/2) to get the quadrant
        phalfsum = np.arccos(cosphalfsum)
        phalfsum *= ac.sgnp(sinratio*sin((phi+lam)/2)*cos(theta/2)/cos(ptheta/2))

    # global phase multipied by power, make result in [-pi, pi],
    # then fix it to [0,pi] unique parameter range
    pbeta = ac.setanglepm(p * beta)
    if pbeta < 0:
        pbeta = ac.movebypi(pbeta)
        phalfsum = ac.movebypi(phalfsum)
        phalfdiff = ac.movebypi(phalfdiff)

    pphi = phalfsum + phalfdiff
    plam = phalfsum - phalfdiff

    # # debugging: print output
    # print('alpha/2:',alpha/2)
    # print('sin ratio:',sinratio)
    # print('ptheta:', ptheta)
    # print('phalfsum:', phalfsum)
    # print('phalfdiff', phalfdiff)

    return (ptheta,pphi,plam,pbeta)

# alhamdolela, alhamdolela, alhamdolela
# it was an annoying problem, took the better part of two weeks to get this going, was totally
# losing it one point ... bifadlillah it worked out with tawakkul 3alaAllah, persistence and confidence
# and the feeling of accomplishment is that much better after the struggle!
