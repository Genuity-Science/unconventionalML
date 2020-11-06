# qtools are specifically quantum functions that are reused often in
# many algorithms

# native libraries
import numpy as np
from itertools import product
import math
from numpy.linalg import norm
from numpy import power, ceil, arccos, sign, pi

# IBM QISkit libraries
from qiskit import QuantumProgram
from qiskit import QuantumCircuit as QC

# our libraries
import quantumtools.accessories as ac
# import qtools.unitary as unt

# from qtools.unitary import u4powerparams, u4conjparams
# from qtools.gates import cnu4
# from qtools.accessories import *

def initializeQenv():
    pass
    return

def numGroveriterations(numstates, numtargets=1):
    """ The theoretical maximum number of Grover iterations needed"""
    return int(np.sqrt(numstates/numtargets) * pi/4)

def numGroverOptimizeiterations(numitems):
    """ The theoretical maximum number of Grover Optimize iterations needed"""
    return int(22.5*np.sqrt(numitems)+1.4*np.log(numitems)**2) if numitems>0 else 0

def getqubitlist(qreg):
    """list of qubits from quantum register, this facilitates manipulating them"""
    qubits = [qreg[i] for i in range(0,qreg.size)]
    return qubits

def multi1qgates(qcirc, qubits, gate, bitstring=None,flag=0):
    """ Apply the same single qubit gate to all qubits[i] where
    bitstring[i] == flag, or to all qubits if bitstring unspecified
    """
    # usually called with bitstring flipped, so it matches the qubits
    nbits = len(qubits)
    looprange = range(0,nbits)

    # indicator shows whether or not to apply gate to each qubit
    applygate = [True for i in looprange]

    if bitstring is not None:

        # check size
        assert len(bitstring) == nbits,"number of qubits must match bitstring length"
        # count bitstring indices in reverse, since indexing of strings starts
        # from the left, while indexing of qubits starts from the right
        applygate = [int(bitstring[i]) == flag for i in range(-1,-1-nbits,-1)]
        #print("applygate ",applygate)

    for i in looprange:
        if applygate[i]:
            gate(qcirc,qubits[i])
    return

def allqubitsX(qcirc, qubits):
    """ Apply X to all qubits
    """
    multi1qgates(qcirc, qubits, QC.x)
    return

def allqubitsH(qcirc, qubits):
    """ Apply H to all qubits
    """
    multi1qgates(qcirc, qubits, QC.h)
    return

###############################################################################
###############################################################################
# the next three functions create the initial superposition from a vector of
# amplitudes. quantizevec(vec) is the main function that uses the other two

def combinehalf(vec):
    """halves the length of the vector by adding adjacent pairs of elements """
    N=len(vec)
    assert N%2 == 0, "vector must have even length"

    return np.array([vec[2*i]+vec[2*i+1] for i in range(0,int(N/2))],dtype=float)

def thetahierarchy(vec):
    """returns the hierarchy of angle rotations needed to get to the vec amplitudes
    from a starting amplitude of (1,0,0,....,0), i.e the |000...0> state """
    #need divide each by parent weight so each pairs add to one, then translate that to an angle
    N=len(vec)
    nbits = np.log2(N)
    assert nbits.is_integer(), "vector length must be power of 2"
    nbits=int(nbits)

    theta_hier=[]
    lasthalfed=vec
    for i in range(0,nbits):
        newhalfed=combinehalf(lasthalfed)
        # double angle formula. arccos means each theta is between 0 and pi inclusive
        # lasthalfed[2*j]/newhalfed[j] = cos^2(theta/2), lasthalfed[2*j+1]/newhalfed[j] = sin^2(theta/2),
        # as per definition of current theta, where newhalfed[j] = lasthalfed[2*j] + lasthalfed[2*j+1]
        thetas=np.array([(arccos((lasthalfed[2*j]-lasthalfed[2*j+1])/newhalfed[j]) if newhalfed[j] !=0 else 0) for j in range(0, len(newhalfed))], dtype=float)
        theta_hier.append(thetas)
        lasthalfed=newhalfed

    # theta_hier.reverse() #don't reverse to stay in same order as qubits
    return theta_hier

# a dictionary that helps put the signs back into the amplitudes
signstophilambda={(1,1):(0,0),(1,-1):(pi,-pi),(-1,1):(pi,pi),(-1,-1):(2*pi,0)}
# this assumes b is real .. in general it can be complex
# to do make it so b is complex, and its phase is applied in the final step (in pairs)

def quantizevec(qcirc, qubits, b,ancillaqubits=None):
    """takes unit vector b and makes a quantum superposition state \sum_i b_i |i>
    assumes register starts all zeros
    Based on https://arxiv.org/abs/quant-ph/0208112v1, Grover and Rudolph"""

    N = len(b)
    assert N>1, "vector b must have at least two entries"
    nbits = int(np.ceil(np.log2(N)))
    assert ac.equaltol(norm(b),1), "b must be a unit vector" #if not, we normalize

    # pad the vector b with zeros to nearest power of 2
    b = np.append(np.array(b), np.zeros(2 ** nbits-N))
    b2 = power(b,2)
    bs = ac.sgnp(b) #signs

    # theta hierarchary
    theta_hier=thetahierarchy(b2)

    # indexing of qubits starts from the right. Then we need to start applying gates with left-most qubit,
    # whose digit represents the highest value in the binary number
    # start with left-most qubit, which is the last qubit, labeled qubit nbits-1 (or -1)
    qcirc.u3(theta_hier[nbits-1][0],0,0,qubits[nbits-1])

    # loop through remaining qubits, left to right (largest to smallest)
    for i in reversed(range(0,nbits-1)): #i = hierarchy level
        for ind,theta in enumerate(theta_hier[i]):
            # there are 2**i thetas in theta_hier[i]
            # apply cu3 or cnu3 (effectively cnu4 with beta=0), with angle theta and binary rep of ind
            # setting the control values of control qubits, which are all qubits higher than (i.e. on left of) current qubit
            cvbitstring=('{0:0' + str(nbits-1-i) + 'b}').format(ind) #
            if i > 0:
                par=(theta,0,0,0) # parameters for simple 1 qubit rotation
            else: # i = 0, last one in the reversed loop (right most qubit), incorporate the signs in the control bits
                angles=signstophilambda[tuple(bs[2*ind:2*ind+2])] #get two signs
                par=(theta,angles[0],angles[1],0) # get sign info from bs to fix parameters

            import quantumtools.gates as gt #import here to avoid error due to circular dependency
            # use all higher qubits as control
            gt.cnu4(qcirc,par,controlqubits=qubits[i+1:],targetqubit=qubits[i], ancillaqubits=ancillaqubits, cvbitstring=cvbitstring)

    return

###############################################################################
###############################################################################
