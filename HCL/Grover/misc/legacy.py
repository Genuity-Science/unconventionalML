# DISCARDED CODE That is too involved to delete
# Jan 2018

# FROM unitary.py      FUNCTION: u4powerparams
# use arctan to get phalfsum, result in [-pi/2, pi/2], then fix it to [-pi, pi] using additional
# information of sign of cos(phalfsum) =  cos(p*alpha/2)/cos(ptheta/2) to get the quadrant
# case of phi+lam =/- pi must b treated differently
sdic = {-2: -1, -1: 1, 0: 1, 1: -1, 2: -1}
tanratio = tan(p * alpha / 2) / tan(alpha / 2)
phalfsum = np.arctan(tanratio * tan((phi+lam)/2))
s=sign(cos(p * alpha/2)/cos(ptheta/2))
q=phalfsum//(pi/2) # will be -2,-1,0,1 depending on quadrant, (or 2 if phalfsum = pi)
if s != sdic[q]:
    phalfsum += sgnp(phalfsum) * pi
    print('flipped s',s,' q ',q)

# ---------------------
# FROM check.py
def checku4powerparams_old(par,p=2.0):
    """this is not a good test, since an nth root of the nth power need not give
    back original input there are multiple possible roots, and the function
    u4powerparams gives just one
    """
    parp=u4powerparams(*par,p)
    parpp=r(u4powerparams(*parp,1/p))
    return parp,parpp,equaltol(par,parpp)
    # mat1=u4m(*par)
    # mat2=u4m(*parpp)
    # return round(mat1,decimals=decimalnum),round(mat2,decimals=decimalnum),(mat1==mat2).all()

# ---------------------
# FROM Grover.py ... no longer use the bitstringmap approach. It comes from pyquils module, where it was
# used to define a whole matrix for a unitary gate for the Grover oracle.... but I just used ctrl-n-z instead
from itertools import product
alphabet = "01"
target_bitstring_phase = -1
nontarget_bitstring_phase = 1

## Formulate the Problem ##
def makebitstringmap(target_bitstring):
    """
    Make a bitstringmap dictionary where bitstrings are keys with value 1,
    and only the target_bitsring has value -1
    """
    # e.g.  = '10001000'

    assert target_bitstring.strip(alphabet) is '', "target bitstring must be made of " + alphabet

    bit = tuple(set(alphabet))
    bitstring_map = {}

    # We construct the bitmap for the oracle
    for bitstring in product(bit, repeat=len(target_bitstring)):
        bitstring = "".join(bitstring)
        if bitstring == target_bitstring:
            bitstring_map[bitstring] = target_bitstring_phase
        else:
            bitstring_map[bitstring] = nontarget_bitstring_phase

    return bitstring_map

# ---------------------
def classicfind(bitstringmap):
    """"
    Use classical loop to return all keys in dictionary bitstringmap whose value is the target
    """
    # may be return only first key, or all keys in a different format?

    return [bitstring for bitstring, phase in bitstringmap.items() if phase == target_bitstring_phase]
classicfind.__doc__ +=  str(target_bitstring_phase)

# ---------------------
# FROM Grover.groverfind()
    # modify for ibmqx5 topology, where there is only a 4 qubit ctrl cascade
    if backendcode==2:
        qubits=qubits[1:]
        # remove qubit [0], then 1->2->3->4 forms ctrl cascade (max is 4)
        # this means nbits has to be larger than usual by 1
    # turns out this doesn't work, because ibmqx5 ctrl-not must only be to
    # adjacent qubit, not any qubit lower on the cascade
# ---------------------

# FROM qtools, support functions for quantizevec
# def combinehalfhierarchy(vec):
#     # takes vector of probabilities (i.e. squared amplitudes)
#     N = len(vec)
#     nbits = int(np.log2(N))
#
#     veccombined=[vec]
#     for i in range(0,nbits):
#         veccombined.append(combinehalf(veccombined[i]))
#
#     #need divide each by parent weight so each pairs add to one, then translate that to an angle
#     for i in range(0,nbits):
#          for j in range(0,len(veccombined[i])):
#              veccombined[i][j] /=veccombined[i+1][int(j/2)]
#
#     #veccombined.reverse()
#     return veccombined
