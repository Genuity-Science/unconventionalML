# script containing many custom controlled gate functions
# the circuit is the first parameter, c=controlled, n=number n of the controls,

# native libraries
import numpy as np
from numpy import sin, cos, tan, exp, pi, round, sign, abs, sqrt, real
from numpy import conjugate as cnj
from random import uniform, shuffle

# IBM QISkit libraries
from qiskit import QuantumProgram
from qiskit import QuantumCircuit as QC

# our libraries
import quantumtools.unitary as un
import quantumtools.qtools as qt

# import qtools.Qconfig as Qconfig
# from .unitary import * #u4powerparams, u4conjparams
# from .accessories import *
# from .qtools import *

decimalnum = 16  # small number introduces error ...
#    ideally these gates should be called the same way built in gates are, i.e. qcirc.cu1 etc.
cnt=0
parsqrtxc=un.u4conjparams(*un.parsqrtx)

# def c2x(qcirc,controlqubits,targetqubit, cvbitstring=None):
#     """ Toffoli gate (CCNOT) https://arxiv.org/pdf/0803.2316.pdf
#     Uses recursion, cvbitstring is optional string of control values for the qubits
#     """
#     assert len(controlqubits)==2, "Toffoli gate must have two control qubits" #number of control qubits
#
#     if cvbitstring is not None:
#         # invert the control qubits whose control value is zero
#         multi1qgates(qcirc, controlqubits, QC.x, cvbitstring, 0)
#
#     # "almost Toffoli: sec 6.2 in Barenco et al., elementary gates for quantum computation
#     # phase of the ctrls: 10 is off, won't see it with naive measurement
#     # qcirc.ry(pi/4,targetqubit)
#     # qcirc.cx(controlqubits[0],targetqubit)
#     # qcirc.ry(pi/4,targetqubit)
#     # qcirc.cx(controlqubits[1],targetqubit)
#     # qcirc.ry(-pi/4,targetqubit)
#     # qcirc.cx(controlqubits[0],targetqubit)
#     # qcirc.ry(-pi/4,targetqubit)
#
#     if cvbitstring is not None:
#         # invert the control qubits whose control value is zero
#         multi1qgates(qcirc, controlqubits, QC.x, cvbitstring, 0)
#     return

def cnx(qcirc,controlqubits,targetqubit, ancillaqubits=None,cvbitstring=None,ancilla_start_0=True):
    """ controlled X (i.e. NOT) with n control qubits
    based on Lemma 7.2 Barenco et al.
    ancilla_start_0: can we assume ancilla bits start at |0> ket? Saves work if true
    """
    global cnt

    if cvbitstring is not None:
        # invert the control qubits whose control value is zero
        qt.multi1qgates(qcirc, controlqubits, QC.x, cvbitstring, 0)

    nc=len(controlqubits) #number of control qubits
    # print("nc ",nc)
    if nc==1: # CNOT
        qcirc.cx(controlqubits[0],targetqubit)
    elif nc==2: # Toffoli
        qcirc.ccx(*controlqubits,targetqubit)
    else: # nc > 2
        if ancillaqubits is not None: # use ancilla qubits
            num_ancilla = nc - 2
            assert len(ancillaqubits) >= nc -2, "{} ancilla qubits is not enough for {} controls".format(num_ancilla,nc)
            # above ensures at least one ancilla qubit
            # chop off extra ancilla qubits
            ancillaqubits=ancillaqubits[0:num_ancilla]

            # 8 pieces gate groups, 4 are single Toffoli and 4 are loops of num_ancilla-1 Toffolis.
            # 3 pieces (first 2 and last 1) executed only if not ancilla_start_0
            if not ancilla_start_0:
                qcirc.ccx(ancillaqubits[0],controlqubits[0],targetqubit)
                for i in range(1,num_ancilla):
                    qcirc.ccx(ancillaqubits[i],controlqubits[i],ancillaqubits[i-1])

            # main block, executed regardless of initial value of ancillas
            # print("controlqubits[-2:] ", controlqubits[-2])
            # print("ancillaqubits[-1] ", type(ancillaqubits[-1]))
            qcirc.ccx(*controlqubits[-2:],ancillaqubits[-1])
            for i in reversed(range(1,num_ancilla)):
                qcirc.ccx(ancillaqubits[i],controlqubits[i],ancillaqubits[i-1])
            qcirc.ccx(ancillaqubits[0],controlqubits[0],targetqubit)
            for i in range(1,num_ancilla):
                qcirc.ccx(ancillaqubits[i],controlqubits[i],ancillaqubits[i-1])
            qcirc.ccx(*controlqubits[-2:],ancillaqubits[-1])

            if not ancilla_start_0: # if we cannot assume ancilla bits start at 0
                for i in reversed(range(1,num_ancilla)):
                    qcirc.ccx(ancillaqubits[i],controlqubits[i],ancillaqubits[i-1])
            # cnt+=1
            # print('linear ',cnt)

        else: # no ancilla qubits
            # TODO replace this with the O(2^n) model, for the case if no ancilla.
            # in which case the algorithm is based on Gray codes

            # Below is inefficient model, O(3^n) ....
            # the last control qubit is used as control, and removed from the next cnu4
            # five gates: by number of controls, 1, nc-1, 1, nc-1, nc-1

            # separate the last control qubit from the rest
            firstcontrolqubit = controlqubits[0]
            othercontrolqubits = controlqubits[1:]

            # The five gates in sequence, one of them recursively calls cnu4
            cu4(qcirc,un.parsqrtx,firstcontrolqubit,targetqubit) #ctrl-sqrt(U)
            cnx(qcirc, othercontrolqubits, firstcontrolqubit,ancillaqubits)
            cu4(qcirc,parsqrtxc,firstcontrolqubit,targetqubit) #ctrl-sqrt(conj(U))
            cnx(qcirc, othercontrolqubits, firstcontrolqubit,ancillaqubits)
            cnu4(qcirc, un.parsqrtx, othercontrolqubits, targetqubit,ancillaqubits) #ctrl-n-sqrt(U)
            #cnt+=1
            #print('exponential ',cnt)

    if cvbitstring is not None:
        # invert the control qubits whose control value is zero
        qt.multi1qgates(qcirc, controlqubits, QC.x, cvbitstring, 0)
    return

def crx(qcirc,theta,controlqubit,targetqubit,cv=1):
    """ controlled Rx (i.e. X rotation)
    cv=control value
    """

    if cv == 0:
        qcirc.x(controlqubit)
    qcirc.u1(beta,controlqubit)
    qcirc.cu3(theta,phi,lam,controlqubit,targetqubit)
    if cv == 0:
        qcirc.x(controlqubit)
    return

def cu4(qcirc,par,controlqubit,targetqubit,cv=1):
    """ controlled U4 (i.e. most general single qubit unitary, with global phase)
    ideally this should be called the same way built in gates are, i.e. qcirc.cu1 etc.
    cv=control value
    """
    (theta,phi,lam,beta)=par
    if cv == 0:
        qcirc.x(controlqubit)
    qcirc.u1(beta,controlqubit)
    qcirc.cu3(theta,phi,lam,controlqubit,targetqubit)
    if cv == 0:
        qcirc.x(controlqubit)
    return

def cnu4(qcirc,par,controlqubits,targetqubit, ancillaqubits=None,cvbitstring=None):
    """ controlled U4 (i.e. most general single qubit unitary, with global phase)
    with n control qubits
    Uses recursion, based on lemma 7.5 in Barenco et al., elementary gates for quantum computation
    ideally this should be called the same way built in gates are, i.e. qcirc.cu1 etc.
    cvbitstring is optional string of control values for the qubits
    """

    if cvbitstring is not None:
        # invert the control qubits whose control value is zero
        qt.multi1qgates(qcirc, controlqubits, QC.x, cvbitstring, 0)

    nc=len(controlqubits) #number of control qubits
    if nc==1: #base case
        cu4(qcirc,par,controlqubits[0],targetqubit)
    else: #nc > 1
        # the last control qubit is used as control, and removed from the next cnu4
        # five gates: by number of controls, 1, nc-1, 1, nc-1, nc-1

        # separate the last control qubit from the rest
        firstcontrolqubit = controlqubits[0]
        othercontrolqubits = controlqubits[1:]

        # find parameters for the square root gate and its conjugate
        parsqrt = un.u4powerparams(*par,0.5)
        parsqrtc = un.u4conjparams(*parsqrt)

        # The five gates in sequence, one of them recursively calls cnu4
        cu4(qcirc,parsqrt,firstcontrolqubit,targetqubit) #ctrl-sqrt(U)
        cnx(qcirc, othercontrolqubits, firstcontrolqubit,ancillaqubits)
        cu4(qcirc,parsqrtc,firstcontrolqubit,targetqubit) #ctrl-sqrt(conj(U))
        cnx(qcirc, othercontrolqubits, firstcontrolqubit,ancillaqubits)
        cnu4(qcirc, parsqrt, othercontrolqubits, targetqubit,ancillaqubits) #ctrl-n-sqrt(U)


    if cvbitstring is not None:
        # invert the control qubits whose control value is zero
        qt.multi1qgates(qcirc, controlqubits, QC.x, cvbitstring, 0)
    return
