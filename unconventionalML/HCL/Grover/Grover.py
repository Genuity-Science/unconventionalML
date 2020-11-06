# native libraries
import numpy as np
import math
import random
import operator

# IBM QISkit libraries
from qiskit import QuantumProgram
import quantumtools.Qconfig as Qconfig

# our libraries
from quantumtools.unitary import u4powerparams, u4conjparams
from quantumtools.gates import *
from quantumtools.qtools import *

# backendcode: 0: local simulator, 1: ibmqx simulator, 2: ibmqx5
backendarray=['local_qasm_simulator_cpp','ibmqx5','local_statevector_simulator','local_qasm_simulator','ibmqx_qasm_simulator','ibmqx_hpc_qasm_simulator', 'ibmqx5']

# make file for Grover oracles

def oracle_mark(qcirc, qubits,ancillaqubits=None,target_bitstrings=None):
    """ execute an artificial oracle for grover's algorithm on qcirc, i.e.
    it flips the phase for states that match one of target (marked) bitstrings, which are
    specified internally"
    The oracle is the 'given' that Grover's algorithm uses to find the target bitstrings
    To artificially construct an oracle, one of course uses the target bitstrings!
    Then finds it using Grover's algorithm. In practical cases, the oracle would
    naturally come out of the problem somehow
    """

    # the target bitstring should ideally be defined in here, to truly keep this
    # oracle a black box. But for convenience for testing purposes, we can make
    # it a parameter

    if target_bitstrings==None:
        target_bitstrings =['000']
    nstringbits = len(target_bitstrings[0])

    #nstrings=len(target_bitstrings)
    nqubits=len(qubits) #make sure we have enough qubits for the string
    assert nqubits >= nstringbits, "bitstring too long for " + str(nqubits) + " qubits"

    # apply ctrl-n-z for each string as the control
    # TODO: make this smarter, e.g. if 110 and 111 both target bitstrings, then
    # it should simply apply to left (last) two qubits 11. if 1100 and 1111, then a
    # shortcut - first a ctrl[0] -> [1], then n-ctrl on left 3 qubits 110, then a
    # ctrl[0] -> [1] again.
    # It is easier when need to flip a power of 2 number of bitstrings, but harder otherwise
    for bitstring in target_bitstrings:
        # all bitstrings must be same length
        assert len(bitstring)==nstringbits, "bitstring " + bitstring + " mismatched length"

        # controlled z gate is symmetric
        # for ctrl-n-z, choice of target qubit is arbitrary, so use first qubit as target
        cvbitstring=bitstring[0:-1]
        tgt=int(bitstring[-1]) # first qubit is last entry in bitstring

        if not tgt: # flip the target qubit if its value in bitstring (control) is 0
            qcirc.x(qubits[0])

        # apply the controlled phase (only changes what matches bitstring
        #replaced
        # cnu4(qcirc, parz, qubits[1:nstringbits], qubits[0],ancillaqubits,cvbitstring)
        # with these three
        qcirc.h(qubits[0])
        cnx(qcirc, qubits[1:nstringbits], qubits[0],ancillaqubits,cvbitstring)
        qcirc.h(qubits[0])

        if not tgt: # flip the target qubit back if its value in bitstring (control) is 0
            qcirc.x(qubits[0])
    return

    # def oracle(qcirc,qubits):
    #     """the oracle to call in grover's algorithm. target specified internally"""
    #     tgts = ["1100"]
    #     artificialoracle(qcirc,qubits,tgts)/
    #    / return
    # the function implementing Grover's algorithm

def groverdiffusion(qcirc,qubits,ancillaqubits=None):
    """The Grover Diffusion operator, 2|s><s|-I, where |s> is the equal
    superposition of all n-bit binary numbers """
    # 2|s><s|-I = H^{oxn} (2|0><0|-I) H^{oxn}
    # 2|0><0|-I = X^{oxn} ctrl-n-z X^{oxn} (up to overall) minus sign
    # overall minus can be applied through Rz, or ignored as irrelevant global phase

    allqubitsH(qcirc,qubits)
    allqubitsX(qcirc,qubits)
    # qcirc.rz(pi,qubits[0]) # i global phase
    # for ctrl-n-z, choice of target qubit is arbitrary, so use first one
    # replaced
    # cnu4(qcirc, parz, qubits[1:], qubits[0],ancillaqubits)
    # with these three
    qcirc.h(qubits[0])
    cnx(qcirc, qubits[1:], qubits[0],ancillaqubits)
    qcirc.h(qubits[0])

    # qcirc.rz(pi,qubits[0]) # i global phase
    allqubitsX(qcirc,qubits)
    allqubitsH(qcirc,qubits)
    return

# the function implementing Grover's algorithm
def groverfind(oracle,nqubits,backendcode=0,numsolutions=1,shots=1024,withancilla=True,target_bitstrings=None):
    """"
        inputs are the oracle (function), and the number of bits
        This function is the one that is called to execute Grover's algorithm
    """
    showancilla=False
    starttime=time.time()
    numiterations = numGroveriterations(2**nqubits,numsolutions)

    if nqubits <= 3:# no need for ancilla if too few qubits
        withancilla=False

    #for numiterations in range(1,numGroveriterations(2**nqubits,numsolutions)):

    # set up quantum environment. (design): is this the right place to define these?
    backend = backendarray[backendcode]
    qprog = QuantumProgram()
    qprog.set_api(Qconfig.APItoken,Qconfig.config['url'])
    qreg = qprog.create_quantum_register("qr", nqubits)
    creg = qprog.create_classical_register("cr", nqubits)

    circ_qregs = [qreg]
    circ_cregs = [creg]

    qubits = getqubitlist(qreg)

    if withancilla:
        qareg = qprog.create_quantum_register("qar", nqubits-3) # ancilla registers
        circ_qregs.append(qareg)

        if showancilla:
            careg = qprog.create_classical_register("car", nqubits-3)
            circ_cregs.append(careg)

        ancillaqubits = getqubitlist(qareg)
    else:
        ancillaqubits = None

    qcirc = qprog.create_circuit("qc", circ_qregs, circ_cregs)

    #looprange = range(0,nqubits)

    # there might be an issue with defining numiterations this way, since nqubits
    # may be greater than the length of the target bitstring. also there will be
    # more solutions, since we get superpositions of all values of extra qubits

    # construct initial state - uniform superposition
    allqubitsH(qcirc,qubits)

    # Iterate
    for i in range(0,numiterations):
        # print("Grover iteration ",i)
        # the oracle
        oracle(qcirc,qubits,ancillaqubits,target_bitstrings=target_bitstrings)
        # Grover groverdiffusion
        groverdiffusion(qcirc,qubits,ancillaqubits)

    qcirc.measure(qreg,creg)
    if withancilla and showancilla:
    # sanity check that ancillas are zero
        qcirc.measure(qareg,careg)
    res=qprog.execute(["qc"],backend,shots=shots)

    executiontime = time.time() - starttime

    # add extra info to the results object
    res.numiterations = numiterations
    res.groverexecutiontime = executiontime
    counts=res.get_counts("qc")

    print("Grover's algorithm executed in ", timestring(executiontime), " Used %s,"
    % backend.replace("_"," ") , " on ", nqubits, " qubits, over ", shots,
    " shots. With ancilla = ",withancilla, ".\nResults:")
    print(counts)

    return res #counts
