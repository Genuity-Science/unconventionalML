# Based on Hooer Durr Algorithm
# native librlm,aries
import numpy as np
import math
import random
import operator
import time

# IBM QISkit libraries
from qiskit import QuantumProgram
import quantumtools.Qconfig as Qconfig

# our libraries
from quantumtools.unitary import u4powerparams, u4conjparams
from quantumtools.gates import *
from quantumtools.qtools import *

from Grover.Grover import groverdiffusion

# backendcode: 0: local simulator, 1: ibmqx simulator, 2: ibm¡qx5
backendarray = ['local_qasm_simulator_cpp','ibmqx5','local_statevector_simulator','local_qasm_simulator','ibmqx_qasm_simulator','ibmqx_hpc_qasm_simulator']

fsin = lambda x: np.round(np.sin(10*x),3) #f(x)... arbitrary oscillating function
fvalues = list(map(fsin,range(0,2**6)))

def num_qubits_from_length(length):
    return int(np.ceil(np.log2(length)))


def opt_oracle_lessthan(qcirc, qubits, last_index,ancillaqubits=None,fvalues=fvalues):
    """ execute an artificial oracle for grover's minimization algorithm on qcirc, i.e.
    it flips the phase for states x that satisfy f(x) < f(y), for f defined in
    this function. (long run, f should be a parameter)
    y is an index (input)
    -in a real quantum implementation, this would (should) have some natural and simple circuit that doesn't
    use knowledge of f ahead of time.
    """

    # for this to be practical, the oracle must be created efficiently from basic
    # knowledge
    nstringbits = num_qubits_from_length(len(fvalues))

    target_bitstrings = [bin(i)[2:].zfill(nstringbits) for i,v in enumerate(fvalues) if v<fvalues[last_index]]
    #print("target_bitstrings: ",target_bitstrings)

    nqubits = len(qubits) #make sure we have enough qubits for the string
    assert nqubits >= nstringbits, "bitstring too long for " + str(nqubits) + " qubits"

    # apply ctrl-n-z for each string as the control - like oracle_mark, need to make it
    # more efficient
    for bitstring in target_bitstrings:
        # controlled z gate is symmetric
        # for ctrl-n-z, choice of target qubit is arbitrary, so use first qubit as target
        cvbitstring=bitstring[0:-1]
        tgt = int(bitstring[-1]) # first qubit is last entry in bitstring

        if not tgt: # flip the target qubit if its value in bitstring (control) is 0
            qcirc.x(qubits[0])

        # apply the controlled phase (only changes what matches bitstring
        # cnu4(qcirc, parz, qubits[1:nstringbits], qubits[0],ancillaqubits,cvbitstring)
        qcirc.h(qubits[0])
        cnx(qcirc, qubits[1:nstringbits], qubits[0],ancillaqubits,cvbitstring)
        qcirc.h(qubits[0])

        if not tgt: # flip the target qubit back if its value in bitstring (control) is 0
            qcirc.x(qubits[0])
    return

def groverminimize_split(opt_oracle,bits_reduce=2,backendcode=0, withancilla=True, fvalues=fvalues, silent=False):
    """recursively split the list of values fvalues to several and minimize over them
    qubits_reduce is number of qubits by which to drop the representation each iteration
    i.e. divide and conquer """
    numvalues=len(fvalues)
    nbits=np.ceil(np.log2(numvalues))-bits_reduce
    # print("numvalues ",numvalues)
    # print("fvalues ",fvalues)
    # print("nbits ",nbits)


    if nbits <= 3: # base case, minimum 3 qubits
        return groverminimize(opt_oracle, backendcode=backendcode, withancilla=withancilla, fvalues=fvalues, silent=silent)
    else: # recursive case
        # best to split it so each piece contains a power of two values - to use all qubit coefficients efficiently

        piece_len = int(2 ** nbits)
        npieces = int(np.ceil(numvalues / piece_len))  # could be 2**bits_reduce or 2**bits_reduce-1

        # index of minimum item in each piece (index relative to overall fvalues)
        split_indices = [i*piece_len+groverminimize_split(opt_oracle, backendcode=backendcode, withancilla=withancilla,
                fvalues=fvalues[i*piece_len:int(min((i+1)*piece_len,numvalues))], silent=silent) for i in range(npieces)]
        # for i in range(npieces):
        #     print("fsplit i:",i, " ",fvalues[i*piece_len:int(min((i+1)*piece_len,numvalues))])
        split_values = [fvalues[j] for j in split_indices]

        # minimize between pieces
        relative_ind = groverminimize(opt_oracle,backendcode=0, withancilla=True, fvalues=split_values, silent=silent)
        # print("npieces ", npieces)
        # print("split_indices ",split_indices)
        # print("split_values ",split_values)
        # print("relative_ind ",relative_ind )
        return split_indices[relative_ind]

def groverminimize(opt_oracle,backendcode=0, withancilla=True, fvalues=fvalues, silent=False):#,shots=1024):
    """"
        inputs are the oracle (function), and the number of bits ...
        each iteration it narrows it further
        https://arxiv.org/abs/quant-ph/9607014v2
    """
    showancilla=False
    starttime=time.time()

    nvalues=len(fvalues)

    nqubits = num_qubits_from_length(nvalues)

    if nqubits <= 3:# no need for ancilla if too few qubits
        withancilla=False

    # set up quantum environment. (design): is this the right place to define these?
    backend = backendarray[backendcode]
    qprog = QuantumProgram()
    # qprog.set_api(Qconfig.APItoken,Qconfig.config['url'])

    numOiterations = numGroverOptimizeiterations(nqubits)
    #numGiterations = numGroveriterations(2 ** nqubits, numsolutions)
    optimized_index = random.randint(0,nvalues-1)
    max_repeats = nqubits # number of iterations
    repeats = 0


    for j in range(numOiterations):

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
        # construct initial state - uniform superposition
        allqubitsH(qcirc,qubits)

        # Iterate
        # use very rough heuristic for number of solutions, should be high to start and gets smaller, halves each time ..., improvement TODO
        numsol = np.ceil(2**(nqubits-j*2/3))#numsol=np.ceil(nqubits/(2*sqrt(j+1)))
        n_iter=numGroveriterations(2 ** nqubits, numsol)
        # print("assumed ",numsol, " solutions to calculate ",n_iter, " Grover oracle+diffusion cycles")
        for i in range(0,n_iter):
            # the oracle
            opt_oracle(qcirc,qubits,optimized_index,ancillaqubits,fvalues=fvalues)
            # Grover groverdiffusion
            groverdiffusion(qcirc,qubits,ancillaqubits)

        qcirc.measure(qreg,creg)
        res=qprog.execute(["qc"],backend,timeout=3000)#,shots=1)#024)
        counts=res.get_counts("qc")
        #print(counts)

        # index with maximum counts
        trial_index = int(max(counts.items(), key=operator.itemgetter(1))[0],2)
        # np.random.choice(range(0, 7), p=[0.1, 0.05, 0.05, 0.2, 0.4, 0.2])
        #trial_index=int(list(counts.keys())[0], 2)
        #trial_index = int(max(counts.items(), key=operator.itemgetter(1))[0],2)
        # print("func(trial_index)", func(trial_index), " and func(optimized_index) ",func(optimized_index))

        # print("trial ", trial_index, "opt ind:", optimized_index)
        # get rid of padded indices - use knowledge of fvalues again here (this must be built into oracle)
        if trial_index<nvalues and fvalues[trial_index] < fvalues[optimized_index]:
            # update index
            optimized_index = trial_index
            repeats =0
        else:
            # same index
            repeats+=1

        #print("Iteration ",j, " with optimal index ", optimized_index, " and result:\n",counts) # test output
        # print("Iteration ", j, " with trial index ", trial_index, " and optimal index ", optimized_index)

        if repeats == max_repeats:
            break

    executiontime = time.time() - starttime

    if not silent:
        print("Grover's Optimization algorithm executed over ", j, " iterations in ", timestring(executiontime), " Used %s,"
        % backend.replace("_"," ") , " on ", nqubits, " qubits")#, over ", shots, " shots.")
        print("Optimized index ", optimized_index)

    return optimized_index
