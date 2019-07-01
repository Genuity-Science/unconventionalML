# script that initializes a quantum circuit for toy test purposes
# Jan 2018

import numpy as np
import math

# IBM QISkit libraries
from qiskit import QuantumProgram
import Qconfig

# our library
from tools.unitary import u4powerparams, u4conjparams
from tools.gates import *
from tools.qtools import *
backendarray=['local_qasm_simulator','ibmqx_qasm_simulator', 'ibmqx5']


nbits=4
backendcode=1

backend = backendarray[backendcode]
qprog = QuantumProgram()
qprog.set_api(Qconfig.APItoken, Qconfig.config['url'])
qreg = qprog.create_quantum_register("qr", nbits)
creg = qprog.create_classical_register("cr", nbits)
qcirc = qprog.create_circuit("qc", [qreg], [creg])

qubits = getqubitlist(qreg)
