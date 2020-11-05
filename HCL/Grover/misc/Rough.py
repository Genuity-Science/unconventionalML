import numpy as np
import math
import random
import operator
# IBM QISkit libraries
from qiskit import QuantumProgram
import qtools.Qconfig as Qconfig
# our libraries
from qtools.unitary import u4powerparams, u4conjparams
#from qtools.gates import *
#from qtools.qtools import *

# Check inter-register interaction

backendcode=0
backendarray = ['local_qiskit_simulator','local_qasm_simulator','ibmqx_qasm_simulator','ibmqx_hpc_qasm_simulator', 'ibmqx5']
backend = backendarray[backendcode]

n=3

qprog = QuantumProgram()
qprog.set_api(Qconfig.APItoken,Qconfig.config['url'])
qreg = qprog.create_quantum_register("qr", n)
areg = qprog.create_quantum_register("ar", n)
creg = qprog.create_classical_register("cr", n)
careg = qprog.create_classical_register("car", n)

circ = qprog.create_circuit("qc", [qreg,areg], [creg,careg])
qcirc = qprog.create_circuit("qc", [qreg,areg], [creg,careg])
qcirc.x(qreg[0])
qcirc.h(qreg[1])

for i in range(0,n):
    qcirc.cx(qreg[i], areg[i])
    qcirc.measure(qreg[i],creg[i])
    qcirc.measure(areg[i],careg[i])

res=qprog.execute(["qc"],backend)
res.get_counts("qc")

# number of datapoints possible for each number of qubits (will need distance between datapoints)
x=[(i,round(2**(i/2),3)) for i in range(0,21)]