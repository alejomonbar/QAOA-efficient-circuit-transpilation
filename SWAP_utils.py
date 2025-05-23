import numpy as np
from qiskit import QuantumCircuit
from copy import copy

def SWAP_pairs(nq):
    qubit_order = list(range(nq))
    list_2q = [[(qubit_order[ii],qubit_order[ii+1]) for ii in range(0, nq-1, 2)]]
    for i in range(1, nq):
        for j in range(i % 2, nq-1, 2):
            qubit_order[j], qubit_order[j+1] = qubit_order[j+1], qubit_order[j]
        list_2q.append([(qubit_order[ii],qubit_order[ii+1]) for ii in range(i%2, nq-1, 2)])
    return list_2q

def qubit_order(nq, depth):
    qubit_order = list(range(nq))
    for i in range(1, depth):
        for j in range(i % 2, nq-1, 2):
            qubit_order[j], qubit_order[j+1] = qubit_order[j+1], qubit_order[j]
    return qubit_order

def layer_1D_full_Graph(gamma, G, list_parity, pi):
    depth = len(list_parity) - 1 
    num_qubits = G.number_of_nodes()
    max_weight = np.max(np.abs([G[i][j]["weight"] for i, j in G.edges()]))
    HC = QuantumCircuit(num_qubits) # Add a layer to the swap network
    
    first_layer = list_parity[0]
    for nn, (i, j) in enumerate(first_layer):
            dd = 2*nn + (depth % 2) * (pi%2)
            if dd < num_qubits-1:             
                if (i,j) in G.edges():
                    HC.rzz(2 * G[i][j]["weight"] * gamma / max_weight, dd, dd+1)
    for kk, layer in enumerate(list_parity[1:-1]):
        for nn, (i, j) in enumerate(layer):
            # Determine which matrix element is required from the current permutation
            # Add the ZZ evolution gate with this matrix element
            dd = 2*nn + ((kk + 1 + (depth%2) * (pi%2)) % 2) 
            if dd < num_qubits-1:            
                HC.cx(dd, dd+1)
                if (i,j) in G.edges():
                    HC.rz(2 * G[i][j]["weight"] * gamma / max_weight, dd+1)
                HC.cx(dd+1, dd)
                HC.cx(dd, dd+1)
    last_layer = list_parity[-1]
    for nn, (i, j) in enumerate(last_layer):
            dd = 2*nn + (depth % 2) * (1 - pi % 2)
            if dd < num_qubits-1:             
                if (i,j) in G.edges():
                    HC.rzz(2 * G[i][j]["weight"] * gamma / max_weight, dd, dd+1)
    return HC

def qaoa_SWAP(gammas, betas, G, list_parity):
    """
    Function that creates a QAOA algorithm with a full swap network for a 1D architecture.
    This is useful if the QUBO problem contains combinations between all variables,
    i.e., all possible qubit pairs need to be connected.
    """
    #Method to design the quantum circuit 

    #obtain the values for our qubo input
    #there are the quadratic and linear coeffs.

    num_qubits = G.number_of_nodes()
    p = len(gammas)
    circ = QuantumCircuit(num_qubits)
    #Apply the initial layer of Hadamard gates to all qubits
    for i in range(num_qubits):
        circ.h(i)
    #Outer loop to create each layer
    for pi in range(p):
        layer = layer_1D_full_Graph(gammas[pi], G, list_parity, pi)
        circ = circ.compose(layer, range(num_qubits))
        circ.rx(-2 * betas[pi], range(num_qubits))
        list_parity = list(reversed(copy(list_parity)))
       
    return circ 