from qiskit import QuantumCircuit
import numpy as np
from copy import copy

def PTC_pairs(nq):
    list_2q = [nq * [0]]
    list_2q[0][0] = (0,)
    for j in range(1, nq):
        if len(list_2q[0][j-1])==1: 
            if j % 2:
                list_2q[0][j] = (j,)
            else:
                list_2q[0][j] = (list_2q[0][j-1][0], nq-1)
        else:
            list_2q[0][j] = (list_2q[0][j-1][0]+j % 2, list_2q[0][j-1][1] - 1 + j % 2)
    
    for i in range(1, nq):
        array_i = copy(list_2q[i-1])
        if i % 2 == 1:
            lista =  list(list_2q[i-1][0]) +  list(list_2q[i-1][1]) 
            unique_list = tuple(sorted([item for item in lista if lista.count(item) == 1]))
            array_i[0] = unique_list
        for j in range(i % 2, nq-2, 2):
            lista =  list(list_2q[i-1][j]) +  list(list_2q[i-1][j+1]) + list(list_2q[i-1][j+2])
            unique_list = tuple(sorted([item for item in lista if lista.count(item) == 1]))
            array_i[j+1] = unique_list
        if nq % 2 == i % 2:
            lista =  list(list_2q[i-1][nq-1]) +  list(list_2q[i-1][nq-2]) 
            unique_list = tuple(sorted([item for item in lista if lista.count(item) == 1]))
            array_i[nq-1] = unique_list
        list_2q.append(array_i)
    return list_2q



def UA(G, gamma, sequence_2q, pi):
    nq = G.number_of_nodes()
    qc = QuantumCircuit(nq)
    max_weight = np.max(np.abs([G[i][j]["weight"] for i, j in G.edges()]))
    depth = len(sequence_2q)
    for k in range(depth):
        for qi, ij in enumerate(sequence_2q[k]):
            if ij in G.edges() or (ij[0],ij[0]) in G.edges():
                if k > 0 and ij == sequence_2q[k-1][qi]:
                    pass
                elif len(ij) == 2:
                    if ij[0] in list(G[ij[1]].keys()):
                        qc.rz(2 * G[ij[0]][ij[1]]["weight"] * gamma / max_weight, qi)
                elif len(ij) == 1:
                    if ij[0] in list(G[ij[0]].keys()):
                        qc.rz(2 * G[ij[0]][ij[0]]["weight"] * gamma / max_weight, qi)
        if k < depth - 1:
            for i in range((k + (depth % 2) * pi) % 2, nq-1, 2):
                qc.cx(i+1, i)
            for j in range(((k + 1) + (depth % 2) * pi) % 2, nq-1, 2):
                qc.cx(j, j+1)
    return qc

def UB(qc, nq, beta, sequence_2q):
    qc.rx(-2*beta, 0)
    last_sequence = sequence_2q[-1]
    for i in range(2, nq, 2):
        for j in last_sequence[i]:
            if j in last_sequence[i-1]:
                qc.cx(i, i-1)
                qc.rx(-2 * beta, i)
                qc.cx(i, i-1)
            elif i < nq - 1 and j not in last_sequence[i+1]:
                qc.rx(-2*beta, i)
    for i in range(1, nq, 2):
        for j in last_sequence[i]:
            if j in last_sequence[i-1]:
                qc.cx(i, i-1)
                qc.rx(-2 * beta, i)
                qc.cx(i, i-1)
            elif i < nq - 1 and j not in last_sequence[i+1]:
                qc.rx(-2*beta, i)
    qc.rx(-2*beta, nq-1)
    

def qaoa_PTC(gammas, betas, G, list_parity):
    nq = G.number_of_nodes()
    p = len(gammas)
    qc = QuantumCircuit(nq)
    qc.h(range(nq))
    cost_layer_0 = lambda gamma: UA(G, gamma, copy(list_parity), 0)
    cost_layer_1 = lambda gamma: UA(G, gamma, copy(list_parity), 1)
    for pi in range(p):
        if pi % 2 == 0:
            qc.compose(cost_layer_0(gammas[pi]), inplace=True)
        else:
            qc.compose(cost_layer_1(gammas[pi]), inplace=True)
        UB(qc, nq, betas[pi], list_parity)
        list_parity = list(reversed(copy(list_parity)))
    return qc 

def decode_parity(sample, p, list_parity):
    depth = len(list_parity)
    spin = {1:"0", -1:"1"}
    sz = [(1 if i == "0" else -1) for i in sample]
    nq = len(sample)
    x = nq * [0]
    if p % 2 == 0:
        list_parity = list_parity[0]
        list_q = range(nq)
    else:
        list_parity = list_parity[-1]
        list_q = list(reversed(range(depth))) + list(range(depth, nq))
    for i in list_q:
        if len(list_parity[i]) == 1:
            x[list_parity[i][0]] = sz[i]
        else:
            if x[list_parity[i][0]] != 0:
                x[list_parity[i][1]] = sz[i]/x[list_parity[i][0]]
            elif x[list_parity[i][1]] != 0:
                x[list_parity[i][0]] = sz[i]/x[list_parity[i][1]]
            else:
                print(list_parity[i])
    return "".join(spin[xi] for xi in x)
