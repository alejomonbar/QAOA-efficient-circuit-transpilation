from docplex.mp.model import Model
import numpy as np
from collections import defaultdict
from qiskit import QuantumCircuit

def random_samples(num_samples, n_qubits):
    random_samples = defaultdict(int)
    for i in range(num_samples):
        random_samples["".join(str(i) for i in np.random.choice([0,1], n_qubits))] += 1
    return random_samples


def WMaxCut(G):
    # MIS model as a QUBO problem
    mdl = Model('MaxCut')
    num_vertices = G.number_of_nodes()
    x = {i: mdl.binary_var(name=f"x_{i}") for i in range(num_vertices)}
    mdl.minimize(-mdl.sum(G[i][j]["weight"]*(x[i] + x[j] - 2 * x[i] * x[j]) for i, j in G.edges))
    return mdl

def cost_maxcut(bitstring, G):
    cost = 0
    for i, j in G.edges():
        if bitstring[i] + bitstring[j] in ["10","01"]:
            cost += G[i][j]["weight"]
    return cost


def objective_MaxCut(samples_dict, G, optimal):
    max_cost = cost_maxcut(optimal, G)
    results = []
    probability = 0
    for bitstring, counts in samples_dict.items():
        cost = cost_maxcut(bitstring, G)
        r  = cost/max_cost
        results.append([cost, r, counts])
        if abs(cost - max_cost) < 1e-4:
            probability += counts
        if cost > max_cost:
            print(f"There is a better cost than that of CPLEX: {cost - max_cost}")
            print(f" optimal: {optimal}")
            print(f" current: {bitstring}")
    results = np.array(results)
    shots = np.sum(results[:,2])
    rT = np.sum(results[:,0] * results[:,2])/(shots*max_cost)
    probability /= shots
    return {"results":np.array(results), "max_cut":max_cost, "r":rT, "probability":probability}

def mitigate(samples_dict, G, random=False):
    change = {"0":"1", "1":"0"}
    nq = G.number_of_nodes()
    new_samples = defaultdict(int)
    for bitstring, counts in samples_dict.items():
            for _ in range(counts):
                best_string = bitstring
                best_cost = cost_maxcut(bitstring, G)
                list_qubits = np.arange(nq)
                if random:
                    np.random.shuffle(list_qubits)
                for qi in list_qubits:
                    new_string = "".join((change[i] if n == qi else i) for n, i in enumerate(best_string))
                    new_cost = cost_maxcut(new_string, G)
                    if new_cost > best_cost:
                        best_string = new_string
                        best_cost = new_cost
                new_samples[best_string] += 1
    return new_samples

def qaoa(gammas, betas,  G):
    nodes = len(G.nodes())
    edges = [(i,j) for i, j in G.edges()]
    max_w = np.max(np.abs([G[i][j]["weight"] for i, j in G.edges()]))
    layers = len(gammas)
    qc = QuantumCircuit(nodes)
    qc.h(range(nodes))
    for p in range(layers):
        permutations = np.arange(nodes) # To decrease the depth of the circuit
        for jj in range(nodes):
            for k in range(jj % 2, nodes - 1, 2):
                qubit_pair = (permutations[k], permutations[k+1])
                if qubit_pair in edges or reversed(qubit_pair) in edges:
                    qc.rzz(2 * gammas[p] * G[qubit_pair[0]][qubit_pair[1]]["weight"]/max_w, *qubit_pair)
                permutations[[k, k+1]] = permutations[[k+1, k]]
        qc.rx(-2*betas[p], range(nodes))
    return qc

def total_cnots_depth(G, list_2q):
    cnots = 0
    depth = 0
    flag = False
    for layer in list_2q[::-1]:
        for edge in layer:
            if edge not in G.edges:
                cnots -= 1
            else:
                flag = True
        if flag:
            return cnots, depth
        depth -= 1
    return cnots, depth