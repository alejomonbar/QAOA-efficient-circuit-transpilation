import random
import math
import networkx as nx


def simulated_annealing_func(G_original: nx.Graph, initial_mapping: dict, list_2q: list, callback: dict,
                             initial_temp=0.01, cooling_rate=0.9999,
                             stop_temp=1e-8, max_iter=10000, verbose=False, max_reheat=5):
    n = G_original.number_of_nodes()
    max_iter = int(max_iter)

    # Use list for O(1) index access (faster than dict)
    mapping = [initial_mapping[k] for k in range(n)]

    # Adjacency matrix for O(1) edge lookup (avoids tuple creation + hashing)
    adj = [[False] * n for _ in range(n)]
    for u, v in G_original.edges():
        adj[u][v] = True
        adj[v][u] = True

    num_layers = len(list_2q)

    # Precompute: for each node, list of (layer_index, partner_node) for all 2q edges
    node_edges = [[] for _ in range(n)]
    for l_idx in range(num_layers):
        for edge in list_2q[l_idx]:
            if len(edge) == 2:
                u, v = edge
                node_edges[u].append((l_idx, v))
                node_edges[v].append((l_idx, u))

    # Per-layer valid and invalid 2q-edge counts
    valid_count = [0] * num_layers
    invalid_count = [0] * num_layers

    def recompute_layer_counts():
        for l_idx in range(num_layers):
            vc = 0
            ic = 0
            for edge in list_2q[l_idx]:
                if len(edge) == 2:
                    u, v = edge
                    if adj[mapping[u]][mapping[v]]:
                        vc += 1
                    else:
                        ic += 1
            valid_count[l_idx] = vc
            invalid_count[l_idx] = ic

    recompute_layer_counts()

    def compute_cost_depth():
        cnots = 0
        depth = 0
        for l in range(num_layers - 1, -1, -1):
            cnots -= invalid_count[l]
            if valid_count[l] > 0:
                return cnots, depth
            depth -= 1
        return cnots, depth

    current_cost, current_depth = compute_cost_depth()
    best_mapping = mapping[:]
    best_cost = current_cost
    best_depth = current_depth

    T = initial_temp
    iteration = 0
    current_iter = callback["iterations"][-1] if callback["iterations"] else 0

    # Local references for hot-path speedup
    _randint = random.randint
    _random = random.random
    _exp = math.exp
    n_minus_1 = n - 1
    n_minus_2 = n - 2

    for _ in range(max_reheat):
        T = initial_temp
        iteration = 0
        current_iter = callback["iterations"][-1] if callback["iterations"] else 0

        while T > stop_temp and iteration < max_iter:
            # Fast random pair selection (avoids range() + sample())
            i = _randint(0, n_minus_1)
            j = _randint(0, n_minus_2)
            if j >= i:
                j += 1

            phys_i = mapping[i]
            phys_j = mapping[j]

            # Incremental cost: only check edges involving swapped nodes i or j
            changes = []

            for l_idx, partner in node_edges[i]:
                if partner == j:
                    continue  # edge (i,j): swapping both endpoints doesn't change validity
                mp = mapping[partner]
                old_v = adj[phys_i][mp]
                new_v = adj[phys_j][mp]
                if old_v != new_v:
                    changes.append((l_idx, 1 if new_v else -1))

            for l_idx, partner in node_edges[j]:
                if partner == i:
                    continue
                mp = mapping[partner]
                old_v = adj[phys_j][mp]
                new_v = adj[phys_i][mp]
                if old_v != new_v:
                    changes.append((l_idx, 1 if new_v else -1))

            # Apply incremental changes
            for l_idx, delta in changes:
                valid_count[l_idx] += delta
                invalid_count[l_idx] -= delta

            neighbor_cost, neighbor_depth = compute_cost_depth()
            delta_cost = neighbor_cost - current_cost

            # Accept or reject the new solution
            if delta_cost < 0 or _random() < _exp(-delta_cost / T):
                mapping[i] = phys_j
                mapping[j] = phys_i
                current_cost = neighbor_cost
                current_depth = neighbor_depth

                if current_cost < best_cost:
                    best_mapping = mapping[:]
                    best_cost = current_cost
                    best_depth = current_depth

                    if verbose:
                        print(f"best depth:{best_depth} | best cost:{best_cost} | iteration:{iteration} | T:{T:.4f}")

                    callback["cost"].append(best_cost)
                    callback["depth"].append(best_depth)
                    callback["iterations"].append(current_iter + iteration)
                    callback["T"].append(T)
            else:
                # Reject: revert incremental changes
                for l_idx, delta in changes:
                    valid_count[l_idx] -= delta
                    invalid_count[l_idx] += delta

            T *= cooling_rate
            iteration += 1

    best_mapping_dict = dict(enumerate(best_mapping))
    return best_mapping_dict, best_cost, callback