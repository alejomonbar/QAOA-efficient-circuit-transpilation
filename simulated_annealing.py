import random
import math
import networkx as nx


def simulated_annealing_func(G_original: nx.Graph, initial_mapping: dict, list_2q: list, callback: dict,
                             initial_temp=0.01, cooling_rate=0.9999,
                             stop_temp=1e-8, max_iter=10000, verbose=False,
                             stagnation_limit=None, perturbation_size=None):
    n = G_original.number_of_nodes()
    max_iter = int(max_iter)

    # Stagnation defaults: reheat after 5% of max_iter without improvement,
    # kick with n//5 random swaps to escape the basin
    if stagnation_limit is None:
        stagnation_limit = max(1000, max_iter // 20)
    if perturbation_size is None:
        perturbation_size = max(2, n // 5)

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
    iters_since_improvement = 0
    current_iter = callback["iterations"][-1] if callback["iterations"] else 0

    # Local references for hot-path speedup
    _randint = random.randint
    _random = random.random
    _exp = math.exp
    _sample = random.sample
    n_minus_1 = n - 1
    n_minus_2 = n - 2
    all_nodes = list(range(n))

    while T > stop_temp and iteration < max_iter:
        # --- Stagnation detection: reheat + perturbation ---
        if iters_since_improvement >= stagnation_limit:
            # Restart from best known solution
            mapping[:] = best_mapping
            # Random multi-swap perturbation to escape the basin
            perturb_nodes = _sample(all_nodes, min(perturbation_size, n))
            for k in range(len(perturb_nodes) - 1):
                a, b = perturb_nodes[k], perturb_nodes[k + 1]
                mapping[a], mapping[b] = mapping[b], mapping[a]
            # Recompute counts from scratch after large perturbation
            recompute_layer_counts()
            current_cost, current_depth = compute_cost_depth()
            # Reheat
            T = initial_temp
            iters_since_improvement = 0
            if verbose:
                print(f"  reheat @ iter {iteration} | perturbed {len(perturb_nodes)} nodes | cost after kick: {current_cost}")

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
                iters_since_improvement = 0

                if verbose:
                    print(f"best depth:{best_depth} | best cost:{best_cost} | iteration:{iteration} | T:{T:.4f}")

                callback["cost"].append(best_cost)
                callback["depth"].append(best_depth)
                callback["iterations"].append(current_iter + iteration)
                callback["T"].append(T)
            else:
                iters_since_improvement += 1
        else:
            # Reject: revert incremental changes
            for l_idx, delta in changes:
                valid_count[l_idx] -= delta
                invalid_count[l_idx] += delta
            iters_since_improvement += 1

        T *= cooling_rate
        iteration += 1

    best_mapping_dict = dict(enumerate(best_mapping))
    return best_mapping_dict, best_cost, callback
