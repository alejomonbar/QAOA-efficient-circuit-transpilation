import random
import math
import networkx as nx

def simulated_annealing_func(G_original:nx.Graph, initial_mapping:dict, list_2q:list, callback:dict, 
                                  initial_temp=0.01, cooling_rate=0.9999,
                                  stop_temp=1e-8, max_iter=10000, verbose=False):
    n = G_original.number_of_nodes()
    mapping = initial_mapping.copy()

    edge_set = set()
    for u, v in G_original.edges():
        edge_set.add((min(u, v), max(u, v)))
    
    def fast_total_cnots_depth(list_2q, edge_set):
        cnots = 0
        depth = 0
        for layer in list_2q[::-1]:
            layer_has_valid_edge = False
            for edge in layer:
                if len(edge) == 2:
                    u, v = edge
                    mapped_u, mapped_v = mapping[u], mapping[v]
                    
                    mapped_edge = (min(mapped_u, mapped_v), max(mapped_u, mapped_v))
                    if mapped_edge not in edge_set:
                        cnots -= 1
                    else:
                        layer_has_valid_edge = True
            if layer_has_valid_edge:
                return cnots, depth
            
            depth -= 1
        
        return cnots, depth
    current_cost, current_depth = fast_total_cnots_depth(list_2q, edge_set)
    best_mapping = mapping.copy()
    best_cost = current_cost
    best_depth = current_depth
    
    T = initial_temp
    iteration = 0
    current_iter = callback["iterations"][-1] if callback["iterations"] else 0
    
    while T > stop_temp and iteration < max_iter:
        # Select two random nodes to swap
        i, j = random.sample(range(n), 2)
        
        # Swap the mapping
        mapping[i], mapping[j] = mapping[j], mapping[i] 

        neighbor_cost, neighbor_depth = fast_total_cnots_depth(list_2q, edge_set)
        delta = neighbor_cost - current_cost
        
        # Accept or reject the new solution
        if delta < 0 or random.random() < math.exp(-delta / T):
            current_cost = neighbor_cost
            current_depth = neighbor_depth
            
            if current_cost < best_cost:
                best_mapping = mapping.copy()
                best_cost = current_cost
                best_depth = current_depth
                
                if verbose:
                    print(f"best depth:{best_depth} | best cost:{best_cost} | iteration:{iteration} | T:{T:.4f}")
                
                callback["cost"].append(best_cost)
                callback["depth"].append(best_depth)
                callback["iterations"].append(current_iter + iteration)
                callback["T"].append(T)
        else:
            # Revert the swap
            mapping[i], mapping[j] = mapping[j], mapping[i]
        T *= cooling_rate
        iteration += 1
    return best_mapping, best_cost, callback
