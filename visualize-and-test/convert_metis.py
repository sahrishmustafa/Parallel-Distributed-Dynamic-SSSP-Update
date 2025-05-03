
import sys
from collections import defaultdict

def convert_to_metis(input_file, output_file):
    # Read the edge list and build adjacency lists
    adjacency = defaultdict(list)
    nodes = set()
    
    with open(input_file, 'r') as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            u, v = map(int, line.strip().split())
            nodes.add(u)
            nodes.add(v)
            adjacency[u].append(v)
            adjacency[v].append(u)  # Assuming undirected graph

    # Sort nodes and create mapping to consecutive integers
    sorted_nodes = sorted(nodes)
    node_map = {node: i+1 for i, node in enumerate(sorted_nodes)}  # METIS uses 1-based
    
    # Count edges (undirected, so each edge appears twice in adjacency lists)
    edge_count = sum(len(neighbors) for neighbors in adjacency.values()) // 2

    # Write METIS format
    with open(output_file, 'w') as f:
        # Header: nodes edges [fmt=0] [ncon=1]
        f.write(f"{len(nodes)} {edge_count} 0 1\n")
        
        # Write adjacency lists
        for node in sorted_nodes:
            neighbors = sorted(adjacency[node])
            # Convert to mapped node IDs
            mapped_neighbors = [node_map[v] for v in neighbors]
            # Write: node_id neighbor1 neighbor2 ...
            f.write(f"{node_map[node]} {' '.join(map(str, mapped_neighbors))}\n")

if __name__ == "__main__":

    import os

    dataset_path = os.path.join("..", "dataset", "facebook", "facebook_combined.txt")
    output_path = os.path.join("..", "dataset", "facebook")
    metis_output = output_path + "\\facebook_combined.graph"

    convert_to_metis(dataset_path, metis_output)
  