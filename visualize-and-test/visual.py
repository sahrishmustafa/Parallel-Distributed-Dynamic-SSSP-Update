import networkx as nx
import matplotlib.pyplot as plt
import os
from math import inf

def read_graph_data(graph_path, insertions_path, deletions_path):
    # Read initial graph edges
    with open(graph_path, 'r') as f:
        initial_edges = [tuple(map(int, line.strip().split())) for line in f if line.strip()]

    # Read insertions
    with open(insertions_path, 'r') as f:
        insertions = [tuple(map(int, line.strip().split())) for line in f if line.strip()]

    # Read deletions
    with open(deletions_path, 'r') as f:
        deletions = [tuple(map(int, line.strip().split())) for line in f if line.strip()]
    
    return initial_edges, insertions, deletions


def run_dijkstra(G, source):
    # Initialize distances and predecessors
    distances = {node: inf for node in G.nodes()}
    predecessors = {node: None for node in G.nodes()}
    distances[source] = 0
    
    # Priority queue
    unvisited = set(G.nodes())
    
    while unvisited:
        # Get node with smallest distance
        current = min(unvisited, key=lambda node: distances[node])
        unvisited.remove(current)
        
        # Stop if remaining nodes are unreachable
        if distances[current] == inf:
            break
            
        for neighbor in G.neighbors(current):
            # Assuming all edges have weight=1 (unweighted graph)
            # For weighted graphs, you'd use G[current][neighbor]['weight']
            alt = distances[current] + 1  
            if alt < distances[neighbor]:
                distances[neighbor] = alt
                predecessors[neighbor] = current
    
    return distances, predecessors

def get_shortest_path_edges(predecessors):
    edges = []
    for node, pred in predecessors.items():
        if pred is not None:
            edges.append((pred, node))
    return edges

def save_dijkstra_results(distances, predecessors, output_path):
    with open(output_path, 'w') as f:
        for node in sorted(distances.keys()):
            dist = distances[node]
            parent = predecessors[node] if predecessors[node] is not None else "-1"
            f.write(f"{node}\t{dist}\t{parent}\n")

def visualize_shortest_path(G, distances, predecessors, source, pos):
    plt.figure(figsize=(14, 12))
    
    # Get all edges in the shortest path tree
    path_edges = get_shortest_path_edges(predecessors)
    
    # Color nodes based on distance
    max_dist = max(d for d in distances.values() if d != inf)
    node_colors = []
    for node in G.nodes():
        if distances[node] == inf:
            node_colors.append('red')  # Unreachable
        else:
            # Color gradient from lightblue (close) to darkblue (far)
            intensity = 0.3 + 0.7 * (distances[node] / max_dist)
            node_colors.append((0.2, 0.4, intensity))
    
    # Draw all nodes with distance-based coloring
    nx.draw_networkx_nodes(G, pos, node_size=80, node_color=node_colors)
    
    # Draw all edges (light gray, faint)
    nx.draw_networkx_edges(G, pos, edge_color='lightgray', width=0.5, alpha=0.3)
    
    # Highlight shortest path edges (bold blue)
    nx.draw_networkx_edges(G, pos, edgelist=path_edges, 
                          edge_color='blue', width=2, alpha=0.7)
    
    # Highlight source node
    nx.draw_networkx_nodes(G, pos, nodelist=[source], 
                         node_size=200, node_color='yellow')
    
    # Label nodes with their distances
    labels = {node: f"{node}\n({distances[node]})" for node in G.nodes() 
             if distances[node] != inf}  # Label every 3rd node
    nx.draw_networkx_labels(G, pos, labels, font_size=8)
    
    plt.title(f"Shortest Path Tree from Node {source}\n"
             f"(Yellow = Source, Color = Distance, Blue = Shortest Path Edges)")
    plt.axis('off')
    plt.tight_layout()
    return plt

if __name__ == "__main__":
    # Configuration
    SOURCE_NODE = 2  # Change this to your desired source node
    
    # Get current script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Construct paths
    graph_path = os.path.join(script_dir, '..', 'dataset', 'dataset1', 'graph.txt')
    insertions_path = os.path.join(script_dir, '..', 'dataset', 'dataset1', 'insertions.txt')
    deletions_path = os.path.join(script_dir, '..', 'dataset', 'dataset1', 'deletions.txt')
    output_dir = os.path.join(script_dir, '..','testing_visuals', 'output')
    os.makedirs(output_dir, exist_ok=True)
    
    # Read graph data
    initial_edges, insertions, deletions = read_graph_data(graph_path, insertions_path, deletions_path)
    
    # Create graph and add initial edges
    G = nx.Graph()
    G.add_edges_from(initial_edges)
    
    # Remove deleted edges
    for edge in deletions:
        if G.has_edge(*edge):
            G.remove_edge(*edge)
    
    # Add insertions
    if insertions:
        G.add_edges_from(insertions)
    
    # Run Dijkstra's algorithm
    distances, predecessors = run_dijkstra(G, SOURCE_NODE)
    
    # Save results to file
    results_path = os.path.join(output_dir, 'output1.txt')
    save_dijkstra_results(distances, predecessors, results_path)
    
    # Generate visualization
    # pos = nx.spring_layout(G, seed=42)  # Consistent layout
    # plt = visualize_shortest_path(G, distances, predecessors, SOURCE_NODE, pos)
    
    # # Save and show visualization
    # visualization_path = os.path.join(output_dir, 'shortest_path_tree.png')
    # plt.savefig(visualization_path, dpi=300)
    #plt.show()
    
    #print(f"Results saved to:\n- {results_path}\n- {visualization_path}")
