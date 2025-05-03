import random
import os
import networkx as nx

NUM_INSERT = 0
NUM_DEL = 0

def load_edge_list(path):
    G = nx.read_edgelist(path, nodetype=int, create_using=nx.DiGraph())
    return G

def generate_insertions(G, num_insertions):
    insertions = []
    nodes = list(G.nodes)
    existing = set(G.edges)

    while len(insertions) < num_insertions:
        u = random.choice(nodes)
        v = random.choice(nodes)
        if u != v and (u, v) not in existing:
            insertions.append((u, v))
            existing.add((u, v))  # avoid duplicates

    return insertions

def generate_deletions(G, num_deletions):
    deletions = random.sample(list(G.edges), min(num_deletions, G.number_of_edges()))
    return deletions

def write_edges(filename, edges):
    with open(filename, 'w') as f:
        for u, v in edges:
            f.write(f"{u} {v}\n")

if __name__ == "__main__":

    # Get current script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))

    graph_path = os.path.join(script_dir, '..', 'dataset', 'facebook', 'facebook_combined.txt')
    graph_dir = os.path.join(script_dir, '..', 'dataset', 'facebook')
    G = load_edge_list(graph_path)

    num_insertions = NUM_INSERT
    num_deletions = NUM_DEL

    insertions = generate_insertions(G, num_insertions)
    deletions = generate_deletions(G, num_deletions)

    ins_path = graph_dir + "\\insertions.txt"
    write_edges(ins_path, insertions)
    del_path = graph_dir + "\\deletions.txt"
    write_edges(del_path, deletions)

    print(f"✔️ Generated {len(insertions)} insertions and {len(deletions)} deletions.")
