import os
import random
from collections import defaultdict

NUM_INSERT = 10000
NUM_DEL = 5000

def parse_partition_files(partitions_dir):
    """Parse partition files and extract internal edges (skip (X) edges)."""
    partition_edges = defaultdict(list)
    
    for filename in os.listdir(partitions_dir):
        if filename.startswith("part_") and filename.endswith(".txt"):
            part_id = int(filename.split("_")[1].split(".")[0])
            with open(os.path.join(partitions_dir, filename), 'r') as f:
                for line in f:
                    if line.startswith("Node"):
                        # Extract node ID and its edges
                        node_str = line.split(":")[0].strip()
                        node_id = int(node_str.split()[1])
                        edges_str = line.split(":")[1].strip()
                        
                        # Process edges (skip those with (X))
                        for edge in edges_str.split():
                            if "(X)" not in edge:
                                neighbor = int(edge.replace("(X)", ""))
                                partition_edges[part_id].append((node_id, neighbor))
    
    return partition_edges

def generate_updates(partition_edges, num_insertions=10, num_deletions=5):
    """Generate insertions/deletions only between nodes in the same partition."""
    insertions = []
    deletions = []
    
    # Generate insertions (add new edges within partitions)
    for part_id, edges in partition_edges.items():
        nodes = list(set([u for u, v in edges] + [v for u, v in edges]))
        for _ in range(num_insertions // len(partition_edges)):
            u, v = random.sample(nodes, 2)
            if u != v and (u, v) not in edges and (v, u) not in edges:
                insertions.append((u, v))
    
    # Generate deletions (remove existing edges within partitions)
    for part_id, edges in partition_edges.items():
        candidates = [e for e in edges if e not in deletions and (e[1], e[0]) not in deletions]
        for _ in range(min(num_deletions // len(partition_edges), len(candidates))):
            if candidates:
                u, v = random.choice(candidates)
                deletions.append((u, v))
    
    return insertions, deletions

def write_updates(insertions, deletions, insertions_file="insertions.txt", deletions_file="deletions.txt"):
    """Write updates to files in the required format."""
    with open(insertions_file, 'w') as f:
        for u, v in insertions:
            f.write(f"I {u} {v}\n")
    
    with open(deletions_file, 'w') as f:
        for u, v in deletions:
            f.write(f"D {u} {v}\n")

if __name__ == "__main__":
    # partitions_dir = os.path.join("..", "MPI", "pre-parts")
    # insert_path = os.path.join("..", "MPI", "insertions.txt")
    # del_path = os.path.join("..", "MPI", "deletions.txt")
    partitions_dir = os.path.join("..", "pre-parts")
    insert_path = os.path.join("..", "dataset", "politicians", "insertions.txt")
    del_path = os.path.join("..", "dataset", "politicians", "deletions.txt")

    partition_edges = parse_partition_files(partitions_dir)
    
    # Adjust these numbers as needed
    num_insertions = NUM_INSERT
    num_deletions = NUM_DEL
    
    insertions, deletions = generate_updates(partition_edges, num_insertions, num_deletions)
    write_updates(insertions, deletions,insert_path,del_path)
    
    print(f"Generated {len(insertions)} insertions and {len(deletions)} deletions:")
    print("Insertions:", insertions[:5], "...")
    print("Deletions:", deletions[:5], "...")