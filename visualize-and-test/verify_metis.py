import os

dataset_path = os.path.join("..", "dataset", "dataset1", "graph.txt")

edges = set()
with open(dataset_path, "r") as f:
    for line in f:
        u, v = map(int, line.strip().split())
        if (u, v) in edges or (v, u) in edges:
            print(f"Duplicate edge: {u} {v}")
        edges.add((u, v))
print(f"Total unique edges: {len(edges)}")