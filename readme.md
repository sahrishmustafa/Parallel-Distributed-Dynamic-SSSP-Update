# Dynamic SSSP Updates in Large Graphs using MPI, OpenMP, and METIS

This project implements a **parallel dynamic Single-Source Shortest Paths (SSSP) algorithm** using:

- MPI (inter-node parallelism)
- OpenMP / OpenCL (intra-node parallelism)
- METIS (graph partitioning)
- Real public graph datasets

---

##  Project Goals

- Implement and analyze a **parallel dynamic SSSP algorithm** as described in the paper:  
  *“A Parallel Algorithm Template for Updating Single-Source Shortest Paths in Large-Scale Dynamic Networks”*
- Support dynamic **edge insertions/deletions**
- Evaluate **scalability** and **performance** on real datasets
- Compare:
  - Sequential version
  - MPI-only
  - MPI + OpenMP hybrid

---

## Repository Structure

| Folder | Contents |
|--------|----------|
| `sequential/` | Pure sequential version with Dijkstra + dynamic updates |
| `mpi/`        | Distributed implementation using MPI and METIS |
| `mpi_openmp/` | Hybrid MPI + OpenMP for node-local parallelism |
| `datasets/`   | Real graph datasets in edge list format |
| `results/`    | Timing results, plots, scalability data |
| `visualize-and-test/`    | Utilities for testing algorithm against scratch implementation and visualization |

---

## Datasets

We use real public graphs:
- [SNAP: BHJ, Orkut](https://snap.stanford.edu/data/)
- Synthetic RMAT graphs (Graph500)
- Custom edge-weighted test graphs (soon)

[Use the following link to download them](https://drive.google.com/drive/folders/1xWA5EUnzGv_pB71-b9JpQfxfZaAHuFDm?usp=sharing)

---
