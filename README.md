# Parallel SSSP Update Algorithm

This repository contains the implementation of a **Parallel Algorithm Template for Updating Single-Source Shortest Paths (SSSP)** in large-scale dynamic networks. The project explores sequential, MPI-based, and hybrid MPI + OpenMP solutions to efficiently update shortest paths after graph modifications.

## 📄 Paper Reference

[A Parallel Algorithm Template for Updating Single-Source Shortest Paths in Large-Scale Dynamic Networks](https://drive.google.com/file/d/1Cj7u6bLfbwfjSwtZhfDwv4V5wIuiGB9y/view)

## 🧠 Problem Overview

Given a large dynamic graph, the task is to update the shortest paths from a single source node after a batch of edge insertions or deletions. The algorithm must be scalable and optimized for performance.

## 🛠 Implementations

### 🔹 Sequential Version
- Standard Dijkstra-based update.
- Serves as a baseline for comparison.

### 🔹 MPI Version
- Distributed the graph using **METIS**.
- Parallel computation across nodes using **MPICH**.
- Communication handled via `MPI_Bcast`, `MPI_Gather`, `MPI_Allgather`, etc.

### 🔹 Hybrid MPI + OpenMP Version
- Each MPI process uses OpenMP for multi-threaded computation.
- Achieves the best speedup by leveraging both distributed and shared memory paradigms.

## 🚀 Performance Insights

- **MPI + OpenMP** version delivered the most significant speedups.
- **MPI-only** provided speedups but suffered from communication overhead at scale.
- **Sequential** time increased linearly with the number of updates.
- Communication time was excluded when measuring performance of MPI-based versions.

## 📊 Graphs and Visualizations

Performance plots comparing:
- Number of updates vs Time for each approach.
- Variations with MPI process count and OpenMP thread count.

(*Check `/Misc/` directory*)

## ⚙️ Requirements

- C++17
- MPICH (or compatible MPI implementation)
- OpenMP
- METIS
- GNU build tools

## 🧪 Build and Run

```bash
# Compile (adjust paths as needed)
mpic++ -std=c++17 -fopenmp -Iinclude main.cpp src/Graph.cpp src/SSSP.cpp -lmetis -o main

# Run (example with 4 processes)
mpiexec -np 4 ./main
