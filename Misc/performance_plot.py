import matplotlib.pyplot as plt

# Sequential: [(updates, time_in_sec), ...]
sequential_data = [
    (100, 0.5),
    (500, 1.3),
    (800, 3.1),
    (1000, 3.4),
]

# MPI: process_count -> list of (updates, time_in_sec)
mpi_data = {
    4: [(100, 0.4), (400, 0.4), (800, 0.42), (1000, 0.5)],
    8: [(100, 0.48), (400, 0.5), (800, 0.65), (1000, 0.8)],
}

# MPI + OpenMP: (mpi_procs, omp_threads) -> list of (updates, time_in_sec)
mpi_omp_data = {
    (4, 2): [(100, 0.2), (400, 0.3), (800, 0.35), (1000, 0.4)],
    (4, 4): [(100, 0.4), (400, 0.4), (800, 0.55), (1000, 0.75)],
    (8, 2): [(100, 0.54), (400, 0.7), (800, 1.2), (1000, 1.7)],
}

# ===============================
# Plotting
# ===============================

plt.figure(figsize=(10, 6))

# Plot Sequential
updates, times = zip(*sequential_data)
plt.plot(updates, times, label="Sequential", linewidth=2, marker='o')

# Plot MPI
for proc_count, data in mpi_data.items():
    u, t = zip(*data)
    plt.plot(u, t, label=f"MPI ({proc_count} procs)", linestyle='--', marker='s')

# Plot MPI + OpenMP
for (mpi_procs, omp_threads), data in mpi_omp_data.items():
    u, t = zip(*data)
    plt.plot(u, t, label=f"MPI ({mpi_procs}) + OMP ({omp_threads})", linestyle='-.', marker='^')


plt.xlabel("Number of Updates")
plt.ylabel("Time (ms)")
plt.title("SSSP Update Performance: Sequential vs MPI vs MPI+OpenMP \n Dataset: 20000 Nodes 70001 Edges")
plt.legend()
plt.grid(True)
plt.tight_layout()


plt.savefig("performance_comparison.png")
plt.show()
