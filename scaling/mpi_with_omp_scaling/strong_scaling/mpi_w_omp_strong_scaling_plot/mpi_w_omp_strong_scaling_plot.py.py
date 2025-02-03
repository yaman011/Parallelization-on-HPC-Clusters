import matplotlib.pyplot as plt

# Data for the configurations
nodes = [1, 1, 2, 2, 4, 4]
mpi_processes = [8, 4, 4, 8, 8, 4]
omp_threads = [2, 4, 4, 2, 2, 4]
execution_times = [203.3, 116.1, 124.8, 218.8, 241.5, 153.2]
config_labels = ["1n_8p2t", "1n_4p4t", "2n_4p4t", "2n_8p2t", "4n_8p2t", "4n_4p4t"]

# Calculate speedup and efficiency relative to the baseline (1 node, 8 MPI processes, 2 OMP threads)
baseline_time = execution_times[0]
speedup = [baseline_time / time for time in execution_times]
efficiency = [speedup[i] / nodes[i] for i in range(len(nodes))]

plt.figure(figsize=(14, 7))

# Execution Time Plot
plt.subplot(1, 3, 1)
for i, label in enumerate(config_labels):
    plt.plot(nodes[i], execution_times[i], 'o-', label=f"{label} ({mpi_processes[i]} MPI, {omp_threads[i]} OMP)")

plt.xlabel('Number of Nodes', fontsize=12)
plt.ylabel('Execution Time (seconds)', fontsize=12)
plt.title('Execution Time vs. Number of Nodes', fontsize=14)
plt.xticks([1, 2, 4], labels=['1 Node', '2 Nodes', '4 Nodes'])
plt.yticks(fontsize=10)
plt.legend(title="Configurations", fontsize=10, title_fontsize='11')
plt.grid(True, linestyle='--', alpha=0.6)

# Speedup Plot
plt.subplot(1, 3, 2)
for i, label in enumerate(config_labels):
    plt.plot(nodes[i], speedup[i], 'o-', label=f"{label} ({mpi_processes[i]} MPI, {omp_threads[i]} OMP)")

plt.xlabel('Number of Nodes', fontsize=12)
plt.ylabel('Speedup', fontsize=12)
plt.title('Speedup vs. Number of Nodes', fontsize=14)
plt.xticks([1, 2, 4], labels=['1 Node', '2 Nodes', '4 Nodes'])
plt.yticks(fontsize=10)
plt.legend(title="Configurations", fontsize=10, title_fontsize='11')
plt.grid(True, linestyle='--', alpha=0.6)

# Efficiency Plot
plt.subplot(1, 3, 3)
for i, label in enumerate(config_labels):
    plt.plot(nodes[i], efficiency[i], 'o-', label=f"{label} ({mpi_processes[i]} MPI, {omp_threads[i]} OMP)")

plt.xlabel('Number of Nodes', fontsize=12)
plt.ylabel('Efficiency', fontsize=12)
plt.title('Efficiency vs. Number of Nodes', fontsize=14)
plt.xticks([1, 2, 4], labels=['1 Node', '2 Nodes', '4 Nodes'])
plt.yticks(fontsize=10)
plt.legend(title="Configurations", fontsize=10, title_fontsize='11')
plt.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.show()
