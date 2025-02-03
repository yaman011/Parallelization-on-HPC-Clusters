import matplotlib.pyplot as plt

# Data from the weak scaling configurations
nodes = [1, 2, 2, 4, 4, 4]
grid_sizes = [
    "50x50x50", 
    "100x100x100", 
    "75x75x75", 
    "200x200x200", 
    "50x50x50", 
    "120x120x120"
]
execution_times = [11.27, 46.99, 18.74, 494.13, 220.63, 280.80]

# Calculate efficiency as (baseline time / current time) * number of nodes
baseline_time = execution_times[0]
efficiency = [(baseline_time / execution_times[i]) * nodes[i] for i in range(len(nodes))]

plt.figure(figsize=(16, 8))

# Execution Time Plot
plt.subplot(1, 2, 1)
for i, label in enumerate(grid_sizes):
    plt.plot(nodes[i], execution_times[i], 'o-', label=f"{label} ({nodes[i]} Nodes)")

plt.xlabel('Number of Nodes', fontsize=14)
plt.ylabel('Execution Time (seconds)', fontsize=14)
plt.title('Execution Time vs. Number of Nodes (Weak Scaling)', fontsize=16)
plt.xticks([1, 2, 4], labels=['1 Node', '2 Nodes', '4 Nodes'], fontsize=12)
plt.yticks(fontsize=12)
plt.legend(title="Grid Sizes", fontsize=12, title_fontsize='13')
plt.grid(True, linestyle='--', alpha=0.7)

# Efficiency Plot
plt.subplot(1, 2, 2)
for i, label in enumerate(grid_sizes):
    plt.plot(nodes[i], efficiency[i], 'o-', label=f"{label} ({nodes[i]} Nodes)")

plt.xlabel('Number of Nodes', fontsize=14)
plt.ylabel('Efficiency', fontsize=14)
plt.title('Efficiency vs. Number of Nodes (Weak Scaling)', fontsize=16)
plt.xticks([1, 2, 4], labels=['1 Node', '2 Nodes', '4 Nodes'], fontsize=12)
plt.yticks(fontsize=12)
plt.legend(title="Grid Sizes", fontsize=12, title_fontsize='13')
plt.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
plt.show()
