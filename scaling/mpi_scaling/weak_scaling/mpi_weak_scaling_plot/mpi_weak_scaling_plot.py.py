import matplotlib.pyplot as plt

# Number of processors
processors = [1, 8, 16, 32, 48, 64]

# Corresponding execution times (in seconds)
times = [159.24, 161.51, 160.71, 207.01, 220.89, 344.84]

# Corresponding grid spacing (delta) for each configuration
grid_spacings = [0.01, 0.005, 0.004, 0.003, 0.0025, 0.002]

# Plotting the weak scaling graph
plt.figure(figsize=(12, 8))
plt.plot(processors, times, marker='o')
plt.title('Weak Scaling: Execution Time vs Number of Processors')
plt.xlabel('Number of Processors')
plt.ylabel('Execution Time (seconds)')
plt.xticks(processors)
plt.grid(True)

# Annotate the grid spacing for each point
for i, txt in enumerate(grid_spacings):
    plt.annotate(f'Grid Spacing: {txt}', (processors[i], times[i]), textcoords="offset points", xytext=(0,10), ha='center')

plt.show()
