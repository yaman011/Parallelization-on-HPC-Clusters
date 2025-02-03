import matplotlib.pyplot as plt

# Number of processors and their corresponding nodes
processors = [1, 8, 16, 24, 32, 48, 64]
nodes = [1, 1, 2, 3, 4, 6, 8]  # Number of nodes used for each configuration

# Corresponding execution times (in seconds)
times = [176.21, 20.71, 11.12, 7.41, 5.57, 4.21, 3.11]

# Calculate speedup
speedup = [times[0] / time for time in times]

# Plotting the strong scaling graph
plt.figure(figsize=(12, 8))
plt.plot(processors, speedup, marker='o')
plt.title('Strong Scaling: Speedup vs Number of Processors (with Nodes)')
plt.xlabel('Number of Processors')
plt.ylabel('Speedup')
plt.xticks(processors)
plt.grid(True)

# Annotate the number of nodes for each point
for i, txt in enumerate(nodes):
    plt.annotate(f'{txt} Nodes', (processors[i], speedup[i]), textcoords="offset points", xytext=(0,10), ha='center')

plt.show()
