import matplotlib.pyplot as plt

# Number of threads
threads = [1, 2, 4, 8, 16, 27, 32, 64]

# Corresponding execution times (in seconds)
times = [163.050, 230.565, 326.589, 458.286, 679.120, 860.719, 931.500, 1368.960]

# Corresponding grid spacing (delta) for each configuration
grid_spacings = [0.01, 0.00707, 0.005, 0.00354, 0.0025, 0.00193, 0.00177, 0.00125]

# Plotting the weak scaling graph
plt.figure(figsize=(12, 8))
plt.plot(threads, times, marker='o')
plt.title('Weak Scaling: Execution Time vs Number of Threads (OMP)')
plt.xlabel('Number of Threads')
plt.ylabel('Execution Time (seconds)')
plt.xticks(threads)
plt.grid(True)

# Annotate the grid spacing for each point
for i, txt in enumerate(grid_spacings):
    plt.annotate(f'Grid Spacing: {txt}', (threads[i], times[i]), textcoords="offset points", xytext=(0,10), ha='center')

plt.show()
