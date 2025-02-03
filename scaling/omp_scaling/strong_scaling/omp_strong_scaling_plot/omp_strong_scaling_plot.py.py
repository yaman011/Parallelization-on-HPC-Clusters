import matplotlib.pyplot as plt

# Strong scaling data with 64 CPUs included
execution_times = {
    1: 161.803089,
    2: 81.296086,
    4: 40.758637,
    8: 20.737949,
    16: 10.447832,
    27: 10.056666,
    32: 7.970795,
    48: 5.794110,
    64: 6.073456
}

# Calculate speedup
speedups = {cpus: execution_times[1] / time for cpus, time in execution_times.items()}

# Data for plotting
cpus = list(speedups.keys())
speedup_values = list(speedups.values())

# Plotting the strong scaling results
plt.figure(figsize=(10, 6))
plt.plot(cpus, speedup_values, marker='o', linestyle='-', color='blue')
plt.title('Strong Scaling: Speedup vs Number of CPUs', fontsize=16)
plt.xlabel('Number of CPUs', fontsize=14)
plt.ylabel('Speedup', fontsize=14)
plt.xticks(cpus, fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True)
plt.show()
