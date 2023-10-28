import matplotlib.pyplot as plt

M_values = {'LAPACK': [], 'L3': [], 'L2': [], 'Basic': []}
time_values = {'LAPACK': [], 'L3': [], 'L2': [], 'Basic': []}

file_names = ['timing_results_LAPACK.txt', 'timing_results_L3.txt', 'timing_results_L2.txt', 'timing_results_Basic.txt']

for i, type_ in enumerate(['LAPACK', 'L3', 'L2', 'Basic']):
    with open(file_names[i], 'r') as file:
        for line in file:
            parts = line.strip().split()
            for j in range(len(parts)):
                if parts[j] == 'M:':
                    M_values[type_].append(int(parts[j + 1]))
                elif parts[j] == 'Time:':
                    time_values[type_].append(float(parts[j + 1]))


plt.figure(figsize=(8, 6))

for type_ in ['LAPACK', 'L3', 'L2', 'Basic']:
    plt.plot(M_values[type_], time_values[type_], marker='o', linestyle='-', label=type_)

plt.title('Matrix Size (M) vs. Time')
plt.xlabel('Matrix Size (M)')
plt.ylabel('Time (seconds)')
plt.grid(True)
plt.legend() 

plt.savefig('matrix_time_plot_comparison.png')
plt.show() 
