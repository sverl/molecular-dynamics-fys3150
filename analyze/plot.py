import numpy as np
import matplotlib.pyplot as plt
import subprocess as subp
import glob
import sys

# Parameters
n_cells = 5   # FCC unit cells in each dimension
b = 5.26      # Lattice constant
steps = 5001  # Integration steps

# Variables
start_temps = float(sys.argv[1])   # Lowest inital temperature
end_temps = float(sys.argv[2])     # Highest initial temperature
num_temps = int(sys.argv[3])       # Number of temperatures

# Start a process for each simulation
temps_init = np.linspace(start_temps, end_temps, num_temps)
processes = []
for i, temp_init in enumerate(temps_init):
    print('\rStarting process %d/%d' % (i + 1, num_temps), end='')
    cmd = ['../mol_dyn', str(n_cells), str(temp_init), str(b), str(steps)]
    processes.append(subp.Popen(cmd, stdout=subp.PIPE))
print()

# Wait for each process to finish
for i, process in enumerate(processes):
    print('\rWaiting for process %d/%d' % (i + 1, num_temps), end='')
    process.wait()
print('\nDone!')

# Read the data files generated (also previous runs) and do some statistics
stat_files = glob.glob('statistics_*.txt')
temps_init_files = np.zeros(len(stat_files))
temps = np.zeros(len(stat_files))
diffs = np.zeros(len(stat_files))
for i, stat_file in enumerate(stat_files):
    stats = np.genfromtxt(stat_file).T
    steps, time, temp, E_kin, E_pot, E_tot, MSD = stats

    temps_init_files[i] = float(stat_file[11:-4]) * 119.735
    temps[i] = np.mean(temp)

    a, b = np.polyfit(time, MSD, 1)
    diffs[i] = a / 6

# Plot
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(temps, diffs, 'o', markerfacecolor='none')
plt.xlabel(r'$T/\mathrm{K}$')
plt.ylabel(r'$D/(\mathrm{\aa^2/ps})$')
plt.show()

plt.plot(temps_init_files, temps, 'o', markerfacecolor='none')
plt.xlabel(r'$T_i/\mathrm{K}$')
plt.ylabel(r'$T/\mathrm{K}$')
plt.show()

plt.plot(temps_init_files, temps / temps_init_files, 'o', markerfacecolor='none')
plt.xlabel(r'$T_i/\mathrm{K}$')
plt.ylabel(r'$T/Ti/\mathrm{K}$')
plt.show()
