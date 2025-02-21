#!/usr/bin/python3

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob

# Load files
groFILE = glob.glob('temp_*.gro')[0]
trjFILE = glob.glob('md*noPBC*xtc')[0]
outfile_png = 'orientational_correlation.png'
outfile_csv = 'orientational_correlation.csv'
outfile_fit_png = 'flory_scaling_fit.png'

# Set up the MDAnalysis universe and select C-alpha atoms
u = mda.Universe(groFILE, trjFILE)
CAatoms = u.select_atoms("name CA")

# Number of C-alpha atoms and vectors
num_atoms = len(CAatoms.indices)
max_s = num_atoms - 1

# Initialize storage for average distances and counts
distances = np.zeros(max_s)
counts = np.zeros(max_s)

# Access only the 10,000th frame
frame_number = 10000

if frame_number < len(u.trajectory):
    u.trajectory[frame_number - 1]  # Frames are 0-indexed, so use frame_number - 1
    for s in range(1, max_s):
        for i in range(num_atoms - s):
            if i < i + s:  # Ensure each pair is calculated only once
                dist = np.linalg.norm(CAatoms[i].position - CAatoms[i + s].position)
                distances[s] += dist
                counts[s] += 1

# Average distances by dividing by the count
avg_distances = distances / counts

# Save average distances to a CSV file
np.savetxt(outfile_csv, np.column_stack((np.arange(1, max_s), avg_distances[1:])), delimiter=',', header="s, <R(s)>")

# Define a power-law function for fitting
def scaling_law(s, R0, nu):
    return R0 * s**nu

# Fit the scaling law to the short (s <= 10) and long (s > 10) distance data
s_values_short = np.arange(1, 11)
avg_distances_short = avg_distances[1:11]
params_short, _ = curve_fit(scaling_law, s_values_short, avg_distances_short)

s_values_long = np.arange(11, max_s)
avg_distances_long = avg_distances[11:]
params_long, _ = curve_fit(scaling_law, s_values_long, avg_distances_long)

# Plot the data and the fits
plt.plot(np.arange(1, max_s), avg_distances[1:], 'o', label='Average Distance <R(s)>')
plt.plot(s_values_short, scaling_law(s_values_short, *params_short), 'r-', label=f'Short fit: R0={params_short[0]:.2f}, nu={params_short[1]:.2f}')
plt.plot(s_values_long, scaling_law(s_values_long, *params_long), 'g-', label=f'Long fit: R0={params_long[0]:.2f}, nu={params_long[1]:.2f}')

plt.xlabel("Residue Separation $s$")
plt.ylabel(r"Average Distance $\langle R(s) \rangle$")
plt.title("Flory Scaling Law Fit for Average Intra-Protein Distance")
plt.legend()
plt.grid()
plt.tight_layout()

# Save the fit plot
plt.savefig(outfile_fit_png, dpi=600)
plt.close("all")

# Print the results
print("Short regime (s <= 10):")
print(f"  R0 = {params_short[0]:.2f}, nu = {params_short[1]:.2f}")

print("Long regime (s > 10):")
print(f"  R0 = {params_long[0]:.2f}, nu = {params_long[1]:.2f}")

