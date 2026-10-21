#!/usr/bin/python3

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob
import os

# Define simulations
selected = ["replica_02/AMBER99SBWS", "replica_03/AMBER99SBWS", "replica_05/AMBER99SBWS"]

# Load files
path = os.getcwd()
groFILES = [glob.glob(f"{path}/{sim}/temp_*.gro")[0] for sim in selected]
trjFILES = [glob.glob(f"{path}/{sim}/md*noPBC*xtc")[0] for sim in selected]

# Outputs
outfile_fit_png = "flory_scaling_fit_combined.png"
outfile_csv = "orientational_correlation_combined.csv"

# Storage for combined results
all_avg_distances = []

for groFILE, trjFILE in zip(groFILES, trjFILES):
    print(f"Processing: {groFILE}, {trjFILE}")
    
    # MDAnalysis universe
    u = mda.Universe(groFILE, trjFILE)
    CAatoms = u.select_atoms("name CA")
    num_atoms = len(CAatoms.indices)
    max_s = num_atoms - 1
    
    distances = np.zeros(max_s)
    counts = np.zeros(max_s)
    
    # Specify frame interval
    frame_interval = 1000
    
    # Loop through the trajectory and calculate distances
    for frame in range(0, len(u.trajectory), frame_interval):
        u.trajectory[frame]
        for s in range(0, max_s):
            for i in range(num_atoms - s):
                dist = np.linalg.norm(CAatoms[i].position - CAatoms[i + s].position)
                distances[s] += dist
                counts[s] += 1

    # Calculate average distances
    avg_distances = distances / counts
    all_avg_distances.append(avg_distances)

# Combine results: average over all simulations
all_avg_distances = np.array(all_avg_distances)
avg_distances_combined = np.nanmean(all_avg_distances, axis=0)  # Average across simulations

# Save combined results to CSV
np.savetxt(outfile_csv, np.column_stack((np.arange(1, max_s), avg_distances_combined[1:])), delimiter=',', header="s, <R(s)>")

# Fit the scaling law to short (s <= 10) and long (s > 90) data
def scaling_law(s, R0, nu):
    return R0 * s**nu

s_values_short = np.arange(1, 11)
avg_distances_short = avg_distances_combined[1:11]
params_short, _ = curve_fit(scaling_law, s_values_short, avg_distances_short)

s_values_long = np.arange(11, max_s)
avg_distances_long = avg_distances_combined[11:]
params_long, _ = curve_fit(scaling_law, s_values_long, avg_distances_long)

R0_short, nu_short = params_short
R0_long, nu_long = params_long

# Plot the combined results
plt.figure(figsize=(8, 6))

# Plot individual results
#for avg_distances, sim in zip(all_avg_distances, selected):
#    plt.plot(np.arange(1, max_s), avg_distances[1:], '--', label=f"{sim}")

# Plot combined results
plt.plot(np.arange(1, max_s), avg_distances_combined[1:], 'o', label="Combined Data")

# Add fits
plt.plot(s_values_short, scaling_law(s_values_short, *params_short), 'r-', label=f"Fit: $s \\leq 10$")
plt.plot(s_values_long, scaling_law(s_values_long, *params_long), 'g-', label=f"Fit: $s > 10$")

# Annotate fit equations
x_center = (plt.xlim()[0] + plt.xlim()[1]) / 2
y_center = (plt.ylim()[0] + plt.ylim()[1]) / 2

plt.text(
    x_center - (x_center * 0.95),
    y_center - (y_center * 0.90),
    f'$R(s) = {R0_short:.2f} \cdot s^{{{nu_short:.2f}}}$\n($s \\leq 10$)',
    color='red',
    fontsize=10,
    ha='center'
)

plt.text(
    x_center + (x_center * 0.25),
    y_center * 0.7,
    f'$R(s) = {R0_long:.2f} \cdot s^{{{nu_long:.2f}}}$\n($s > 10$)',
    color='green',
    fontsize=10,
    ha='center'
)

# Finalize plot
plt.xscale('log')
plt.yscale('log')
plt.title("SNAP25")
plt.xlabel("Residue Separation $s$")
plt.ylabel(r"Average Distance $R_s$ (nm)")
plt.tight_layout()

# Save and show
plt.savefig(outfile_fit_png, dpi=600)
plt.show()

# Print results
print("Combined Short regime (s <= 10):")
print(f"  R0 = {R0_short:.2f}, v = {nu_short:.2f}")
print("Combined Long regime (s > 90):")
print(f"  R0 = {R0_long:.2f}, v = {nu_long:.2f}")
