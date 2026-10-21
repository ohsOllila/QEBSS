#!/usr/bin/python3

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob
import os

# Load files
path=os.getcwd()
groFILE = glob.glob('temp_*.gro')[0]
trjFILE = glob.glob('md*noPBC*xtc')[0]
name=path.split("/")[-3]
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

# Specify the frame interval (every 1000th frame)
frame_interval = 1000

# Loop through the trajectory and process every 1000th frame
for frame in range(0, len(u.trajectory), frame_interval):
    u.trajectory[frame]  # This moves to the specified frame number
    for s in range(0, max_s):
        for i in range(num_atoms - s):
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

s_values_long = np.arange(90, max_s)
avg_distances_long = avg_distances[90:]
params_long, _ = curve_fit(scaling_law, s_values_long, avg_distances_long)

# Plot the data and the fits
plt.plot(np.arange(1, max_s), avg_distances[1:], 'o')
plt.plot(s_values_short, scaling_law(s_values_short, *params_short), 'r-')
plt.plot(s_values_long, scaling_law(s_values_long, *params_long), 'g-')


R0_short = params_short[0]
nu_short = params_short[1]

# After fitting the long-range (s > 10) data
R0_long = params_long[0]
nu_long = params_long[1]

print(R0_long)

x_center = (plt.xlim()[0] + plt.xlim()[1]) / 2
y_center = (plt.ylim()[0] + plt.ylim()[1]) / 2

# Add the short-range equation slightly above the center
plt.text(
    x_center - (x_center * 0.9), 
    y_center - (y_center * 0.7),  # Offset slightly upward
    f'$R(s) = {R0_short:.2f} \cdot s^{{{nu_short:.2f}}}$\n($s \leq 10$)', 
    color='red', 
    fontsize=12, 
    ha='center'  # Horizontal alignment
)

# Add the long-range equation slightly below the center
plt.text(
    x_center + (x_center * 0.15), 
    y_center - (y_center * 0.1),  # Offset slightly downward
    f'$R(s) = s^{{{nu_long:.2f}}}$\n($s > 90$)',
    color='green', 
    fontsize=12, 
    ha='center'
)


# Finalize plot
plt.xscale('log')
plt.yscale('log')
plt.title(name)
plt.xlabel("Residue Separation $s$")
plt.ylabel(r"Average Distance $R_s$ (nm)")
plt.tight_layout()

# Save the fit plot
plt.savefig(outfile_fit_png, dpi=600)
plt.close("all")

# Print the results
print("Short regime (s <= 10):")
print(f"  R0 = {params_short[0]:.2f}, v = {params_short[1]:.2f}")

print("Long regime (s > 10):")
print(f"  R0 = {params_long[0]:.2f}, v = {params_long[1]:.2f}")

# Print the total number of frames in the trajectory
num_frames = len(u.trajectory)
print(f"Total number of frames in the trajectory: {num_frames}")
