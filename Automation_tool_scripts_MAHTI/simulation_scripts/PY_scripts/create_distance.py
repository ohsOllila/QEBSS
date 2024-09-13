import os
import numpy as np
import pandas as pd
import glob
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import matplotlib.pyplot as plt



# Function to calculate contact probabilities and save as CSV
def calculate_contact_probabilities(gro_file, xtc_file):
	u = mda.Universe(gro_file, xtc_file)
	CAatoms = u.select_atoms("name CA")
            
	num_residues = len(CAatoms)
	sum_distances = np.zeros((num_residues, num_residues))
	num_frames = 0
            
	# Iterate through all frames
	for ts in u.trajectory:
		distances_array = distances.distance_array(CAatoms.positions, CAatoms.positions)
		sum_distances += distances_array
		num_frames += 1
           
	# Calculate average distances
	avg_distances = sum_distances / num_frames
            
	# Create a DataFrame from the average distances
	df = pd.DataFrame(avg_distances)
            
	# Add column and index names
	df.columns = [f"{i+1}" for i in range(num_residues)]
	df.index = [f"{i+1}" for i in range(num_residues)]
            
	# Save the DataFrame as a CSV file
	csv_filename = "CA_avg_distances.csv"
	df.to_csv(csv_filename)
	print(f"Saved {csv_filename}")
       
	return df

# Function to create and save contact map plot
def create_contact_map_plot(probabilities_df, output_file):
#    vmin = np.floor(np.min(values.probabilities_df) / 10) * 10
#    vmax = np.ceil(np.max(values.probabilities_df) / 10) * 10
    vmin = 0
    vmax = 20

    color_map = plt.imshow(probabilities_df.values, vmin=vmin, vmax=vmax, cmap='jet', origin='lower')

    cbar = plt.colorbar(color_map)
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label('Distance (Å)', rotation=270, labelpad=25, fontsize=18)

    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.xlabel("Residue", fontsize=16)
    plt.ylabel("Residue", fontsize=16)

    plt.tight_layout()

    plt.savefig(output_file, dpi=600)
    plt.close()
    print(f"Distance map saved as {output_file}.")

# Get GRO and XTC files from the current directory
gro_file = glob.glob("*temp*.gro")[0]
xtc_file = glob.glob("*noPBC*.xtc")[0]

# Step 1: Calculate contact probabilities and save the CSV
probabilities_df = calculate_contact_probabilities(gro_file, xtc_file)

# Step 2: Create and save the contact map plot
create_contact_map_plot(probabilities_df, "Distance_map.png")

