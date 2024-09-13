### calc_distancemaps.py

import os
import numpy as np
import glob
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import distances


base_dir = os.path.dirname(os.getcwd()) + "/"
system = glob.glob(base_dir + "Unst_alpha*/model*/*/")

# Loop through all directories
for path in system:
	os.chdir(path)
	print(path)
	# Find the xtc and gro files
	xtc_file = glob.glob(path + "md_*ns_noPBC.xtc")[0]
	gro_file = glob.glob(path + "temp*gro")[0]
            
            
	# Create MDAnalysis universe
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
	print(f"Saved {csv_filename} in {path}")
            
	# Change back to the base directory
	os.chdir(base_dir)

print("All average distances calculated and saved as CSV files.")

