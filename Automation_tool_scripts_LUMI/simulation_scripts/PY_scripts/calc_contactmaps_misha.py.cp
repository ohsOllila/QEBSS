### calc_contactmaps.py

import os
import numpy as np
import pandas as pd
import glob
import MDAnalysis as mda
from MDAnalysis.analysis import distances

base_dir = os.path.dirname(os.getcwd()) + "/"
systems = glob.glob(base_dir + "Unst_alph*/model*/*/")

# Define cutoff distance (in Angstroms)
cutoffs = [5.0, 10.0, 15.0, 20.0]  # You can adjust this value as needed

# Loop through all directories
for cutoff in cutoffs:
	for path in systems:

		os.chdir(path)

		xtc_file = glob.glob(path + "md_*ns_noPBC.xtc")[0]
		gro_file = glob.glob(path + "temp*gro")[0]

		u = mda.Universe(gro_file, xtc_file)
		CAatoms = u.select_atoms("name CA")
		num_residues = len(CAatoms)

		within_cutoff_count = np.zeros((num_residues, num_residues))
		num_frames = 0

		for ts in u.trajectory:
			distances_array = distances.distance_array(CAatoms.positions, CAatoms.positions)
			within_cutoff_count += (distances_array <= cutoff)
			num_frames += 1

		probabilities = within_cutoff_count / num_frames

		df = pd.DataFrame(probabilities)

		df.columns = [f"{i+1}" for i in range(num_residues)]
		df.index = [f"{i+1}" for i in range(num_residues)]

		csv_filename = f"CA_prob_within_{cutoff}A.csv"
		df.to_csv(csv_filename)
		print(f"Saved {csv_filename} in {path}")

		os.chdir(base_dir)

print(f"All probabilities within {cutoff} Angstroms calculated and saved as CSV files.")
