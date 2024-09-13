import os
import numpy as np
import pandas as pd
import glob
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import matplotlib.pyplot as plt


# Define cutoff distance (in Angstroms)
cutoff = 15.0  # Adjust as needed

# Function to calculate contact probabilities and save as CSV
def calculate_contact_probabilities(gro_file, xtc_file, cutoff):
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
    print(f"Saved {csv_filename} in current directory.")

    return df

# Function to create and save contact map plot
def create_contact_map_plot(probabilities_df, output_file):
    color_map = plt.imshow(probabilities_df.values, cmap='jet', vmin=0, vmax=1, origin='lower')

    cbar = plt.colorbar(color_map)
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label('Probability of contact (≤ 15.0 Å)', rotation=270, labelpad=25, fontsize=18)

    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.xlabel("Residue", fontsize=16)
    plt.ylabel("Residue", fontsize=16)

    plt.tight_layout()
    plt.savefig(output_file, dpi=600)
    plt.close()

    print(f"Contact map saved as {output_file}.")

# Get GRO and XTC files from the current directory
gro_file = glob.glob("*temp*.gro")[0]
xtc_file = glob.glob("*noPBC*.xtc")[0]

# Step 1: Calculate contact probabilities and save the CSV
probabilities_df = calculate_contact_probabilities(gro_file, xtc_file, cutoff)

# Step 2: Create and save the contact map plot
create_contact_map_plot(probabilities_df, "Contact_map.png")

