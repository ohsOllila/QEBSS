#!/usr/bin/python3

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np

# Load only specific frames (e.g., every 1000th frame between frames 1000 and 10000)
traj = md.load('md_2000ns_noPBC.xtc', top='temp_md_2000ns.gro', stride=2)

# Map residues to single-letter codes
three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}
residue_names = [three_to_one[residue.name] for residue in traj.topology.residues]

# Calculate DSSP secondary structure and coil probability
dssp = md.compute_dssp(traj, simplified=True)
coil_counts = np.sum(dssp == 'C', axis=0)
total_sampled_frames = dssp.shape[0]
coil_probabilities = coil_counts / total_sampled_frames


print(coil_counts)
# Print coil probabilities for each residue
for i, prob in enumerate(coil_probabilities):
    print(f"Residue {i+1} ({residue_names[i]}): Probability of coil = {prob:.2f}")

average_coil_probability = np.mean(coil_probabilities)
print(f"\nAverage probability of coil structure across all residues: {average_coil_probability:.2f}")

# Plot coil probabilities as a bar chart
plt.figure(figsize=(18, 4))
plt.bar(range(len(coil_probabilities)), coil_probabilities, color='grey')
plt.xlabel('Residue Index')
plt.ylabel('Coil Probability')
plt.title('Coil Probability for Each Residue')
plt.xticks(ticks=range(len(coil_probabilities)), labels=[f"{i+1} {residue_names[i]}" for i in range(len(coil_probabilities))], rotation=90)
plt.tight_layout()
plt.savefig("order_probability.png")

