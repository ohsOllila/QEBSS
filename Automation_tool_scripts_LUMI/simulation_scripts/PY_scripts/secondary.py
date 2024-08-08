#!/usr/bin/python3

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np

traj = md.load('md_2000ns_smooth.xtc', top='temp_md_2000ns.gro')

three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

residue_names = [three_to_one[residue.name] for residue in traj.topology.residues]

dssp = md.compute_dssp(traj, simplified=True)

initial_secondary_structure = dssp[0]
final_secondary_structure = dssp[-1]

print(initial_secondary_structure)
print(final_secondary_structure)
color_map = {'C': 'grey', 'H': 'red', 'E': 'blue'}

initial_colors = [color_map[structure] for structure in initial_secondary_structure]
final_colors = [color_map[structure] for structure in final_secondary_structure]

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(18, 3))

for i, color in enumerate(initial_colors):
	ax1.bar(i, 1, color=color)

for i, color in enumerate(final_colors):
	ax2.bar(i, 1, color=color)

for ax in [ax1, ax2]:
	ax.set_xlim(0, len(dssp[0]))
	ax.set_ylim(0, 1)
	ax.set_yticks([])

xticks = [0]

for i in range(9, len(dssp[0]), 10):
	xticks.append(i)

print(xticks)
ax1.set_xticks(xticks)
ax1.set_xticklabels([f"{i+1}" for i in xticks], rotation=90)

ax2.set_xticks(range(len(dssp[0])))
ax2.set_xticklabels([f"{i+1} {residue_names[i]}" for i in range(len(dssp[0]))], rotation=90)
ax2.set_xlabel('Residue Index')

ax1.set_title('Initial Secondary Structure')
ax2.set_title('Final Secondary Structure')

# Add a legend
legend_labels = ['Coil (C)', 'Helix (H)', 'Extended (E)']
handles = [plt.Line2D([0], [0], color=color_map[code], lw=4) for code in color_map]
ax2.legend(handles, legend_labels, loc='upper right', bbox_to_anchor=(1.1, 1.0))

plt.tight_layout()
plt.show()

# Save the plot to a file
plt.savefig("secondary_with_residues.png")
