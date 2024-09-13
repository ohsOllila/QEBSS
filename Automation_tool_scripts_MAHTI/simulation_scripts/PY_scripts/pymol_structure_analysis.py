#!/usr/bin/python3

from pymol import cmd
import pymol
import glob
import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np


GRO='temp_md_1500ns.gro'
XTC='md_1500ns_smooth.xtc'
traj = md.load(XTC, top=GRO)

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

cmd.set("ray_opaque_background", 1)

rep_name = path.split('/')[-3]
forcefield = path.split('/')[-2]
selected = color_map.get(rep_name, 'black')

md = glob.glob(path + 'md*smooth*xtc')
temp = glob.glob(path + 'temp*gro')
if len(md) > 0 and len(temp) > 0:
	cmd.load(temp[0], rep_name + forcefield)
	cmd.load_traj(md[0], rep_name + forcefield, state=1, interval=2000)
	cmd.color(selected, rep_name + forcefield)
	cmd.set('all_states', 'on')
	cmd.ray(300, 300)
obj = cmd.get_object_list('all')
for i in obj[1:]:
	cmd.align(obj[0], i, object='aln', transform=0)
cmd.png(path + 'Ensemble_' + rep_name + '_' +  forcefield + '_aligned_fig.png')	
cmd.delete('all')
cmd.quit()

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

