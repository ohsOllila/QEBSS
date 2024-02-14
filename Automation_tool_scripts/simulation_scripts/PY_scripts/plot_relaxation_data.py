#!/usr/bin/python3


import matplotlib.pyplot as plt
import os

ff_dir = os.path.basename(os.getcwd())
rep_dir = os.path.basename(os.path.dirname(os.getcwd()))
prot_dir = os.path.basename(os.path.dirname(os.path.dirname(os.getcwd())))

# Initialize lists to store data
R1_values_exp = []
R2_values_exp = []
NOE_values_exp = []

# Read data from the text file
with open('/scratch/project_2003809/cmcajsa/MD-stabilization/structures/forcefield_compare/simulation_scripts/MD_scripts/snear_exp_data.txt', 'r') as file:
	lines = file.readlines()
	for line in lines[1:]:
		# Split the line into T1, T2, and NOE values
		parts = line.split()
		print(parts)
		R1_values_exp.append(parts[1])
		R2_values_exp.append(parts[3])
		NOE_values_exp.append(parts[5])


# Create a sequence (x-axis) for plotting

sequence = list(range(1, len(R1_values_exp) + 1))



plt.plot(sequence, R1_values_exp, marker='o', linestyle='-')
plt.xlabel('residue number')
plt.ylabel('R_1 values')
plt.title(prot_dir + ' ' + rep_dir + ' ' + ff_dir)
plt.savefig('relaxation_R1_plot.png')
plt.close()

plt.plot(sequence, R2_values, marker='o', linestyle='-')
plt.xlabel('residue number')
plt.ylabel('R_2 values')
plt.title(prot_dir + ' ' + rep_dir + ' ' + ff_dir)
plt.savefig('relaxation_R2_plot.png')
plt.close()

plt.plot(sequence, NOE_value_INV, marker='o', linestyle='-')
plt.xlabel('residue number')
plt.ylabel('NOE values')
plt.title(prot_dir + ' ' + rep_dir + ' ' + ff_dir)
plt.savefig('relaxation_NOE_plot.png')
plt.close()
