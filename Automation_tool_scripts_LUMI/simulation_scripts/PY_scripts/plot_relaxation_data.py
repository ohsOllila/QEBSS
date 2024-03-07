#!/usr/bin/python3


import matplotlib.pyplot as plt
import os

ff_dir = os.path.basename(os.getcwd())
rep_dir = os.path.basename(os.path.dirname(os.getcwd()))
prot_dir = os.path.basename(os.path.dirname(os.path.dirname(os.getcwd())))

# Initialize lists to store data
T1_values = []
T2_values = []
NOE_values = []

# Read data from the text file
with open('relaxation_data.txt', 'r') as file:
	lines = file.readlines()
	for line in lines:
		# Split the line into T1, T2, and NOE values
		parts = line.split()
		T1_values.append(float(parts[1]))
		T2_values.append(float(parts[3]))
		NOE_values.append(float(parts[5]))

# Create a sequence (x-axis) for plotting
sequence = list(range(1, len(T1_values) + 1))

R1_values = [1 / value for value in T1_values]
R2_values = [1 / value for value in T2_values]
NOE_value_INV = [1 / value for value in NOE_values]


plt.plot(sequence, R1_values, marker='o', linestyle='-')
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
