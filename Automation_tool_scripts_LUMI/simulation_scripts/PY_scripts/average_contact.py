#!/usr/bin/python3



import matplotlib.pyplot as plt
import os
import shutil
import fileinput
import glob
import mdtraj as md
import math
import csv
import subprocess
import MDAnalysis as mda
import sys
import yaml
import re
import numpy as np
from collections import Counter
from matplotlib.image import imread
import pandas as pd
import statistics
from pymol import cmd
import pymol


SIM_DIR = os.getcwd()
BASE_DIR=os.path.dirname(SIM_DIR)
RESULTS=BASE_DIR + '/results/'



Best_rog_files=glob.glob(RESULTS + 'U*/Best_rog_landscape.txt')
Best_relax_files=glob.glob(RESULTS + 'U*/relaxation_times.txt')
Ctimes_files=glob.glob(RESULTS + 'U*/Ctimes_Coeffs.txt')

color_list=['red', 'blue', 'green', 'purple', 'orange']

rog_data=[[], [], []]

for i in Best_rog_files:
	values=[]
	counts=[]
	with open(i, 'r') as file:
		lines = file.readlines()
		rog_data[0].append((i.split("/")[-2], float(lines[0].split()[1])))
		for line in range(1, len(lines)):
			# Split the line into T1, T2, and NOE values
			parts = lines[line].split()
			values.append(float(parts[0]))
			counts.append(float(parts[1]))
	rog_data[1].append(values)
	rog_data[2].append(counts)
	
names=[item[0] for item in rog_data[0]]
rog_avg=[item[1] for item in rog_data[0]]



for i in range(len(Best_rog_files)):			
	plt.plot(rog_data[1][i], rog_data[2][i], color=color_list[i], label=names[i])
	plt.axvline(x=rog_avg[i], color=color_list[i], linestyle='--', label='Mean Rog value')
	plt.legend()

plt.xlabel('Count')
plt.ylabel('Radius of Gyration (nm)')
plt.title("Radius of Gyration landscape best cases")
plt.tight_layout()
plt.savefig(RESULTS + 'best_rog_landscape.png')
plt.close()

fig, axs = plt.subplots(3, 1, figsize=(11, 8))
for i, data in enumerate(sorted(Best_relax_files)):
	PROTEIN = data.split("/")[-2]
	print(PROTEIN)
	existing_res=[[], [], []]
	R1_values=[]
	R2_values=[]
	NOE_values=[]
	with open(data, 'r') as file:
		lines = file.readlines()
		for line in range(0, len(lines)):
			parts = lines[line].split()
			try:
				R1_values.append(1 / float(parts[1]))
				existing_res[0].append(line+1)
			except:
				pass
			try:
				R2_values.append(1 / float(parts[3]))
				existing_res[1].append(line+1)
			except:	
				pass
			try:
				NOE_values.append(float(parts[5]))
				existing_res[2].append(line+1)
			except:
				pass

	axs[0].plot(existing_res[0], R1_values, label=PROTEIN, marker='o', linestyle='-', lw=1.0, markersize=2, color=color_list[i])
	axs[0].set_xlabel('Residue number')
	axs[0].set_ylabel('R1_values (1/s)')
	axs[0].set_title('R1_data')
	axs[0].legend()
	axs[1].plot(existing_res[1], R2_values, label=PROTEIN, marker='o', linestyle='-', lw=1.0, markersize=2, color=color_list[i])
	axs[1].set_xlabel('Residue number')
	axs[1].set_ylabel('R2_values (1/s)')
	axs[1].set_title('R2_data')
	axs[1].legend()
	axs[2].plot(existing_res[2], NOE_values, label=PROTEIN, marker='o', linestyle='-', lw=1.0, markersize=2, color=color_list[i])
	axs[2].set_xlabel('Residue number')				
	axs[2].set_ylabel('NOE_values (1/s)')
	axs[2].set_title('NOE_data')
	axs[2].legend()
	
plt.tight_layout()
plt.savefig(RESULTS + 'Avg_relax_compair.png')
plt.close('all')


for i, data in enumerate(sorted(Best_relax_files)):
	PROTEIN = data.split("/")[-2]
	with open(data, 'r') as file:
		lines = file.readlines()
		y_values=[]
		for line in range(0, len(lines)):
			parts = lines[line].split()
			try:
				y_values.append(float(parts[7])*10**10)
			except:
				pass
		plt.plot(range(1, len(y_values)+1), y_values, label='R2 Data_' + PROTEIN, marker='o', linestyle='-', lw=1.0, markersize=2, color=color_list[i])
		plt.xlabel('Residue number')
		plt.ylabel('Effective correlation time (10^(-10) s)')
		plt.title('Effective correlation')
		plt.legend()
plt.tight_layout()
plt.savefig(RESULTS + 'Avg_tau_effective_area.png')
plt.close()

for i, data in enumerate(sorted(Ctimes_files)):
	PROTEIN = data.split("/")[-2]    
	x_vals, y_vals, weights = [], [], []
	with open(data, 'r') as file:
		lines = file.readlines()
		for line_idx in range(len(lines)):
			parts = lines[line_idx].split()
			if "Res" in str(parts[0]):
				x = float(parts[0].replace('Res_nr_', '')) 
				y_num = 1
				while line_idx + y_num < len(lines):
					next_parts = lines[line_idx + y_num].split()
					if "Res" not in next_parts[0]:
						try:
							y = float(next_parts[2].rstrip(','))
							weight = float(next_parts[3])
							if y < 89.125:
								x_vals.append(x)
								y_vals.append(y)
								weights.append(weight * 100)
							y_num += 1
						except (IndexError, ValueError):
							break
					else:
						break

	plt.scatter(x_vals, y_vals, s=weights, label=PROTEIN, zorder=5, color=color_list[i])
	plt.xlabel('Residue number')
	plt.ylabel('Timescales (ns)')
	plt.title('Timescale Scatter Plot')
	plt.legend()
plt.tight_layout()
plt.savefig(RESULTS + 'Avg_timescales.png')
plt.close()
