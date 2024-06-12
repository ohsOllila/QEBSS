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
Best_relax_files=glob.glob(RESULTS + 'U*/relaxation_times.csv')
Ctimes_files=glob.glob(RESULTS + 'U*/Ctimes_Coeffs.csv')

color_list=['red', 'blue', 'green', 'purple', 'orange']

def extract_values_pandas(file_path, nr, column_number, include_header=False):
        df = pd.read_csv(file_path[nr])
        if include_header:
                header = df.columns[column_number]
                column = df.iloc[:, column_number]
                column = pd.to_numeric(column, errors='coerce')
                return header, column
        else:
                column = df.iloc[:, column_number]
                column = pd.to_numeric(column, errors='coerce')
                return column

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



for i, data in enumerate(sorted(Best_relax_files)):
	tau, tau_values=extract_values_pandas(Best_relax_files, i, 10, include_header=True)
	res, res_values=extract_values_pandas(Best_relax_files, i, 0, include_header=True)

	df = pd.DataFrame({tau: tau_values, res: res_values})
	df.dropna(inplace=True)

	PROTEIN = data.split("/")[-2]

	
	df.plot(x=res, y=tau, label=PROTEIN, marker='o', linestyle='-', lw=1.0, markersize=2, color=color_list[i])
	plt.xlabel('Residue number')
	plt.ylabel('Effective correlation time (ns)')
	plt.title('Effective correlation')
	plt.legend()
plt.tight_layout()
plt.savefig(RESULTS + 'Avg_tau_effective_area.png')
plt.close()

def plot_images(input, output):
	col=len(input)
	if 'relaxation_compaired' in input[0]:
		fig, axs = plt.subplots(col, 1, figsize=(8.27, 11.69))
	else:
		fig, axs = plt.subplots(1, col, figsize=(12, 4))
	for i, data in enumerate(input):
		data=input[i]
		PROTEIN = data.split("/")[-4]
		img = imread(data)
		if col == 1:
			plt.imshow(img)
			plt.title(PROTEIN)
			plt.axis('off')
		else:
			axs[i].imshow(img)
			axs[i].set_title(PROTEIN)
			axs[i].axis('off')
	plt.tight_layout()
	plt.savefig(RESULTS + output + '.png')
	plt.close()


contact_png = sorted(glob.glob(RESULTS + 'U*/rep_to_exp_data/Accepted_cases/Avg_corr.png'))
plot_images(contact_png, "Avg_correlation_map")

contact_png = sorted(glob.glob(RESULTS + 'U*/rep_to_exp_data/Accepted_cases/Avg_contact.png'))
plot_images(contact_png, "Avg_contact_map")

contact_png = sorted(glob.glob(RESULTS + 'U*/rep_to_exp_data/Accepted_cases/Ensemble_*.png'))
plot_images(contact_png, "Avg_structure")


relaxation_png = sorted(glob.glob(RESULTS + 'U*/rep_to_exp_data/Accepted_cases/average_relaxation_compaired_plot.png'))
plot_images(relaxation_png, "Avg_relaxation_times")

contact_png = sorted(glob.glob(RESULTS + 'U*/rep_to_exp_data/Accepted_cases/Timescale_plot_best.png'))
plot_images(contact_png, "Avg_timescales_times")
