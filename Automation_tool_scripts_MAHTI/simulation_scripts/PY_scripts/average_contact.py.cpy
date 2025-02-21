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
from matplotlib.table import Table
import pandas as pd
import statistics
from pymol import cmd
import pymol


SIM_DIR = os.getcwd()
BASE_DIR=os.path.dirname(SIM_DIR)
RESULTS=BASE_DIR + '/results/'

PROTEINS = ["SNARE", "ChiZ", "KRS", "aSyn", "ICL2"]

Best_rog_files=glob.glob(RESULTS + '*/Best_rog_landscape.txt')

Best_relax_files=glob.glob(RESULTS + '*/relaxation_times.csv')
Ctimes_files=glob.glob(RESULTS + '*/Ctimes_Coeffs.csv')

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

for i in sorted(Best_rog_files):
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

print(rog_data)

for i in range(len(Best_rog_files)):			
	plt.plot(rog_data[1][i], rog_data[2][i], color=color_list[i], label=names[i])
	plt.axvline(x=rog_avg[i], color=color_list[i], linestyle='--', label='Mean Rog value')
	plt.legend()

plt.xlabel('Count')
plt.ylabel('Radius of Gyration (nm)')
plt.title("Radius of Gyration landscape best cases")
plt.tight_layout()
plt.savefig(RESULTS + 'best_rog_landscape.png', dpi=300)
plt.close()



fig, ax = plt.subplots()
for i, data in enumerate(Best_relax_files):
	tau, tau_values=extract_values_pandas(Best_relax_files, i, 10, include_header=True)
	res, res_values=extract_values_pandas(Best_relax_files, i, 0, include_header=True)

	df = pd.DataFrame({tau: tau_values, res: res_values})
	df.dropna(inplace=True)

	PROTEIN = data.split("/")[-2]

	df.plot(x=res, y=tau, label=PROTEIN, marker='o', linestyle='-', lw=1.0, markersize=2, color=color_list[i], ax=ax)

plt.xlabel('Residue number')
plt.ylabel('Effective correlation time (ns)')
plt.title('Effective correlation')
plt.legend()
plt.tight_layout()


plt.savefig(RESULTS + 'Avg_tau_effective_area.png', dpi=2000)
plt.close()


def plot_images(input, output, title=True):
	col=len(input)
	fig, axs = plt.subplots(col, 1, figsize=(5, 11.69))
	for i, data in enumerate(input):
		data=input[i]
		PROTEIN = data.split("/")[-4]
		img = imread(data)
		if col == 1:
			plt.imshow(img)
			if title:
				plt.title(PROTEIN, fontweight='bold')
			plt.axis('off')
		else:
			axs[i].imshow(img)
			if title:
				axs[i].set_title(PROTEIN, fontweight='bold')
			axs[i].axis('off')
	plt.tight_layout()
	plt.savefig(RESULTS + output + '.pdf', dpi=350)
	plt.close()

def plot_images_row(input, output, title=True):
	col=len(input)
	fig, axs = plt.subplots(1, col, figsize=(12, 6))
	for i, data in enumerate(input):
		data=input[i]
		PROTEIN = data.split("/")[-4]
		img = imread(data)
		if col == 1:
			plt.imshow(img)
			if title:
				plt.title(PROTEIN, fontweight='bold')
			plt.axis('off')
		else:
			axs[i].imshow(img)
			if title:
				axs[i].set_title(PROTEIN, fontweight='bold')
			axs[i].axis('off')
	plt.tight_layout()
	plt.savefig(RESULTS + output + '.png', dpi=350)
	plt.close()

def plot_images_extrem(input, output, title=True):
	col=len(input)
	fig, axs = plt.subplots(col, 1, figsize=(6, 12))
	for i, data in enumerate(input):
		data=input[i]
		PROTEIN = data.split("/")[-4]
		img = imread(data)
		if col == 1:
			plt.imshow(img)
			if title:
				plt.title(PROTEIN, fontweight='bold')
			plt.axis('off')
		else:
			axs[i].imshow(img)
			if title:
				axs[i].set_title(PROTEIN, fontweight='bold', fontsize=16, loc='left', rotation=90, pad=20, verticalalignment='top', y=+0.57, x=-0.03)
			axs[i].axis('off')

	labels_top = ["R1", "R2", "hetNOE"]
	positions_top = [0.21, 0.52, 0.84]	
	for label, pos in zip(labels_top, positions_top):
		fig.text(pos, 0.95, label, fontsize=16, fontweight='bold', ha='center')

	fig.text(0.15, 0.06, "red", ha='center', fontsize=12, fontweight='bold', color='red')
	fig.text(0.19, 0.06, " = worst", ha='left', fontsize=12, fontweight='bold', color='black')

	fig.text(0.45, 0.06, "yellow", ha='center', fontsize=12, fontweight='bold', color='gold')
	fig.text(0.50, 0.06, " = middle", ha='left', fontsize=12, fontweight='bold', color='black')

	fig.text(0.77, 0.06, "green", ha='center', fontsize=12, fontweight='bold', color='green')
	fig.text(0.81, 0.06, " = best", ha='left', fontsize=12, fontweight='bold', color='black')

	fig.text(0.52, 0.1, "RESIDUES", ha='center', fontsize=16, fontweight='bold', color='black')

	plt.tight_layout(rect=[0.03, 0.03, 1, 0.95])
#	plt.tight_layout()
	plt.subplots_adjust(bottom=0.15, hspace = 0.07)
	plt.savefig(RESULTS + output + '.pdf', dpi=350)
	plt.close()

def plot_images_multiple_sets(input_sets, output):
    num_sets = len(input_sets)
    num_images = max(len(images) for images in input_sets)

    fig, axs = plt.subplots(num_sets, num_images, figsize=(8.27, 11.69))
    
    for row, image_set in enumerate(input_sets):
        for col, data in enumerate(image_set):
            PROTEIN = data.split("/")[-4]
            img = imread(data)

            ax = axs[row, col] if num_sets > 1 else axs[col]
            ax.imshow(img)
            if "align" in data:
                ax.set_title(PROTEIN, fontweight='bold')
            ax.axis('off')

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.05)
    plt.savefig(RESULTS + output + '.pdf', dpi=300)
    plt.close()

avg_structure = [glob.glob(RESULTS + i + '/rep_to_exp_data/Accepted_cases/Ensemble_*.png')[0] for i in PROTEINS]
avg_contact_map = [glob.glob(RESULTS + i + '/rep_to_exp_data/Accepted_cases/Avg_contact.png')[0] for i in PROTEINS]
avg_distance_map = [glob.glob(RESULTS + i + '/rep_to_exp_data/Accepted_cases/Avg_dist.png')[0] for i in PROTEINS]
avg_correlation_map = [glob.glob(RESULTS + i + '/rep_to_exp_data/Accepted_cases/Avg_corr.png')[0] for i in PROTEINS]
#avg_timescales_times = [glob.glob(RESULTS + i + '/rep_to_exp_data/Accepted_cases/Timescale_plot_avg.png')[0] for i in PROTEINS]

plot_images_multiple_sets([avg_structure, avg_contact_map, avg_distance_map, avg_correlation_map], "All_Images_Combined")

print([glob.glob(RESULTS + i + '/rep_to_exp_data/Accepted_cases/Ensemble_*.png')[0] for i in PROTEINS])

relaxation_png = [glob.glob(RESULTS + i + '/rep_to_exp_data/Accepted_cases/average_relaxation_compaired_plot.png')[0] for i in PROTEINS]
plot_images(relaxation_png, "Avg_relaxation_times")

Extreme_png = [glob.glob(RESULTS + i + '/rep_to_exp_data/Accepted_cases/Extreme_cases_plot.png')[0] for i in PROTEINS]
plot_images_extrem(Extreme_png, "Extremes")


#Dens_png = [glob.glob(RESULTS + i + '/rep_to_exp_data/density_landscape_plot.png')[0] for i in PROTEINS]
#plot_images(Dens_png, "Density_all", title=False)

Dens_png = [glob.glob(RESULTS + i + '/rep_to_exp_data/Average_rog_plot.png')[0] for i in PROTEINS]
plot_images_row(Dens_png, "Density_all", title=False)

