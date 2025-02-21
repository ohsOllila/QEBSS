#!/usr/bin/python3


import matplotlib.pyplot as plt
import os
import shutil
import fileinput
import glob
import gc
import math
import csv
import subprocess
import re
import numpy as np
from collections import Counter
from matplotlib.image import imread
from itertools import combinations
import pandas as pd
import statistics
import random
import MDAnalysis as mda
from matplotlib.ticker import FuncFormatter
from pymol import cmd
import io
from PIL import Image
from matplotlib.backends.backend_pdf import PdfPages
import mdtraj as md
from matplotlib.patches import Rectangle
from math import log10, floor
from matplotlib.lines import Line2D
import pymol

three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}


FORCEFIELDS=["AMBER03WS", "AMBER99SB-DISP", "AMBER99SBWS", "DESAMBER", "CHARMM36M"]
REPLICAS=["replica_01", "replica_02", "replica_03", "replica_04", "replica_05"]

row_header =["Replica 1", "Replica 2", "Replica 3", "Replica 4", "Replica 5"]
color_list=['red', 'blue', 'green', 'purple', '#F08000']

color_map = {
    'AMBER03WS': 'red',
    'AMBER99SB-DISP': 'blue',
    'AMBER99SBWS': 'green',
    'CHARMM36M': 'purple',
    'DESAMBER': '#F08000'
}

color_map_sec = {
    'H': 'red',
    'E': 'blue',
    'C': 'gray',
}

SIM_DIR = os.getcwd() + '/'
BASE_DIR = os.path.dirname(os.path.dirname(SIM_DIR)) + '/'
exp_data = os.path.join(BASE_DIR, os.path.basename(os.path.normpath(SIM_DIR)) + '_exp_data.txt')
Unst_folder = os.path.join(SIM_DIR.replace(os.path.basename(os.path.normpath(SIM_DIR)), 'results'), os.path.basename(os.path.normpath(SIM_DIR))) + '/'
relax_folder = os.path.join(Unst_folder, 'rep_to_exp_data/')
avg_script = os.path.join(BASE_DIR, "simulation_scripts", "PY_scripts", "average_contact.py")
py_script = os.path.join(BASE_DIR, "simulation_scripts", "PY_scripts", "Old_Relaxations_for_Samuli.py")
pymol_analysis = os.path.join(BASE_DIR, "simulation_scripts", "PY_scripts", "pymol_structure_analysis.py")
best_cases_folder = os.path.join(relax_folder, "Accepted_cases/")
relax_data = sorted(glob.glob(SIM_DIR + 'replica*/*/relaxation_times.csv'))

os.makedirs(relax_folder, exist_ok=True)
os.makedirs(best_cases_folder, exist_ok=True)

with open(exp_data, 'r') as file:
	lines = file.readlines()
	first_row =lines[0].split()
	magn_field=first_row[5]
	res_nr=len(lines)-1

def convert_to_one_letter(three_letter_code):
	return three_to_one.get(three_letter_code, 'X')

seq_string=[]


with open("replica_01.pdb", 'r') as file:
	lines = file.readlines()
	for i in range(len(lines)):
		parts = lines[i].split()
		if len(parts) > 3 and parts[0] == "ATOM" and parts[2] == "CA": 
			one_letter_sequence = convert_to_one_letter(parts[3])
			seq_string.append(one_letter_sequence)

protein_sequence = ''.join(seq_string)



xticks = []
for i in range(9, res_nr, 25):
	xticks.append(i)
'''
xticks = []
for i in range(9, res_nr, 50):
	xticks.append(i)
'''

Names=['/'.join(i.split("/")[-3:-1]) for i in relax_data]


def create_amino_acid_pairs(sequence):
	sequence = sequence.upper()
	pairs = [(sequence[i], sequence[i + 1]) for i in range(len(sequence) - 1)]
	
	return pairs


def extract_columns(csv_paths, sim_case_nr, variable):
	result = {}
	for filename in csv_paths:
		with open(filename, 'r') as csvfile:
			sim_spec = '/'.join(filename.split("/")[-3:-1])
			reader = csv.reader(csvfile)
			header = next(reader)
			columns = [[] for _ in range(len(header))]
            
			for row in reader:
				for i, value in enumerate(row):
						columns[i].append(value)
			file_data = {}
			for col_name, col_values in zip(header, columns):
				file_data[col_name] = col_values
            
			result[sim_spec] = file_data
	
	spec_key=list(result.keys())[sim_case_nr]
	variable_list = result.get(spec_key, {}).get(variable, [])

	return variable_list

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

def extract_values_pandas_list(file_path, nr, column_number):
	df = pd.read_csv(file_path[nr])
	column = df.iloc[:, column_number]
	column_list = column.tolist()
	
	return column_list

def deviating_points(lst):
	temp_list = [(float(x) if x != "n" else float('nan'), idx) for idx, x in enumerate(lst)]
	valid_points = [(val, idx) for val, idx in temp_list if not (val != val)]
	valid_points.sort(key=lambda x: abs(x[0]), reverse=True)
	max_idx = [x[1] for x in valid_points[:4]]
    
	return max_idx






		

def axs_plot_avg(data_path, data_idx, column_nr, include_dev_points=False, include_header=True, axs=None, color_coord=True):
	x_header, x_list = extract_values_pandas(data_path, data_idx, 0, include_header)
	y_header, y_list = extract_values_pandas(data_path, data_idx, column_nr, include_header)

	plot_settings = {
		"name": Names[data_idx],
		"rep_name": Names[data_idx].split('/')[0],
		"forcefield": Names[data_idx].split('/')[1],
		"parameter": y_header.split('_')[0],
		"ax_label": None,
		"label": None,
		"color": "black" if ("exp" in y_header or not color_coord) else color_map.get(Names[data_idx].split('/')[1], "black")
	}

	if "sim" in y_header:
		plot_settings["color"] = "red"
		#plot_settings["label"] = f"{plot_settings['parameter']} simulation average"

	df = pd.DataFrame({x_header: x_list, y_header: y_list}).dropna()
	df.plot(x=x_header, y=y_header, ax=axs, label=plot_settings["label"], marker='o', linestyle="-", lw=5.0, markersize=7, color=plot_settings["color"])

	ticks = list(range(9, res_nr + 1, 20))
	if "hetNOE" in y_header:
		plot_settings["ax_label"] = f"{plot_settings['parameter']}"
		
		axs.set_xticks(ticks)
		axs.set_xticklabels([f"{i+1}" for i in ticks])
		#axs.set_xticklabels([f"{i+1+412}" for i in ticks])
		if "ICL2" in SIM_DIR:
			axs.set_xlabel('Residue', fontsize = 35)
		else:
			axs.set_xlabel("")

		axs.set_ylabel(plot_settings["ax_label"], fontsize=35)
		axs.set_yticks([-2, 0])
	
		axs.tick_params(axis='x', labelsize=30, direction='in', length=10, width=2)
		axs.tick_params(axis='y', labelsize=30, length=10, width=2)
		axs.set_ylim(-3, 0.8)
		axs.set_xlim(0, len(seq_string)+1)


	elif "R1" in y_header:
		plot_settings["ax_label"] = f"{plot_settings['parameter']}"
		axs.set_xticks(ticks)
		axs.tick_params(axis='x', direction='in', length=10, width=2)
		axs.set_xlabel("")
		axs.set_xticklabels([]) 
		axs.tick_params(axis='y', labelsize=30, length=10, width=2)
		axs.set_ylabel(plot_settings["ax_label"], fontsize=35, labelpad=30)
		axs.set_yticks([1, 2])
		axs.set_ylim(0, 3)
		axs.set_xlim(0, len(seq_string)+1)

	else:
		plot_settings["ax_label"] = f"{plot_settings['parameter']}"
		axs.set_xticks(ticks)
		axs.tick_params(axis='x', direction='in', length=10, width=2)
		axs.set_xlabel("")
		axs.set_xticklabels([]) 


		axs.tick_params(axis='y', labelsize=30, length=10, width=2)
		axs.set_ylabel(plot_settings["ax_label"], fontsize=35, labelpad=10)
		axs.set_yticks([5, 15])
		axs.set_ylim(0, 20)

		axs.set_xlim(0, len(seq_string)+1)

	axs.legend().remove()






avg_path=glob.glob(Unst_folder + "relaxation_times.csv")



fig, axs = plt.subplots(3, 1, figsize=(13, 6))
axs_plot_avg(avg_path, 0, 2, include_header=True, axs=axs[0])
axs_plot_avg(avg_path, 0, 5, include_header=True, axs=axs[1])
axs_plot_avg(avg_path, 0, 8, include_header=True, axs=axs[2])
axs_plot_avg(avg_path, 0, 1, include_header=True, axs=axs[0])
axs_plot_avg(avg_path, 0, 4, include_header=True, axs=axs[1])
axs_plot_avg(avg_path, 0, 7, include_header=True, axs=axs[2])
plt.tight_layout()
for ax in axs.flatten():  # Iterate through all subplots
    for spine in ax.spines.values():
        spine.set_linewidth(4)
plt.subplots_adjust(left=0.15, hspace=0)
plt.savefig(relax_folder + 'Accepted_cases/average_relaxation_compaired_plot.png', dpi=300)
plt.close('all')











