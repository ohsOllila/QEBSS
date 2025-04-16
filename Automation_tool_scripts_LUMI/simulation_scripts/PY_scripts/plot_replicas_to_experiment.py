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
from pymol import cmd
import io
from PIL import Image
from matplotlib.backends.backend_pdf import PdfPages
import mdtraj as md
from matplotlib.patches import Rectangle
from math import log10, floor
from matplotlib.lines import Line2D
import pymol
from sklearn.preprocessing import MinMaxScaler

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

xticks = []
for i in range(9, res_nr, 50):
	xticks.append(i)


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

def remove_largest(lst):
	var_type=type(lst[0])
	temp_list = [float(x) for x in lst if x != "n"]
	temp_list.sort(key=abs, reverse=True)
	max_val = temp_list[:4]

	for val in max_val:
		lst.remove(var_type(val))

	return lst


def find_largest_value(exp, sim):
	ymax_list=[]
	ymin_list=[]
	for i in range(len(Names)):
		exp_data_clean = [float(i) for i in extract_columns(relax_data, i, exp) if i != 'n']
		sim_data_clean = [float(j) for j in extract_columns(relax_data, i, sim) if j != 'n']

		ymax_list.append(max(exp_data_clean))
		ymax_list.append(max(sim_data_clean))

		ymin_list.append(min(exp_data_clean))
		ymin_list.append(min(sim_data_clean))
	ymax=max(ymax_list)
	ymin=min(ymin_list)

	return ymin, ymax



def calculate_rmsd(idx, diff):
	values_diff = extract_columns(relax_data, idx, diff)
	cleaned_list = [float(i) for i in values_diff if i != 'n']
	
	return math.sqrt(sum(float(x) ** 2 for x in remove_largest(cleaned_list)) / len(remove_largest(cleaned_list)))


def calculate_rmsre(idx, exp, diff):
	exp_data=extract_columns(relax_data, idx, exp)
	diff_data=extract_columns(relax_data, idx, diff)

	index_dev=deviating_points(diff_data)

	exp_data = [exp_data[i] for i in range(len(exp_data)) if i not in index_dev] 
	diff_data = [diff_data[j] for j in range(len(diff_data)) if j not in index_dev]

	relative_errors = [(float(diff)) / float(a) for a, diff in zip(exp_data, diff_data) if a != "n" and diff != "n"]
	squared_relative_errors = [e ** 2 for e in relative_errors]
	mean_squared_relative_error = sum(squared_relative_errors) / len(squared_relative_errors)
	rmsre = math.sqrt(mean_squared_relative_error)

	return rmsre


def ranking_value(diff):
	RMSD_lists=[]
	for i in range(len(Names)):
		RMSD_lists.append(calculate_rmsd(i, diff))
		#RMSD_lists.append(calculate_rmsre(i, exp, diff))
	rank=min(float(i) for i in RMSD_lists)

	return rank

def plot_all(file_path, output_file_name, execution):
	fig, axs = plt.subplots(5, 5, figsize=(15, 15))
	for ff_idx, ff in enumerate(FORCEFIELDS):
		for idx in range(len(file_path)):
			rep_name=file_path[idx].split("/")[-3]
			rep_idx=int(rep_name.split("/")[0][-1])-1
			if ff in file_path[idx]:
				execution(file_path, idx, axs[rep_idx, ff_idx])
	plt.tight_layout()
	plt.savefig(relax_folder + output_file_name)
	plt.close()


def avg_csv(input, use_header=False, use_index=False):
	header_option = 'infer' if use_header else None
	index_option = 0 if use_header else None

	dataframes = [pd.read_csv(file, header=header_option, index_col=index_option) for file in input]
	stacked_array = np.array([df.values for df in dataframes])
	mean_array = np.mean(stacked_array, axis=0)
	matrix = pd.DataFrame(mean_array)

	return matrix

def round_sig(x, sig=2):
	scale = sig - int(floor(log10(abs(x)))) - 1
	rounded_value = round(x, scale)
    
	return rounded_value

def read_rog_values(data_path):
	rog_values = []
	with open(data_path, 'r') as file:
		for line in file.readlines()[27:]:
			try:
				value = float(line.split()[1])
				rog_values.append(value)
			except IndexError:
				pass
	rounded_values = [round(value, 1) for value in rog_values]
	rounded_counts = {value: rounded_values.count(value) for value in set(rounded_values)}
	sorted_items = sorted(rounded_counts.items())
	values = [value[0] for value in sorted_items]
	counts = [value[1] for value in sorted_items]


	return rog_values, counts, values


def axs_plot(data_path, data_idx, column_nr, include_dev_points=False, include_header=True, axs=None, color_coord=True):
	x_header, x_list = extract_values_pandas(data_path, data_idx, 0, include_header)
	y_header, y_list = extract_values_pandas(data_path, data_idx, column_nr, include_header)

	name = Names[data_idx]
	rep_name = Names[data_idx].split('/')[0]
	forcefield = Names[data_idx].split('/')[1]

	parameter = y_header.split('_')[0]

	if color_coord == True:
		color = color_map.get(forcefield, "black")
		if "exp" in y_header:
			color = "black"
	elif "results" in data_path[data_idx]:
		color = "red"
	else:
		color = color_list[int(rep_name[-1]) - 1]

	df = pd.DataFrame({x_header: x_list, y_header: y_list}).dropna()
	df.plot(y=y_header, ax=axs, marker='o', linestyle="-", lw=5, markersize=6, color=color)

	axs.legend().remove()
	if "sim" in y_header:
			rank_label=Best_case_sum[data_idx][2 if parameter == "R1" else 3 if parameter == "R2" else 4 if parameter == "hetNOE" else None]
			legend_item = Line2D([0], [0], lw=0, color='w', label=f'Rank value = $\\mathbf{{{str(rank_label)}\\%}}$')
			if "R2" in y_header:
				axs.legend(handles=[legend_item], fontsize=25, loc='upper right', frameon=False)
			else:
				axs.legend(handles=[legend_item], fontsize=25, loc='lower right', frameon=False)
			#axs.set_title(f'Rank value = {str(rank_label)}%', fontweight='bold', fontsize=25)
			if include_dev_points==True:
				max_idx= deviating_points(extract_columns(data_path, data_idx, parameter + "_diff"))
				y_values= [extract_columns(data_path, data_idx, y_header)[i] for i in max_idx]
				x_values= [extract_columns(data_path, data_idx, x_header)[i] for i in max_idx]
				#axs.scatter(x_values, y_values, color='black', s=50)
	
	if "hetNOE" in y_header:
		y_label=parameter
	elif "Tau" in y_header:
		y_label="Eff.Corr (ns)"
	else:
		y_label=parameter + " (1/s)"
	
	y_min, y_max = find_largest_value(y_header, y_header.replace("exp", "sim"))
	if "hetNOE" in y_header:
		yticks = np.linspace(round(y_min, 1), round(y_max, 1), 4)
	else:
		yticks = np.linspace(0, round(y_max, 1), 4)
	

	print(yticks)

	yticks = np.round(yticks, 1)
	axs.set_xticks([])
	axs.set_yticks(yticks)
	axs.tick_params(axis='y', labelsize=25)
	axs.set_xlim(-2, res_nr+1)

	if "hetNOE" not in y_header:
		axs.set_ylim(-y_max*0.1, y_max+y_max*0.1)
		
	

	return y_label
		

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
	df.plot(x=x_header, y=y_header, ax=axs, label=plot_settings["label"], marker='o', linestyle="-", lw=3.0, markersize=5, color=plot_settings["color"])

	if "hetNOE" in y_header:
		plot_settings["ax_label"] = f"{plot_settings['parameter']}"
		ticks = list(range(10, res_nr + 1, 10))
		axs.set_xticks(ticks)
		#axs.set_xticklabels([f"{i+1}" for i in ticks])
		axs.set_xticklabels([f"{i+1}" for i in ticks])
		if "ICL2" in SIM_DIR:
			axs.set_xlabel('Residue', fontsize = 30)
		else:
			axs.set_xlabel("")

		axs.set_ylabel(plot_settings["ax_label"], fontsize=25)
		axs.tick_params(axis='x', labelsize=20)
		axs.tick_params(axis='y', labelsize=20)

	else:
		plot_settings["ax_label"] = f"{plot_settings['parameter']} (1/s)"
		axs.set_xlabel("")
		axs.set_xticks([])
		axs.tick_params(axis='y', labelsize=25)
		axs.set_ylabel(plot_settings["ax_label"], fontsize=25)

	ymin, ymax = axs.get_ylim()
	axs.set_yticks([round(i, 1) for i in np.linspace(ymin, ymax, 4)])
	axs.set_xlim(0-0.5, res_nr+0.5)


def axs_plot_extremes(data_path, data_idx, column_nr, color, include_header=True, axs=None):
	x_header, x_list = extract_values_pandas(data_path, data_idx, 0, include_header)
	y_header, y_list = extract_values_pandas(data_path, data_idx, column_nr, include_header)

	name = Names[data_idx],
	rep_name = Names[data_idx].split('/')[0]
	forcefield = Names[data_idx].split('/')[1]



	df = pd.DataFrame({x_header: x_list, y_header: y_list}).dropna()
	df.plot(x=x_header, y=y_header, ax=axs, marker='o', linestyle="-", lw=3, markersize=3.5, color=color)

	axs.tick_params(which="both", labelsize=16)
	axs.get_legend().remove()

	axs.set_xlabel("")
	axs.set_ylabel("")

def calculate_coil_probability(path, output, stride=20):
    # Mapping of residue names to single-letter codes
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    # List to accumulate the coil counts and total frame count
    total_coil_counts = None
    total_frames = 0
    residue_names = None


    for file in path:
        traj_file = glob.glob(file + "/md*noPBC*xtc")[0]
        print(traj_file)
        gro_file = glob.glob(file + "/temp*gro")[0]

        # Load the trajectory file with the corresponding topology
        traj = md.load(traj_file, top=gro_file, stride=stride)

        # Map residues to single-letter codes
        if residue_names is None:
            residue_names = [three_to_one[residue.name] for residue in traj.topology.residues]

        # Compute DSSP (secondary structure) for the trajectory
        dssp = md.compute_dssp(traj, simplified=True)

        # Count coil (C) residues across frames
        coil_counts = np.sum(dssp == 'C', axis=0)

        # Update the total coil counts and frame count
        if total_coil_counts is None:
            total_coil_counts = coil_counts
        else:
            total_coil_counts += coil_counts  # Accumulate coil counts across frames

        total_frames += dssp.shape[0]  # Add the number of frames in this trajectory

    # Calculate the average coil probability across all frames
    average_coil_probabilities = total_coil_counts / total_frames

    # Plot the combined coil probabilities as a bar chart
    plt.figure(figsize=(14, 4), dpi=300)
    plt.bar(range(len(average_coil_probabilities)), average_coil_probabilities, color='grey', width=0.5)
    plt.xlabel('Residue Index')
    plt.ylabel('Average Coil Probability')
    #plt.xticks(ticks=range(len(average_coil_probabilities)), labels=[f"{i+1} {residue_names[i]}" for i in range(len(average_coil_probabilities))], rotation=90, fontsize=5)
    plt.gca().set_xticks(range(len(average_coil_probabilities)))
    plt.gca().set_xticklabels([f"{i + 1} {residue_names[i]}" for i in range(len(average_coil_probabilities))], rotation=90, fontsize=6)
    plt.tight_layout()

    plt.savefig(output)

def merge_images_vertically(pymol_path, second_path, output_path):
    # Open the images
    img1 = Image.open(pymol_path)
    img2 = Image.open(second_path)

    # Ensure both images have the same width
    width = min(img1.width, img2.width)  # Use the smaller width
    img1 = img1.resize((width, int(img1.height * (width / img1.width))))  # Resize maintaining aspect ratio
    img2 = img2.resize((width, int(img2.height * (width / img2.width))))  # Resize maintaining aspect ratio

    # Calculate the combined height
    total_height = img1.height + img2.height

    # Create a new blank image with the combined dimensions
    merged_image = Image.new('RGB', (width, total_height))

    # Paste the images on top of each other
    merged_image.paste(img1, (0, 0))
    merged_image.paste(img2, (0, img1.height))

    # Save the resulting image
    merged_image.save(output_path)

def csv_to_pdf_with_conditional_formatting(csv_file, value_min_threshold=None, value_max_threshold=None):
	df = pd.read_csv(csv_file)
	pdf_file=csv_file.replace("csv", "pdf")

	fig, ax = plt.subplots(figsize=(12, 8))
	ax.axis('off')

	ax.axis('tight')
	ax.axis('off')
	ax.set_title(SIM_DIR.split("/")[-2], fontsize=16, pad=20)

	table = ax.table(cellText=df.values,
	colLabels=df.columns,
	cellLoc='center',
	loc='center',
	colColours=["lightgray"] * df.shape[1])

	for (i, j), val in np.ndenumerate(df.values):
		try:
			if value_min_threshold <= int(val) <= value_max_threshold:
				table[(i+1, j)].set_facecolor('lightgray')  # +1 to skip the header row
		except:
			pass
	table.auto_set_font_size(False)
	table.set_fontsize(10)
	table.scale(1.2, 1.2) 
	
	with PdfPages(pdf_file) as pdf:
		pdf.savefig(fig, bbox_inches='tight')



RMSD_R1=[]
RMSD_R2=[]
RMSD_hetNOE=[]
for i in range(len(Names)):	
	RMSD_R1.append(calculate_rmsd(i, "R1_diff"))
	RMSD_R2.append(calculate_rmsd(i, "R2_diff"))
	RMSD_hetNOE.append(calculate_rmsd(i, "hetNOE_diff"))




Best_cases=[]
Best_case_sum = []

ranking_file=relax_folder + 'ranking_table.csv'

with open(ranking_file, 'w', newline="") as csvfile:
    ranking_results = []
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["Force field", "Replica", "R1 RMSD", "R1 (%)", "R2 RMSD", "R2 (%)", "hetNOE RMSD", "hetNOE (%)", "Sum (%)"])

    for i in range(len(Names)):
        case = Names[i]
        rep_name = case.split("/")[0]
        ff_name = case.split("/")[1]

        RMSD_R1 = calculate_rmsd(i, "R1_diff")
        RMSD_R2 = calculate_rmsd(i, "R2_diff")
        RMSD_hetNOE = calculate_rmsd(i, "hetNOE_diff")

        R1_rank_value = (float(RMSD_R1) / ranking_value("R1_diff")) * 100
        R2_rank_value = (float(RMSD_R2) / ranking_value("R2_diff")) * 100
        hetNOE_rank_value = (float(RMSD_hetNOE) / ranking_value("hetNOE_diff")) * 100
        
        Ranking_sum = R1_rank_value + R2_rank_value + hetNOE_rank_value
        
        ranking_results.append([
            ff_name, rep_name, 
            round_sig(RMSD_R1), int(round_sig(R1_rank_value)), 
            round_sig(RMSD_R2), int(round_sig(R2_rank_value)), 
            round_sig(RMSD_hetNOE), int(round_sig(hetNOE_rank_value)), 
            int(round_sig(Ranking_sum))
        ])
        
        Best_case_sum.append([case, int(round_sig(Ranking_sum)), int(round_sig(R1_rank_value)), int(round_sig(R2_rank_value)), int(round_sig(hetNOE_rank_value))])
        
        if round_sig(R1_rank_value) / 100 <= 1.5 and round_sig(R2_rank_value) / 100 <= 1.5 and round_sig(hetNOE_rank_value) / 100 <= 1.5:
            Best_cases.append(case)
            print(case)
        else:
            print(case + " do not meet criteria")

    ranking_results.sort(key=lambda x: x[-1])

    for result in ranking_results:
        csvwriter.writerow(result)



print(Best_case_sum)
if len(Best_cases)==0:
	lowest_case = min(Best_case_sum, key=lambda x: x[1])
	Best_cases.append(lowest_case[0])

print(Best_cases)

worst_idx=Names.index(max(Best_case_sum, key=lambda x: x[1])[0])
middle_idx = Names.index(sorted(Best_case_sum, key=lambda x: x[1])[len(Best_case_sum) // 2][0])

if len(Best_cases) > 0:
	best_idx = Names.index(min(Best_case_sum, key=lambda x: x[1])[0])
else:
	best_idx=Names.index(Best_cases[0])
print(Names[best_idx])


best_worst_middle_cases_text = best_cases_folder + "best_worst_middle_cases_list.csv"

with open(best_worst_middle_cases_text, 'w', newline='') as file:
	writer = csv.writer(file)
	writer.writerow(['Color', 'Directory'])

	writer.writerow(['#029f4cff', SIM_DIR + Names[best_idx]])
	writer.writerow(['#0024f0ff', SIM_DIR + Names[middle_idx]])
	writer.writerow(['#ff0202ff', SIM_DIR + Names[worst_idx]])

csv_to_pdf_with_conditional_formatting(ranking_file, value_min_threshold = 100, value_max_threshold=150)

os.makedirs(Unst_folder + "correlation_functions/", exist_ok=True)

i = 1
while i <= res_nr :
	corr_lists=[]
	for j in Best_cases:
		try:
			corr_list=[]
			file_name = SIM_DIR + j + "/correlation_functions/NHrotaCF_"+ str(i) + ".xvg"
			with open(file_name, 'r') as file:
				lines = file.readlines()
				for line in range(0, len(lines)-1):
					parts =lines[line].split()
					corr_list.append(float(parts[1]))
			corr_lists.append(corr_list)
		except:
			pass
	try:
		average_list=[sum(x) / len(x) for x in zip(*corr_lists)]
		times = [i * 10.000 for i in range(0, len(corr_lists[0]))]
		new_file=Unst_folder + 'correlation_functions/NHrotaCF_' + str(i) + '.xvg'
		with open(new_file, 'w') as f:
			for x, y in zip(times, average_list):
				f.write(f"{x:.3f}\t{y:.5f}\n")
	
	except:
		pass
	corr_lists=[]			
	i += 1

shutil.copy(py_script, Unst_folder)

with fileinput.FileInput(Unst_folder + "Old_Relaxations_for_Samuli.py", inplace=True) as file:
	for line in file:
		line = line.replace('magn_field=magn_field', 'magn_field=' + str(magn_field))
		print(line, end="")

os.chdir(Unst_folder)
subprocess.run(["python3", Unst_folder + "Old_Relaxations_for_Samuli.py"])


avg_path=glob.glob(Unst_folder + "relaxation_times.csv")


os.chdir(SIM_DIR)

for item in FORCEFIELDS:
	fig, axs = plt.subplots(1, 3, figsize=(15, 6))
	for i, name in enumerate(Names):
		if name.split("/")[1]==item:
			axs_plot(relax_data, i, 2, include_header=True, axs=axs[0], color_coord=False)
			axs_plot(relax_data, i, 5, include_header=True, axs=axs[1], color_coord=False)
			axs_plot(relax_data, i, 8, include_header=True, axs=axs[2], color_coord=False)
	else:
		axs_plot(relax_data, i, 1, include_header=True, axs=axs[0], color_coord=False)
		axs_plot(relax_data, i, 4, include_header=True, axs=axs[1], color_coord=False)
		axs_plot(relax_data, i, 7, include_header=True, axs=axs[2], color_coord=False)
	plt.tight_layout()
	plt.savefig(relax_folder + item + '_plot.png')
	plt.close('all')

fig, axs = plt.subplots(1, 3, figsize=(15, 6))
for j in Best_cases:
	i=Names.index(j)
	axs_plot(relax_data, i, 2, include_header=True, axs=axs[0])
	axs_plot(relax_data, i, 5, include_header=True, axs=axs[1])
	axs_plot(relax_data, i, 8, include_header=True, axs=axs[2])
axs_plot(relax_data, i, 1, include_header=True, axs=axs[0])
axs_plot(relax_data, i, 4, include_header=True, axs=axs[1])
axs_plot(relax_data, i, 7, include_header=True, axs=axs[2])
plt.tight_layout()
plt.savefig(relax_folder + 'Accepted_cases/Best_cases_plot')
plt.close('all')

fig, axs = plt.subplots(1, 3, figsize=(8.27, 2.338))


axs_plot_extremes(relax_data, best_idx, 2, "#029f4cff", axs=axs[0])
axs_plot_extremes(relax_data, best_idx, 5, "#029f4cff", axs=axs[1])
axs_plot_extremes(relax_data, best_idx, 8, "#029f4cff", axs=axs[2])

axs_plot_extremes(relax_data, middle_idx, 2, "#0024f0ff", axs=axs[0])
axs_plot_extremes(relax_data, middle_idx, 5, "#0024f0ff", axs=axs[1])
axs_plot_extremes(relax_data, middle_idx, 8, "#0024f0ff", axs=axs[2])

axs_plot_extremes(relax_data, worst_idx, 2, "#ff0202ff", axs=axs[0])
axs_plot_extremes(relax_data, worst_idx, 5, "#ff0202ff", axs=axs[1])
axs_plot_extremes(relax_data, worst_idx, 8, "#ff0202ff", axs=axs[2])


axs_plot_extremes(relax_data, i, 1, "black", axs=axs[0])
axs_plot_extremes(relax_data, i, 4, "black", axs=axs[1])
axs_plot_extremes(relax_data, i, 7, "black", axs=axs[2])
plt.tight_layout()
plt.savefig(relax_folder + 'Accepted_cases/Extreme_cases_plot.png', dpi=300)
plt.savefig(relax_folder + 'Accepted_cases/Extreme_cases_plot.pdf', dpi=300)
plt.close('all')


#fig, axs = plt.subplots(3, 1, figsize=(12, 6))

#axs_plot_extremes(relax_data, best_idx, 2, "Simulation", "blue", axs=axs[0])
#axs_plot_extremes(relax_data, best_idx, 5, "Simulation ", "blue", axs=axs[1])
#axs_plot_extremes(relax_data, best_idx, 8, "Simulation", "blue", axs=axs[2])

#axs_plot_extremes(relax_data, best_idx, 1, "NMR data", "black", axs=axs[0])
#axs_plot_extremes(relax_data, best_idx, 4, "NMR data", "black", axs=axs[1])
#axs_plot_extremes(relax_data, best_idx, 7, "NMR data", "black", axs=axs[2])
#plt.tight_layout()
#plt.savefig(relax_folder + 'Accepted_cases/Article_relaxation_plot', dpi=300)
#plt.close('all')


fig, axs = plt.subplots(3, 1, figsize=(13, 6))
axs_plot_avg(avg_path, 0, 2, include_header=True, axs=axs[0])
axs_plot_avg(avg_path, 0, 5, include_header=True, axs=axs[1])
axs_plot_avg(avg_path, 0, 8, include_header=True, axs=axs[2])
axs_plot_avg(avg_path, 0, 1, include_header=True, axs=axs[0])
axs_plot_avg(avg_path, 0, 4, include_header=True, axs=axs[1])
axs_plot_avg(avg_path, 0, 7, include_header=True, axs=axs[2])
plt.tight_layout()
plt.subplots_adjust(left=0.10, hspace=0.14)
plt.savefig(relax_folder + 'Accepted_cases/average_relaxation_compaired_plot.png', dpi=300)
plt.close('all')


#axs_plot(avg_path, 0, 10, include_header=True, axs=None)
#plt.savefig(relax_folder + 'Accepted_cases/Tau_effective_area_avg.png')
#plt.close('all')


def relaxation_combined(sim, exp, output, execution):
	fig, axs = plt.subplots(5, 5, figsize=(25, 25))
	fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.05, wspace=0.07, hspace=0.07)
	for idx, item in enumerate(Names):
		rep_name = item.split("/")[0]
		forcefield = item.split("/")[1]
		i = int(rep_name[-1]) - 1
		j=FORCEFIELDS.index(forcefield)
		y_label = execution(relax_data, idx, sim, include_dev_points=True, include_header=True, axs=axs[i, j])
		if exp is not None:
			execution(relax_data, idx, exp, include_dev_points=True, include_header=True, axs=axs[i, j])
		if rep_name + "/" + forcefield in Best_cases:
			autoAxis = axs[i, j].axis()
			print(autoAxis)
			rec = Rectangle((autoAxis[0]-(autoAxis[1] - autoAxis[0]) * 0.05, autoAxis[2]- (autoAxis[3] - autoAxis[2]) * 0.05), 
					(autoAxis[1] - autoAxis[0])*1.1, 
					(autoAxis[3] - autoAxis[2])*1.1, 
					fill=False, linewidth=8, edgecolor='black')
			axs[i, j].add_patch(rec)
			rec.set_clip_on(False)  # Ensure the rectangle is not clipped by axes limits
		#if item == "replica_01/AMBER03WS":
		#	axs[i, j].imshow(np.zeros((10, 10)), cmap='gray')
		#	axs[i, j].set_xticks([])
		#	axs[i, j].set_yticks([])
	for ax1, replica in zip(axs[:,0], row_header):
		ax1.set_ylabel(f"$\\bf{{{replica}}}$\n\n{y_label}", fontsize=30, rotation=90)
	for ax2, forcefield in zip(axs[0], FORCEFIELDS):
		color = color_map.get(forcefield.split("/")[-1], "black")
		ax2.set_title(forcefield, weight='bold', fontsize=30, pad=22, color = color)
	for ax3 in axs[-1]:
		ax3.set_xlabel('Residue', fontsize=30)
		ax3.tick_params(axis='x', labelsize=27)
		#axs.set_xticklabels([f"{i+1}" for i in xticks])
		ax3.set_xticks(xticks)
		ax3.set_xticklabels([f"{i+1}" for i in xticks])
	for row in axs:
		for a in row[1:]:  # Skip the first column
			a.set_yticks([])
	plt.savefig(relax_folder + output +'.png')
	plt.close()

relaxation_combined(1, 2, 'R1_relaxation_combined_plot', axs_plot)

relaxation_combined(4, 5, 'R2_relaxation_combined_plot', axs_plot)
relaxation_combined(7, 8, 'hetNOE_relaxation_combined_plot', axs_plot)

relaxation_combined(10, None, "Tau_effective_area_all", axs_plot)

def distance_map(input, output):
    all_data = []
    for rep in REPLICAS:
        for ff in FORCEFIELDS:
            file = os.path.join(SIM_DIR, rep, ff, "CA_avg_distances.csv")
            if os.path.exists(file):
                df = pd.read_csv(file, index_col=0)
                all_data.append(df.values)
    all_data = np.array(all_data)  # Convert to a NumPy array
  
    if "dist" in input:
        vmin = np.floor(np.min(all_data) / 10) * 10
        vmax = 20 
    else:
        vmin=0
        vmax=1
    fig, axs = plt.subplots(5, 5, figsize=(15, 15))  # 5x5 grid, adjust as needed
    fig.subplots_adjust(left=0.1, right=0.85, top=0.85, bottom=0.1, wspace=0.07, hspace=0.07)


    for i, rep in enumerate(REPLICAS):
        for j, ff in enumerate(FORCEFIELDS):
            if i >= 5 or j >= 5:  # Avoid index errors with fixed grid size
                continue
            
            file = os.path.join(SIM_DIR, rep, ff, input)
            if os.path.exists(file):
                df = pd.read_csv(file, index_col=0)
                
                color_map = axs[i, j].imshow(df.values, cmap='jet', vmin=vmin, vmax=vmax, origin='lower')
                axs[i, j].tick_params(axis='both', which='major', labelsize=16)
                axs[i, j].set_xticks([])
                axs[i, j].set_yticks(xticks)
                axs[i, j].set_yticklabels([f"{i+1}" for i in xticks])
                if rep + "/" + ff in Best_cases:
                    autoAxis = axs[i, j].axis()
                    print(autoAxis)
                    rec = Rectangle((autoAxis[0]-(autoAxis[1] - autoAxis[0]) * 0.04, autoAxis[2]- (autoAxis[3] - autoAxis[2]) * 0.04), 
                        (autoAxis[1] - autoAxis[0])*1.08, 
                        (autoAxis[3] - autoAxis[2])*1.08, 
                        fill=False, linewidth=5, edgecolor='black')
                    axs[i, j].add_patch(rec)
                    rec.set_clip_on(False)

    cbar_ax = fig.add_axes([0.87, 0.1, 0.02, 0.75])  # Position for the colorbar
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(plt.cm.ScalarMappable(cmap='jet', norm=norm), cax=cbar_ax)
    cbar.ax.tick_params(labelsize=22)
    if "dist" in input:
        cbar.set_label('Distance (Å)', rotation=270, labelpad=35, fontsize=25)
    else:
        cbar.set_label('Probability of contact (≤ 15.0 Å)', rotation=270, labelpad=35, fontsize=25)
    for ax1, replica in zip(axs[:,0], row_header):
        ax1.set_ylabel(f"$\\bf{{{replica}}}$\n\n Residue", fontsize=16, rotation=90)
    for ax2, forcefield in zip(axs[0], FORCEFIELDS):
        ax2.set_title(forcefield, weight='bold', fontsize=16, pad=19, color = "black")
    for ax3 in axs[-1]:
        ax3.set_xlabel('Residue', fontsize=16)
        ax3.tick_params(axis='x', labelsize=16)
        #axs.set_xticklabels([f"{i+1}" for i in xticks])
        ax3.set_xticks(xticks)
        ax3.set_xticklabels([f"{i+1}" for i in xticks])
    for row in axs:
        for a in row[1:]:  # Skip the first column
            a.set_yticks([])


    output_file = os.path.join(relax_folder, output)
    plt.savefig(output_file, dpi=600)
    plt.close()

distance_map("CA_avg_distances.csv", "Distance_map_combined.png")
distance_map("CA_prob_within_15.0A.csv", "Contact_map_combined.png")

def correlation_map(input, output):
    fig, axs = plt.subplots(5, 5, figsize=(15, 15))  # 5x5 grid, adjust as needed
    fig.subplots_adjust(left=0.1, right=0.85, top=0.85, bottom=0.1, wspace=0.07, hspace=0.07)

    vmin=-1
    vmax=1
    for i, rep in enumerate(REPLICAS):
        for j, ff in enumerate(FORCEFIELDS):
            if i >= 5 or j >= 5:  # Avoid index errors with fixed grid size
                continue

            
            file = os.path.join(SIM_DIR, rep, ff, input)
            if os.path.exists(file):
                df = pd.read_csv(file, index_col=0)
                
                color_map = axs[i, j].imshow(df.values, cmap='seismic', vmin=-1, vmax=1, origin='lower')
                axs[i, j].tick_params(axis='both', which='major', labelsize=16)
                axs[i, j].set_xticks([])
                axs[i, j].set_yticks(xticks)
                axs[i, j].set_yticklabels([f"{i+1}" for i in xticks])
                if rep + "/" + ff in Best_cases:
                    autoAxis = axs[i, j].axis()
                    print(autoAxis)
                    rec = Rectangle((autoAxis[0]-(autoAxis[1] - autoAxis[0]) * 0.04, autoAxis[2]- (autoAxis[3] - autoAxis[2]) * 0.04), 
                        (autoAxis[1] - autoAxis[0])*1.08, 
                        (autoAxis[3] - autoAxis[2])*1.08, 
                        fill=False, linewidth=5, edgecolor='black')
                    axs[i, j].add_patch(rec)
                    rec.set_clip_on(False)

    cbar_ax = fig.add_axes([0.87, 0.1, 0.02, 0.75])  # Position for the colorbar
    cbar = plt.colorbar(plt.cm.ScalarMappable(cmap='seismic'), cax=cbar_ax)
    cbar.ax.tick_params(labelsize=22)
    cbar.set_label('Correlation value', rotation=270, labelpad=35, fontsize=25)
    for ax1, replica in zip(axs[:,0], row_header):
        ax1.set_ylabel(f"$\\bf{{{replica}}}$\n\n Residue", fontsize=16, rotation=90)
    for ax2, forcefield in zip(axs[0], FORCEFIELDS):
        ax2.set_title(forcefield, weight='bold', fontsize=16, pad=20, color = "black")
    for ax3 in axs[-1]:
        ax3.set_xlabel('Residue', fontsize=16)
        ax3.tick_params(axis='x', labelsize=16)
        #axs.set_xticklabels([f"{i+1}" for i in xticks])
        ax3.set_xticks(xticks)
        ax3.set_xticklabels([f"{i+1}" for i in xticks])
    for row in axs:
        for a in row[1:]:  # Skip the first column
            a.set_yticks([])

    output_file = os.path.join(relax_folder, output)
    plt.savefig(output_file, dpi=600)
    plt.close()

correlation_map("LRAEcorrelationHELICALstart.csv", "Correlation_combined.png")



def plot_rog_density_landscape_all(data_path, output, same=False, axs=None):
	rog_values_all=[]
	count_max=[]

	if same == False:
		fig, axs = plt.subplots(5, 5, figsize=(25, 25))
		fig.subplots_adjust(left=0.05, right=0.95, top=0.90, bottom=0.05, wspace=0.09, hspace=0.09)

	for path in data_path:
		name='/'.join(path.split("/")[-3:-1])
		rep_name=name.split("/")[0]
		forcefield=name.split("/")[1]

		rog_values, counts, values = read_rog_values(path)
		counts=[(x - min(counts)) / (max(counts) - min(counts)) for x in counts]
		count_max.append(max(counts))

		if same==False:

			i = int(rep_name[-1]) - 1
			j = FORCEFIELDS.index(forcefield)


			axs[i, j].set_xticks([])
			axs[i, j].set_yticks([])
			axs[i, j].set_xlim(0.5, 7)	
			axs[i, j].plot(values, counts, color=color_map.get(forcefield, "black"))
			axs[i, j].fill_between(values, counts, color=color_map.get(forcefield, "black"), alpha=0.3)
			axs[i, j].axvline(x=statistics.mean(rog_values), color="black", linestyle='--')

			if name in Best_cases:
				autoAxis = axs[i, j].axis()
				print(autoAxis)
				rec = Rectangle((autoAxis[0]-(autoAxis[1] - autoAxis[0]) * 0.04, autoAxis[2]- (autoAxis[3] - autoAxis[2]) * 0.04), 
					(autoAxis[1] - autoAxis[0])*1.08, 
					(autoAxis[3] - autoAxis[2])*1.08, 
					fill=False, linewidth=5, edgecolor='black')
				axs[i, j].add_patch(rec)
				rec.set_clip_on(False)

			title_avg = str(round_sig(statistics.mean(rog_values))) + " nm"
			legend_item = Line2D([0], [0], lw=0, color='w', label=f'Avg = $\\mathbf{{{str(title_avg)}}}$')
			axs[i, j].legend(handles=[legend_item], fontsize=25, loc='upper right', frameon=False)
			for ax1, replica in zip(axs[:,0], row_header):
				ax1.set_ylabel(f"$\\bf{{{replica}}}$", fontsize=30, rotation=90)
			for ax2, forcefield in zip(axs[0], FORCEFIELDS):
				ax2.set_title(forcefield, weight='bold', fontsize=30, pad=20, color = "black")
			for ax3 in axs[-1]:
				ax3.set_xlabel('RoG (nm)', fontsize=30)
				ax3.tick_params(axis='x', labelsize=27)
				#axs.set_xticklabels([f"{i+1}" for i in xticks])
				ax3.set_xticks(np.arange(1, 7, 2))
			for row in axs:
				for a in row[1:]:  # Skip the first column
					a.set_yticks([])
			
		else:
			rog_values_all.append(rog_values)
		
	
	if same == True:
		flattened_list=[item for sublist in rog_values_all for item in sublist]
		rounded_values = [round(value, 1) for value in flattened_list]
		rounded_counts = {value: rounded_values.count(value) for value in set(rounded_values)}
		sorted_items = sorted(rounded_counts.items())
		values = [value[0] for value in sorted_items]
		counts = [value[1] for value in sorted_items]


		plt.plot(values, counts)
		plt.axvline(x=statistics.mean(flattened_list), color="black", linestyle='--', label="Average RoG")
		plt.xlabel('Radius of Gyration (nm)')
		plt.title('Radius of gyration landscape for selected simulations')
		plt.legend()

	plt.savefig(relax_folder + output +'.png')
	plt.close()


plot_rog_density_landscape_all(glob.glob(SIM_DIR + "*/*/md*gyrate.xvg"), "density_landscape_plot", same=False)
plot_rog_density_landscape_all([glob.glob(SIM_DIR + i + "/md*gyrate.xvg")[0] for i in Best_cases], "Accepted_cases/best_rog_density_landscape_plot", same=True)


def dif_to_exp(input, output):

	used_labels = [[], [], []]
	fig, axs = plt.subplots(1, 3, figsize=(15, 6))
	for i, ff in enumerate(FORCEFIELDS):
		list=[[], [], []]
		for num, item in enumerate(input):
			if ff in item and item in input:
				R1_diff_unfiltered=extract_values_pandas(relax_data, num, 3, include_header=False)
				R1_diff = [value for value in R1_diff_unfiltered if not math.isnan(value)]
				rep_name=item.replace("/" + ff, "")
				try:
					if rep_name in used_labels[0]:
						axs[0].scatter(i, statistics.mean(R1_diff), zorder=5, color=color_list[int(rep_name[-1]) - 1])
					else:
						axs[0].scatter(i, statistics.mean(R1_diff), zorder=5, label=rep_name, color=color_list[int(rep_name[-1]) - 1])
						used_labels[0].append(rep_name)
				except:
					pass
				list[0].extend(R1_diff)	

				R2_diff_unfiltered=extract_values_pandas(relax_data, num, 6, include_header=False)
				R2_diff = [value for value in R2_diff_unfiltered if not math.isnan(value)]
				try:
					if rep_name in used_labels[1]:
						axs[1].scatter(i, statistics.mean(R2_diff), zorder=5, color=color_list[int(rep_name[-1]) - 1])
					else:
						axs[1].scatter(i, statistics.mean(R2_diff), zorder=5, label=rep_name, color=color_list[int(rep_name[-1]) - 1])
						used_labels[1].append(rep_name)
				except:
					pass
				list[1].extend(R2_diff)
			
				hetNOE_diff_unfiltered=extract_values_pandas(relax_data, num, 9, include_header=False)
				hetNOE_diff = [value for value in hetNOE_diff_unfiltered if not math.isnan(value)]
				try:
					if rep_name in used_labels[2]:
						axs[2].scatter(i, statistics.mean(hetNOE_diff), zorder=5, color=color_list[int(rep_name[-1]) - 1])
					else:
						axs[2].scatter(i, statistics.mean(hetNOE_diff), zorder=5, label=rep_name, color=color_list[int(rep_name[-1]) - 1])
						used_labels[2].append(rep_name)
				except:
					pass
				list[2].extend(hetNOE_diff)
		try:
			R1_avg=statistics.mean(list[0])
			axs[0].bar(i, R1_avg, label=ff)
		except: 
			pass
		try:
			R2_avg=statistics.mean(list[1])
			axs[1].bar(i, R2_avg, label=ff)
		except:
			pass
		try:
			hetNOE_avg=statistics.mean(list[2])
			axs[2].bar(i, hetNOE_avg, label=ff)
		except: 
			pass

	axs[0].set_title('R1 average differences', fontsize=16)
	axs[0].legend()
	axs[1].set_title('R2 average differences', fontsize=16)
	axs[1].legend()
	axs[2].set_title('hetNOE average differences', fontsize=16)
	axs[2].legend()


	fig.suptitle("Average differences between simulations and experimental data", fontsize=16)
	plt.tight_layout()
	plt.savefig(relax_folder + output + '.png')
	plt.close()

dif_to_exp(Names, 'Difference_to_experiment_plot')



timescale_data_avg=glob.glob(Unst_folder + "Ctimes_Coeffs.csv")


def create_timescale_scatter_plot(input, output):
    fig, axs = plt.subplots(5, 5, figsize=(15, 15))  # 5x5 grid, adjust as needed
    fig.subplots_adjust(left=0.1, right=0.85, top=0.95, bottom=0.05, wspace=0.07, hspace=0.07)


    for i, rep in enumerate(REPLICAS):
            for j, ff in enumerate(FORCEFIELDS):


                file = os.path.join(SIM_DIR, rep, ff, input)
                if os.path.exists(file):
                    df = pd.read_csv(file, header=None)

                    C_times_ns = df.iloc[:, 0].tolist()
                    Coeffs = df.iloc[:, 1].tolist()

                    x_vals, y_vals, weights = [], [], []
                    for line_idx in range(len(C_times_ns)):
                        if "Res" in str(C_times_ns[line_idx]):
                            x = [int(num) for num in re.findall(r'\d+', str(C_times_ns[line_idx]))][0]
                            y_num = 1
                            while line_idx + y_num < len(C_times_ns):
                                next_parts = C_times_ns[line_idx + y_num]
                                if "Res" not in next_parts:
                                    try:
                                        y = float(C_times_ns[line_idx + y_num])
                                        weight = float(Coeffs[line_idx + y_num])
                                        if y < 89.125:
                                            x_vals.append(x)
                                            y_vals.append(y)
                                            weights.append(weight * 300)
                                        y_num += 1
                                    except (IndexError, ValueError):
                                        break
                                else:
                                    break
                    
                    scaler = MinMaxScaler()
                    weights = scaler.fit_transform(np.array(weights).reshape(-1, 1)).flatten()

                    scatter = axs[i, j].scatter(x_vals, y_vals, c=weights, cmap='viridis', zorder=3)
                    axs[i, j].set_xticks([])
                    axs[i, j].tick_params(axis='y', labelsize=16, pad=10)
                    axs[i, j].set_xlim(0, res_nr)
                    axs[i, j].set_yscale('log')
                    axs[i, j].set_ylim(10**(-3), 10**2)
                    axs[i, j].set_yticks([10**(-3), 10**(-1), 10**1])


                    if rep + "/" + ff in Best_cases:
                        autoAxis = axs[i, j].axis()
                        print(autoAxis)
                        rec = Rectangle((autoAxis[0]-(autoAxis[1] - autoAxis[0]) * 0.03, autoAxis[2]*1.5**(-1)), 
                            (autoAxis[1] - autoAxis[0])*1.06, 
                            (autoAxis[1] - autoAxis[0])*2, 
                            fill=False, linewidth=5, edgecolor='black')
                        axs[i, j].add_patch(rec)
                        rec.set_clip_on(False)

                        

                        

    for ax1, replica in zip(axs[:,0], row_header):
        y_label="Timescales (ns)"
        ax1.set_ylabel(f"$\\bf{{{replica}}}$\n\n{y_label}", fontsize=16, rotation=90)
    for ax2, forcefield in zip(axs[0], FORCEFIELDS):
        color = color_map.get(forcefield.split("/")[-1], "black")
        ax2.set_title(forcefield, weight='bold', fontsize=16, pad=22, color = "black")
    for ax3 in axs[-1]:
        ax3.set_xlabel('Residue', fontsize=16)
        ax3.tick_params(axis='x', labelsize=16)
			#axs.set_xticklabels([f"{i+1}" for i in xticks])
        ax3.set_xticks(xticks)
        ax3.set_xticklabels([f"{i+1}" for i in xticks])
    for row in axs:
        for a in row[1:]:  # Skip the first column
            a.set_yticks([])


    cbar_ax = fig.add_axes([0.87, 0.1, 0.02, 0.75])  # Position for the colorbar
    norm = plt.Normalize(vmin=0, vmax=1)
    cbar = plt.colorbar(plt.cm.ScalarMappable(cmap='jet', norm=norm), cax=cbar_ax)
    cbar.ax.tick_params(labelsize=22)
    cbar.set_label('Weight Coefficient', rotation=270, labelpad=35, fontsize=25)



    output_file = os.path.join(relax_folder, output)
    plt.savefig(output_file, dpi=600)
    plt.close()	
		


#create_timescale_scatter_plot(timescale_data_avg, 0, axs=None)
#plt.savefig(best_cases_folder + 'Timescale_plot_avg.png')
#plt.close()

create_timescale_scatter_plot("Ctimes_Coeffs.csv", "Timescale_plot_all.png")





def plot_avg_rog_bar(input, output):
	used_labels_rog = []
	avg_rog_list = []

	forcefield_data = {ff: [] for ff in FORCEFIELDS}

	for item in input:
		for ff in FORCEFIELDS:
			if ff in item:
				rog_values, counts, values = read_rog_values(item)
				name = '/'.join(item.split("/")[-3:-1])
				rep_name = name.split("/")[0]
				forcefield_data[ff].append((rep_name, statistics.mean(rog_values), rog_values))

	rog_values_by_forcefield = {ff: [] for ff in FORCEFIELDS}

	for ff, data in forcefield_data.items():
		for entry in data:
			rog_values_by_forcefield[ff].extend(entry[2])


	sorted_forcefields = []

	for ff in FORCEFIELDS:
		if rog_values_by_forcefield[ff]:
			avg_rog = statistics.mean(rog_values_by_forcefield[ff])
			sorted_forcefields.append((ff, avg_rog))

#	sorted_forcefields.sort(key=lambda x: x[1], reverse=True)

	idx=0
	best_rog_avg = Unst_folder + "best_rog_avg.csv"
	with open(best_rog_avg, mode='w', newline='') as csv_file:
		writer = csv.writer(csv_file)
		writer.writerow(["Forcefield/Replica", "Average_RoG"])
		for ff, rog_avg in sorted_forcefields:
			#plt.bar(idx, rog_avg, label=ff, color=color_map.get(ff, "black"))

			for rep_name, data, _ in forcefield_data[ff]:
				plt.scatter(idx, data, zorder=5, color=color_map.get(ff, "black"))
				writer.writerow([rep_name + "/" + ff, data])
			idx += 1


	plt.legend()
	plt.title(SIM_DIR.split("/")[-2])
	plt.ylabel('Average Gyration Radius (nm)')
	plt.tight_layout()
	plt.savefig(relax_folder + output + '.png')
	plt.close()

plot_avg_rog_bar(sorted(glob.glob(SIM_DIR + "*/*/md*gyrate.xvg")), "Average_rog_plot")

rog_values_best_all=[]

for i in [glob.glob(SIM_DIR + i + "/md*gyrate.xvg")[0] for i in Best_cases]:
	rog_values, counts, values = read_rog_values(i)
	rog_values_best_all.extend(rog_values)

rounded_values = [round(value, 1) for value in rog_values_best_all]
rounded_counts = {value: rounded_values.count(value) for value in set(rounded_values)}
sorted_items = sorted(rounded_counts.items())
values = [value[0] for value in sorted_items]
counts = [value[1] for value in sorted_items]

with open(Unst_folder + "Best_rog_landscape.txt", 'w') as f:
	f.write(f"Rog_avg:\t{statistics.mean(rog_values_best_all):.5f}\n")
	for x, y in zip(values, counts):
		f.write(f"{x:.3f}\t{y:.5f}\n")




def plot_ensembles_images(input, output):
	fig, axs = plt.subplots(5, 5, figsize=(15, 15))
	for j, ff in enumerate(FORCEFIELDS): 
		for k, item in enumerate(input):	
			if ff in item:
				rep_name = item.split("/")[-3]		
				i = int(rep_name[-1]) - 1
				img = imread(item)
				axs[i, j].imshow(img)
				del img
				gc.collect()

				axs[i, j].set_title(f'{rep_name}/{ff}')
				axs[i, j].axis('off')
				
				if rep_name + "/" + ff in Best_cases:
					rect = Rectangle((0, 0), 1, 1, transform=axs[i, j].transAxes,
						linewidth=5, edgecolor='black', facecolor='none')
					axs[i, j].add_patch(rect) 

	plt.tight_layout()
	plt.savefig(f'{relax_folder}/{output}.png', dpi=300)
	plt.close('all')

def contact_map_avg(input, output):
	matrix = avg_csv(input, use_header=True, use_index=True)
	print(matrix)

	color_map = plt.imshow(matrix, cmap='jet', vmin=0, vmax=1, origin='lower')
	cbar = plt.colorbar(color_map)
	cbar.ax.tick_params(labelsize=16)
	cbar.set_label('Probability of contact (≤ 15.0 Å)', rotation=270, labelpad=25, fontsize=18)

	plt.tick_params(axis='both', which='major', labelsize=12)
	plt.xlabel("Residue", fontsize=16)
	plt.ylabel("Residue", fontsize=16)

	plt.tight_layout()
	plt.savefig(output + "Avg_contact.png", dpi=600)
	plt.close()

contact_map_avg([glob.glob(SIM_DIR + i + "/CA_prob_within_15.0A.csv")[0] for i in Best_cases], best_cases_folder)


def corr_map_avg(input, output):
	matrix = avg_csv(input)

	csv_output_path = os.path.join(output, "Avg_corr_matrix.csv")
	headers=create_amino_acid_pairs(protein_sequence)

	df = pd.DataFrame(matrix)
	print(df)
	df.columns = headers
	df.index = headers

	df.to_csv(csv_output_path, index=True)

	fig, ax = plt.subplots()
	color_map = plt.imshow(matrix,vmin=-1, vmax=1, origin='lower')
	color_map.set_cmap("seismic")

	plt.colorbar()
	plt.savefig(output + "Avg_corr.png", dpi=600)
	plt.close()

#corr_map_avg([glob.glob(SIM_DIR + i + "/*correlation*.csv")[0] for i in Best_cases], best_cases_folder)


def create_avg_distances_plot(input_file, output_file):
	matrix = avg_csv(input_file, use_header=True, use_index=True) 

	vmin = 0
	vmax = 20

	color_map = plt.imshow(matrix, vmin=vmin, vmax=vmax, cmap='jet', origin='lower')

	cbar = plt.colorbar(color_map)
	cbar.ax.tick_params(labelsize=16)
	cbar.set_label('Distance (Å)', rotation=270, labelpad=25, fontsize=18)

	plt.tick_params(axis='both', which='major', labelsize=12)
	plt.xlabel("Residue", fontsize=16)
	plt.ylabel("Residue", fontsize=16)

	plt.tight_layout()
	plt.savefig(output_file + "Avg_dist.png", dpi=600)
	plt.close()

#create_avg_distances_plot([glob.glob(SIM_DIR + i + "/CA_avg_distances.csv")[0] for i in Best_cases], best_cases_folder)

def run_pymol_operations():
	pdb_data = sorted(glob.glob(SIM_DIR + "replica*/*/"))	
	cmd.bg_color("white")
	cmd.set("ray_opaque_background", 1)
	path=[]
	for i in pdb_data:
		if i.split('/')[-3]+ "/" + i.split('/')[-2] in Best_cases:
			path.append(i)		
			name = i.split('/')[-4]

			md = glob.glob(i + 'md*noPBC.xtc')[0]
			temp = glob.glob(i + 'temp*gro')[0]

			cmd.load(temp, name)
			cmd.load_traj(md, name, state=1, interval=20000)
			cmd.hide('all')
			cmd.show('ribbon', name)
			cmd.color('green', name)
			cmd.ray(300, 300)
	cmd.set('all_states', 'on')
	cmd.intra_fit(f"{name}")
	cmd.center()
	cmd.zoom()

	pymol_path = best_cases_folder + 'Ensemble_' + SIM_DIR.split('/')[-2] + '_aligned_fig.png'
	second_path = best_cases_folder + 'Secondary probability_' + SIM_DIR.split('/')[-2] + '.png'
	merged_img = best_cases_folder + 'Ensemble_secondary_' + SIM_DIR.split('/')[-2] + '.png'	

	cmd.png(pymol_path, width=3000, height=1200, dpi=300)

	
	calculate_coil_probability(path, second_path)	
	
	merge_images_vertically(pymol_path, second_path, merged_img)
	
	cmd.delete('all')

	for case in pdb_data:	
		try:
			cmd.bg_color("white")
			cmd.set("ray_opaque_background", 1)

			rep_name = case.split('/')[-3]
			forcefield = case.split('/')[-2]
			selected = color_map.get(forcefield, 'black')

			md = glob.glob(case + 'md*ns_noPBC.xtc')[0]
			temp = glob.glob(case + 'temp*gro')[0]
	
			cmd.load(temp, rep_name + forcefield)
			cmd.load_traj(md, rep_name + forcefield, state=1, interval=20000)
			cmd.color(selected, rep_name + forcefield)
			cmd.set('all_states', 'on')
			cmd.ray(300, 300)
			cmd.intra_fit(f"{rep_name + forcefield}")
			cmd.center()
			cmd.zoom()

			pymol_path = case + 'Ensemble_' + rep_name + '_' +  forcefield + '_aligned_fig.png'
			cmd.png(pymol_path, width=3000, height=1200, dpi=300)     
			second_path = best_cases_folder + 'Secondary probability ' + rep_name + '_' +  forcefield + '.png'
			merged_img = best_cases_folder + 'Ensemble_secondary_' + rep_name + '_' +  forcefield + '.png'

			calculate_coil_probability(case, second_path)
			merge_images_vertically(pymol_path, second_path, merged_img)

			cmd.delete('all')
		except:
			pass
	
	
#run_pymol_operations()


plot_ensembles_images(sorted(glob.glob(SIM_DIR + "replica*/*/Contact_map.png")), "Contact_map_combined")
plot_ensembles_images(sorted(glob.glob(SIM_DIR + "replica*/*/Distance_map.png")), "Distance_map_combined")
plot_ensembles_images(sorted(glob.glob(SIM_DIR + "replica*/*/*correlation*.png")), "Correlation_combined")
plot_ensembles_images(sorted(glob.glob(SIM_DIR + "replica*/*/Ensemble*replica*aligned_fig.png")), "Ensembles_aligned_combined")


subprocess.run(["python3", avg_script])


best_cases_text = best_cases_folder + "best_cases_list.txt"
with open(best_cases_text, 'w') as file:
	for i in Best_cases:
		file.write(SIM_DIR + i + '\n')

	



cmd.quit()
































