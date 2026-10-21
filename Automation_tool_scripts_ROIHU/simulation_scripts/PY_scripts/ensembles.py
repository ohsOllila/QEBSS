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

three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}


FORCEFIELDS=["AMBER03WS", "AMBER99SB-DISP", "AMBER99SBWS", "DESAMBER", "CHARMM36M"]
REPLICAS=["replica_01", "replica_02", "replica_03", "replica_04", "replica_05"]

row_header =["Replica 1", "Replica 2", "Replica 3", "Replica 4", "Replica 5"]
color_list=['red', 'blue', 'green', 'purple', 'orange']

color_map = {
    'AMBER03WS': 'red',
    'AMBER99SB-DISP': 'blue',
    'AMBER99SBWS': 'green',
    'CHARMM36M': 'purple',
    'DESAMBER': 'orange'
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



def round_sig(x, sig=2):
	scale = sig - int(floor(log10(abs(x)))) - 1
	rounded_value = round(x, scale)
    
	return rounded_value



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
    plt.figure(figsize=(14, 3), dpi=300)
    plt.bar(range(1, len(average_coil_probabilities)+1), average_coil_probabilities, color='grey', width=0.9)
    #plt.xlabel('Residue', fontsize = 35)
    #plt.ylabel('Avg C', fontsize = 35)
    plt.xlim(0, len(seq_string)+1)
    plt.ylim(0, 1)
    #plt.xticks(ticks=range(len(average_coil_probabilities)), labels=[f"{i+1} {residue_names[i]}" for i in range(len(average_coil_probabilities))], rotation=90, fontsize=5)
    #plt.gca().set_xticks(xticks)
    plt.gca().set_xticks([])
    plt.gca().set_yticks([])
    plt.gca().tick_params(axis='both', labelsize=28)
    #plt.gca().set_xticklabels([f"{i+1}" for i in xticks], fontsize=28)
    #plt.gca().set_xticklabels([f"{i+1+412}" for i in xticks], fontsize=28)
    plt.tight_layout()

    plt.savefig(output)
    plt.close()

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
        
        
    Best_case_sum.append([case, int(round_sig(Ranking_sum)), int(round_sig(R1_rank_value)), int(round_sig(R2_rank_value)), int(round_sig(hetNOE_rank_value))])
        
    if round_sig(R1_rank_value) / 100 <= 1.5 and round_sig(R2_rank_value) / 100 <= 1.5 and round_sig(hetNOE_rank_value) / 100 <= 1.5:
        Best_cases.append(case)
        print(case)
    else:
        print(case + " do not meet criteria")




print(Best_case_sum)
if len(Best_cases)==0:
	lowest_case = min(Best_case_sum, key=lambda x: x[1])
	Best_cases.append(lowest_case[0])

print(Best_cases)


def plot_ensembles_images(input, output):
	fig, axs = plt.subplots(5, 5, figsize=(15, 15))
	fig.subplots_adjust(left=0.05, right=0.98, top=0.95, bottom=0.02, wspace=0.02, hspace=0.02)
	for j, ff in enumerate(FORCEFIELDS): 
		for k, item in enumerate(input):	
			if ff in item:
				rep_name = item.split("/")[-3]		
				i = int(rep_name[-1]) - 1
				img = imread(item)
				axs[i, j].imshow(img)
				del img
				gc.collect()

				axs[i, j].set_xticks([])
				axs[i, j].set_yticks([])
				
				if rep_name + "/" + ff in Best_cases:
					autoAxis = axs[i, j].axis()
					rec = Rectangle((autoAxis[0]-(autoAxis[1] - autoAxis[0]) * 0.01, autoAxis[2]- (autoAxis[3] - autoAxis[2]) * 0.01), 
						(autoAxis[1] - autoAxis[0])*1.02, 
						(autoAxis[3] - autoAxis[2])*1.02, 
						fill=False, linewidth=5, edgecolor='black')
					axs[i, j].add_patch(rec)
					rec.set_clip_on(False) 

	#data = np.zeros((10, 10))
	#axs[0, 0].imshow(data, cmap='gray', extent=(-5.0, 207.0, -10.765, 139.946), aspect='auto')
	#axs[0, 0].set_xticks([])
	#axs[0, 0].set_yticks([])
	for ax1, replica in zip(axs[:,0], row_header):
		ax1.set_ylabel(f'$\\bf{{{replica}}}$', fontsize=18, rotation=90)
	for ax2, forcefield in zip(axs[0], FORCEFIELDS):
		ax2.set_title(forcefield, weight='bold', fontsize=18, pad=22, color = "black")
	plt.savefig(f'{relax_folder}/{output}.png', dpi=300)
	plt.close('all')





def run_pymol_operations():
	pdb_data = sorted(glob.glob(SIM_DIR + "replica*/*/"))	
	cmd.bg_color("white")
	cmd.set("ray_opaque_background", 1)
	path=[]
	for i in pdb_data:
		if i.split('/')[-3]+ "/" + i.split('/')[-2] in Best_cases:
			path.append(i)		
			name = i.split('/')[-4]
			print(name)

			md = glob.glob(i + 'md*noPBC.xtc')[0]
			temp = glob.glob(i + 'temp*gro')[0]

			#cmd.load_traj(md, name, state=1, interval=20000)
			#cmd.hide('all')
			#cmd.show('ribbon', name)
			#cmd.color('green', name)
			#cmd.ray(300, 300)
	#cmd.set('all_states', 'on')
	#cmd.intra_fit(f"{name}")
	#cmd.center()
	#cmd.zoom()

	pymol_path = best_cases_folder + 'Ensemble_' + SIM_DIR.split('/')[-2] + '_aligned_fig.png'
	second_path = best_cases_folder + 'Secondary_probability_' + SIM_DIR.split('/')[-2] + '.png'
	merged_img = best_cases_folder + 'Ensemble_secondary_' + SIM_DIR.split('/')[-2] + '.png'	

	#cmd.png(pymol_path, width=3000, height=1200, dpi=300)

	
	calculate_coil_probability(path, second_path)	
	
	merge_images_vertically(pymol_path, second_path, merged_img)
	
	cmd.delete('all')

	for case in pdb_data:	
		try:
			#cmd.bg_color("white")
			#cmd.set("ray_opaque_background", 1)

			rep_name = case.split('/')[-3]
			forcefield = case.split('/')[-2]
			selected = color_map.get(forcefield, 'black')

			md = glob.glob(case + 'md*ns_noPBC.xtc')[0]
			temp = glob.glob(case + 'temp*gro')[0]
			print(md, temp)
	
			#cmd.load(temp, rep_name + forcefield)
			#cmd.load_traj(md, rep_name + forcefield, state=1, interval=20000)
			#cmd.color(selected, rep_name + forcefield)
			#cmd.set('all_states', 'on')
			#cmd.ray(300, 300)
			#cmd.intra_fit(f"{rep_name + forcefield}")
			#cmd.center()
			#cmd.zoom()

			pymol_path = case + 'Ensemble_' + rep_name + '_' +  forcefield + '_aligned_fig.png'
			#cmd.png(pymol_path, width=2000, height=1500, dpi=300)    
			#cmd.png(pymol_path, dpi=300)   
			second_path = case + 'Secondary_probability ' + rep_name + '_' +  forcefield + '.png'
			merged_img = case + 'Ensemble_secondary_' + rep_name + '_' +  forcefield + '.png'

			calculate_coil_probability([case], second_path)
			merge_images_vertically(pymol_path, second_path, merged_img)

			cmd.delete('all')
		except:
			pass
	
	
run_pymol_operations()



plot_ensembles_images(sorted(glob.glob(SIM_DIR + "replica*/*/Ensemble_secondary*.png")), "Ensembles_aligned_combined")



cmd.quit()
































