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
import pymol

three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}


FORCEFIELDS=["AMBER03WS", "AMBER99SB-DISP", "AMBER99SBWS", "CHARMM36M", "DESAMBER"]
color_list=['red', 'blue', 'green', 'purple', 'orange']

color_map = {
    'AMBER03WS': 'red',
    'AMBER99SB-DISP': 'blue',
    'AMBER99SBWS': 'green',
    'CHARMM36M': 'purple',
    'DESAMBER': 'orange'
}

color_map_sec = {'C': 'grey', 'H': 'red', 'E': 'blue'}



SIM_DIR = os.getcwd() + '/'
BASE_DIR = os.path.dirname(os.path.dirname(SIM_DIR)) + '/'
exp_data = os.path.join(BASE_DIR, os.path.basename(os.path.normpath(SIM_DIR)) + '_exp_data.txt')
Unst_folder = os.path.join(SIM_DIR.replace(os.path.basename(os.path.normpath(SIM_DIR)), 'results'), os.path.basename(os.path.normpath(SIM_DIR))) + '/'
relax_folder = os.path.join(Unst_folder, 'rep_to_exp_data/')
avg_script = os.path.join(BASE_DIR, "simulation_scripts", "PY_scripts", "average_contact.py")
py_script = os.path.join(BASE_DIR, "simulation_scripts", "PY_scripts", "Old_Relaxations_for_Samuli.py")
pymol_analysis = os.path.join(BASE_DIR, "simulation_scripts", "PY_scripts", "pymol_structure_analysis.py")
best_cases_folder = os.path.join(relax_folder, "Accepted_cases/")
relax_data = sorted(glob.glob(SIM_DIR + 'model*/*/relaxation_times.csv'))

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


with open("model_01.pdb", 'r') as file:
	lines = file.readlines()
	for i in range(len(lines)):
		parts = lines[i].split()
		if len(parts) > 3 and parts[0] == "ATOM" and parts[2] == "CA": 
			one_letter_sequence = convert_to_one_letter(parts[3])
			seq_string.append(one_letter_sequence)

protein_sequence = ''.join(seq_string)



xticks = []
for i in range(9, res_nr, 20):
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
	var_type=type(lst[0])

	temp_list = [float(x) for x in lst if x != "n"]
	temp_list.sort(key=abs, reverse=True)
	max_val = temp_list[:4] 

	max_idx=[lst.index(var_type((i))) for i in max_val]
    
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


def axs_plot(data_path, data_idx, column_nr, include_dev_points=False, include_header=True, axs=None, color=True):
	x_header, x_list = extract_values_pandas(data_path, data_idx, 0, include_header)
	y_header, y_list = extract_values_pandas(data_path, data_idx, column_nr, include_header)

	plot_settings = {
		"name": Names[data_idx],
		"rep_name": Names[data_idx].split('/')[0],
		"forcefield": Names[data_idx].split('/')[1],
		"parameter": y_header.split('_')[0],
		"ax_label": None,
		"label": None,
		"color": "black" if ("exp" in y_header or not color) else color_map.get(Names[data_idx].split('/')[1], "black")
	}

	if "Tau" in y_header:
		plot_settings.update({"ax_label": 'Effective correlation (ns)', "label": False})
	elif "exp" in y_header:
		plot_settings.update({"ax_label": f"{plot_settings['parameter']} {'relaxation rates (1/s)' if 'R1' in y_header or 'R2' in y_header else 'values'}", "label": False})
	elif "sim" in y_header:
		plot_settings.update({"ax_label": f"{plot_settings['parameter']} values", "label": False})
		if "results" in data_path[data_idx]:
			plot_settings["label"] = f"{plot_settings['parameter']} simulation average"
		else:
			dev_idx = deviating_points(extract_columns(data_path, data_idx, plot_settings["parameter"] + "_diff"))
			rank_label=Best_case_sum[data_idx][2 if plot_settings['parameter'] == "R1" else 3 if plot_settings['parameter'] == "R2" else 4 if plot_settings['parameter'] == "hetNOE" else None]

	if not color:
		if "exp" in y_header:
			plot_settings["label"] = f"{plot_settings['parameter']} experimental values"
		else:
			plot_settings["color"] = color_list[int(plot_settings["rep_name"][-1]) - 1]
			plot_settings["label"] = plot_settings["name"]

	df = pd.DataFrame({x_header: x_list, y_header: y_list}).dropna()
	df.plot(x=x_header, y=y_header, ax=axs, label=plot_settings["label"], marker='o', linestyle="-", lw=3.0, markersize=5, color=plot_settings["color"])

	if include_dev_points==True:
		try:
#			dev_idx = [10, 11, 12, 13, 22, 23, 24, 25, 26, 27]
#			max_points_x = [ x_list[i] for i in dev_idx ]
#			max_points_y = [ y_list[i] for i in dev_idx ]
			max_points_x = [ x_list[i] for i in dev_idx ]
			max_points_y = [ y_list[i] for i in dev_idx ]
			axs.scatter(max_points_x, max_points_y, color='black', s=40)
		except:
			pass

	if axs:
		axs.set_xlim(0-0.5, res_nr+0.5)
		axs.set_xticks(xticks)
		axs.set_xticklabels([f"{i+1}" for i in xticks])
		axs.set_xlabel('Residue number', fontsize = 18 if "results" not in data_path[data_idx] or not color else 10)
		axs.set_ylabel(plot_settings["ax_label"], fontsize=18 if "results" not in data_path[data_idx] or not color else 10)

		axs.tick_params(axis='x', labelsize=16 if "results" not in data_path[data_idx] or not color else 10)
		axs.tick_params(axis='y', labelsize=16 if "results" not in data_path[data_idx] or not color else 10)

		if plot_settings["label"]==False:
			axs.legend().remove()
		else:
			axs.legend(fontsize=15 if "results" not in data_path[data_idx] and color else 11)
		

		if "results" not in data_path[data_idx]:
			try:
				axs.set_title(f'{plot_settings["name"]}\nRank value = {str(rank_label)}%', fontweight='bold', fontsize=18)
			except:
				axs.set_title(plot_settings["name"], fontweight='bold', fontsize=18)

	else:
		plt.legend('', frameon=False)
		plt.xticks(xticks)
		plt.xlabel('Residue number')
		plt.ylabel(plot_settings["ax_label"])

	if any(keyword in y_header for keyword in ["hetNOE_exp", "R1_exp", "R2_exp"]) and "results" not in data_path[data_idx]:
		y_min, y_max = find_largest_value(y_header, y_header.replace("exp", "sim"))
		if "hetNOE" in y_header or "R1" in y_header:
			interval = 0.5
		elif "R2" in y_header:
			interval = 5
		else:
			interval = 1
    
		yticks = np.arange(round_sig(y_min) if "hetNOE" in y_header else 0, y_max, interval)
		axs.set_yticks(yticks)
		axs.set_ylim(y_min, y_max)

def secondary_structure_bar(GRO, XTC):
	traj = md.load(XTC, top=GRO, stride = 1000)

	initial_secondary_structure = md.compute_dssp(traj, simplified=True)[0]
	final_secondary_structure = md.compute_dssp(traj, simplified=True)[-1]

	del traj
	gc.collect()

	print(initial_secondary_structure)
	initial_colors = [color_map_sec[structure] for structure in initial_secondary_structure]
	final_colors = [color_map_sec[structure] for structure in final_secondary_structure]

	fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(18, 3))

	for i, color in enumerate(initial_colors):
		ax1.bar(i, 1, color=color)

	for i, color in enumerate(final_colors):
		ax2.bar(i, 1, color=color)

	for ax in [ax1, ax2]:
		ax.set_xlim(0, len(initial_secondary_structure))
		ax.set_ylim(0, 1)
		ax.set_yticks([])
		ax.set_xticks([])

	ax1.set_title('Initial Secondary Structure')
	ax2.set_title('Final Secondary Structure')

	buf = io.BytesIO()
	fig.savefig(buf, format='png', bbox_inches='tight')
	plt.close(fig)
	buf.seek(0)
    

	secondary_img = Image.open(buf)
	secondary_img_array = np.array(secondary_img)
	
    
	return secondary_img_array

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
	csvwriter = csv.writer(csvfile)
	#csvwriter.writerow(["Force field", "Replica", "R1 RMSRE", "R1 (%)", "R2 RMSRE", "R2 (%)", "hetNOE RMSRE", "hetNOE (%)", "Sum (%)"])
	csvwriter.writerow(["Force field", "Replica", "R1 RMSD", "R1 (%)", "R2 RMSD", "R2 (%)", "hetNOE RMSD", "hetNOE (%)", "Sum (%)"])
	for i in range(len(Names)):
			case=Names[i]
			rep_name=case.split("/")[0]
			ff_name=case.split("/")[1]
			#RMSRE_R1=calculate_rmsre(i, "R1_exp", "R1_diff")
			#RMSRE_R2=calculate_rmsre(i, "R2_exp", "R2_diff")
			#RMSRE_hetNOE=calculate_rmsre(i, "hetNOE_exp", "hetNOE_diff")

			#R1_rank_value=(float(RMSRE_R1)/ranking_value("R1_exp", "R1_diff"))*100
			#R2_rank_value=(float(RMSRE_R2)/ranking_value("R2_exp", "R2_diff"))*100
			#hetNOE_rank_value=(float(RMSRE_hetNOE)/ranking_value("hetNOE_exp", "hetNOE_diff"))*100
			RMSD_R1=calculate_rmsd(i, "R1_diff")
			RMSD_R2=calculate_rmsd(i, "R2_diff")
			RMSD_hetNOE=calculate_rmsd(i, "hetNOE_diff")


			R1_rank_value=(float(RMSD_R1)/ranking_value("R1_diff"))*100
			R2_rank_value=(float(RMSD_R2)/ranking_value("R2_diff"))*100
			hetNOE_rank_value=(float(RMSD_hetNOE)/ranking_value("hetNOE_diff"))*100
			Ranking_sum=R1_rank_value+R2_rank_value+hetNOE_rank_value
			csvwriter.writerow([ff_name, rep_name, round_sig(RMSD_R1), int(round_sig(R1_rank_value)), round_sig(RMSD_R2), int(round_sig(R2_rank_value)), round_sig(RMSD_hetNOE), int(round_sig(hetNOE_rank_value)), int(round_sig(Ranking_sum))])
			#csvwriter.writerow([ff_name, rep_name, round(RMSRE_R1, 2), round(R1_rank_value), round(RMSRE_R2, 2), round(R2_rank_value), round(RMSRE_hetNOE, 2), round(hetNOE_rank_value), round(Ranking_sum)])
			Best_case_sum.append([case, int(round_sig(Ranking_sum)), int(round_sig(R1_rank_value)), int(round_sig(R2_rank_value)), int(round_sig(hetNOE_rank_value))])
			if round_sig(R1_rank_value)/100 <= 1.5 and round_sig(R2_rank_value)/100 <= 1.5 and round_sig(hetNOE_rank_value)/100 <= 1.5:
				Best_cases.append(case)
				print(case)
			else:
				print(case + " do not meet criteria")
print(Best_case_sum)
if len(Best_cases)==0:
	lowest_case = min(Best_case_sum, key=lambda x: x[1])
	Best_cases.append(lowest_case[0])

print(Best_cases)
csv_to_pdf_with_conditional_formatting(ranking_file, value_min_threshold = 100, value_max_threshold=150)

os.makedirs(Unst_folder + "correlation_functions/", exist_ok=True)
'''
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
'''
avg_path=glob.glob(Unst_folder + "relaxation_times.csv")


os.chdir(SIM_DIR)

for item in FORCEFIELDS:
	fig, axs = plt.subplots(1, 3, figsize=(15, 6))
	for i, name in enumerate(Names):
		if name.split("/")[1]==item:
			axs_plot(relax_data, i, 2, include_header=True, axs=axs[0], color=False)
			axs_plot(relax_data, i, 5, include_header=True, axs=axs[1], color=False)
			axs_plot(relax_data, i, 8, include_header=True, axs=axs[2], color=False)
	else:
		axs_plot(relax_data, i, 1, include_header=True, axs=axs[0], color=False)
		axs_plot(relax_data, i, 4, include_header=True, axs=axs[1], color=False)
		axs_plot(relax_data, i, 7, include_header=True, axs=axs[2], color=False)
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



fig, axs = plt.subplots(3, 1, figsize=(15, 6))
axs_plot(avg_path, 0, 2, include_header=True, axs=axs[0])
axs_plot(avg_path, 0, 5, include_header=True, axs=axs[1])
axs_plot(avg_path, 0, 8, include_header=True, axs=axs[2])
axs_plot(avg_path, 0, 1, include_header=True, axs=axs[0])
axs_plot(avg_path, 0, 4, include_header=True, axs=axs[1])
axs_plot(avg_path, 0, 7, include_header=True, axs=axs[2])
plt.tight_layout()
plt.savefig(relax_folder + 'Accepted_cases/average_relaxation_compaired_plot.png', dpi=300)
plt.close('all')


axs_plot(avg_path, 0, 10, include_header=True, axs=None)
plt.savefig(relax_folder + 'Accepted_cases/Tau_effective_area_avg.png')
plt.close('all')


def relaxation_combined(sim, exp, output, execution):
	fig, axs = plt.subplots(5, 5, figsize=(25, 25))
	for j, ff in enumerate(FORCEFIELDS):
		for idx, item in enumerate(Names):
			if ff in item:
				rep_name=item.split("/")[0]
				i = int(rep_name[-1]) - 1
				execution(relax_data, idx, sim, include_dev_points=True, include_header=True, axs=axs[i, j])
				if exp is not None:
					execution(relax_data, idx, exp, include_dev_points=True, include_header=True, axs=axs[i, j])
				if rep_name + "/" + ff in Best_cases:
					rect = Rectangle((0, 0), 1, 1, transform=axs[i, j].transAxes,
						linewidth=5, edgecolor='green', facecolor='None')
					axs[i, j].add_patch(rect)
	plt.tight_layout(pad=3.0)
	plt.savefig(relax_folder + output +'.png')
	plt.close()

relaxation_combined(1, 2, 'R1_relaxation_combined_plot', axs_plot)

relaxation_combined(4, 5, 'R2_relaxation_combined_plot', axs_plot)
relaxation_combined(7, 8, 'hetNOE_relaxation_combined_plot', axs_plot)

relaxation_combined(10, None, "Tau_effective_area_all", axs_plot)




def plot_rog_density_landscape_all(data_path, output, same=False, axs=None):
	rog_values_all=[]
	count_max=[]

	if same == False:
		fig, axs = plt.subplots(5, 5, figsize=(20, 20))

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


			axs[i, j].set_xticks(np.arange(0.5, 7, 1))
			axs[i, j].tick_params(axis='x', labelsize=18)
			axs[i, j].set_yticks([])
			axs[i, j].set_xlim(0.5, 7)	
			axs[i, j].plot(values, counts, color=color_map.get(forcefield, "black"))
			axs[i, j].fill_between(values, counts, color=color_map.get(forcefield, "black"), alpha=0.3)
			axs[i, j].axvline(x=statistics.mean(rog_values), color="black", linestyle='--')

			if rep_name + "/" + forcefield in Best_cases:
				rect = Rectangle((0, 0), 1, 1, transform=axs[i, j].transAxes,
					linewidth=5, edgecolor='black', facecolor='none')
				axs[i, j].add_patch(rect) 
			axs[i, j].set_xlabel('Radius of Gyration (nm)', fontsize=18)
			title_avg = "avg = " + str(round_sig(statistics.mean(rog_values)))
			axs[i, j].set_title(f"{name}\n{title_avg}", fontweight='bold', fontsize=18)
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

	plt.tight_layout(pad=4.0)
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


timescale_data=sorted(glob.glob(SIM_DIR + 'model*/*/Ctimes_Coeffs.csv'))
timescale_data_avg=glob.glob(Unst_folder + "Ctimes_Coeffs.csv")

def create_timescale_scatter_plot(path, idx, axs=None):
	C_times_ns=extract_values_pandas_list(path, idx, 0)
	Coeffs=extract_values_pandas_list(path, idx, 1)
	name='/'.join(path[idx].split("/")[-3:-1])
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
							weights.append(weight * 100)
						y_num += 1
					except (IndexError, ValueError):
						break
				else:
					break


	if axs is None:
		plt.scatter(x_vals, y_vals, s=weights, zorder=5, color='red')
		plt.xlabel('Residue number', fontsize=14)
		plt.ylabel('Timescales (ns)', fontsize=14)
		plt.title('Timescale Scatter Plot', fontsize=14)
	else:
		axs.scatter(x_vals, y_vals, s=weights, zorder=5, color='red')
		axs.set_xlabel('Residue number', fontsize=14)
		axs.set_ylabel('Timescales (ns)', fontsize=14)
		if name in Best_cases:
			rect = Rectangle((0, 0), 1, 1, transform=axs.transAxes,
				linewidth=5, edgecolor='black', facecolor='None')
			axs.add_patch(rect)
		axs.set_yticks(np.arange(0, 81, 20)) 
		axs.set_ylim(0, 80)
		#axs.set_xticks(xticks) 
		axs.set_xlim(0, res_nr)
		axs.set_title('/'.join(path[idx].split("/")[-3:-1]), fontweight='bold')


create_timescale_scatter_plot(timescale_data_avg, 0, axs=None)
plt.savefig(best_cases_folder + 'Timescale_plot_avg.png')
plt.close()


plot_all(timescale_data, "Timescale_plot_all.png", create_timescale_scatter_plot)





def plot_avg_rog_bar(input, output):
	used_labels_rog = []
	for i, ff in enumerate(FORCEFIELDS):
		list=[]
		for item in input:
			if ff in item:
				rog_values, counts, values = read_rog_values(item)

				name='/'.join(item.split("/")[-3:-1])
				rep_name=name.split("/")[0]


				try:
					if rep_name in used_labels_rog:
						plt.scatter(i, statistics.mean(rog_values), zorder=5, color=color_list[int(rep_name[-1]) - 1])
					else:
						plt.scatter(i, statistics.mean(rog_values), zorder=5, label=rep_name, color=color_list[int(rep_name[-1]) - 1])
						used_labels_rog.append(rep_name)
				except:
					pass
				list.extend(rog_values)
		try:
			Rog_avg=statistics.mean(list)
			plt.bar(i, Rog_avg, label=ff)
			plt.legend()
		except:
			pass
	plt.title(SIM_DIR.split("/")[-2] )
	plt.ylabel('Average Gyration Radius (nm)')
	plt.tight_layout()
	plt.savefig(relax_folder + output + '.png')
	plt.close()

plot_avg_rog_bar(glob.glob(SIM_DIR + "*/*/md*gyrate.xvg"), "Average_rog_plot")

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
				if "align" in item:
					path='/'.join(item.split("/")[0:-1])
					
					XTC = glob.glob(path + '/md*noPBC.xtc')[0]
					GRO = glob.glob(path + '/temp*gro')[0]
					print(GRO)
					ensemble_img = imread(item)
					axs[i, j].imshow(ensemble_img)

					secondary_img_array = secondary_structure_bar(GRO, XTC)
					bar_height = secondary_img_array.shape[0] / ensemble_img.shape[0] * 0.3  
					inset_ax = axs[i, j].inset_axes([0, -bar_height, 1, bar_height], transform=axs[i, j].transAxes)
					inset_ax.imshow(secondary_img_array)
					inset_ax.axis('off')
					del secondary_img_array
					del ensemble_img
					gc.collect()

				else:
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
	if "align" in input[0]:
		legend_labels = ['Coil (C)', 'Helix (H)', 'Extended (E)']
		handles = [plt.Line2D([0], [0], color=color_map_sec[code], lw=4) for code in color_map_sec]
		fig.subplots_adjust(bottom=0.3, wspace=0.33)
		axs[4, 2].legend(handles, legend_labels, loc='upper center', 
			bbox_to_anchor=(0.5, -0.2),fancybox=False, shadow=False, ncol=3)

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

corr_map_avg([glob.glob(SIM_DIR + i + "/*correlation*.csv")[0] for i in Best_cases], best_cases_folder)


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

create_avg_distances_plot([glob.glob(SIM_DIR + i + "/CA_avg_distances.csv")[0] for i in Best_cases], best_cases_folder)

def run_pymol_operations():
	pdb_data = sorted(glob.glob(SIM_DIR + "model*/*/"))	
	cmd.bg_color("white")
	cmd.set("ray_opaque_background", 1)

	for i in pdb_data:
		if i.split('/')[-3]+ "/" + i.split('/')[-2] in Best_cases:
			name = i.split('/')[-4]

			md = glob.glob(i + 'md*noPBC.xtc')[0]
			temp = glob.glob(i + 'temp*gro')[0]

			cmd.load(temp, name)
			cmd.load_traj(md, name, state=1, interval=20000)
			cmd.color('green', name)
			cmd.ray(300, 300)
	cmd.set('all_states', 'on')
	cmd.intra_fit(f"{name}")
	cmd.center()
	cmd.zoom()

	cmd.png(best_cases_folder + 'Ensemble_' + SIM_DIR.split('/')[-2] + '_aligned_fig.png')	
	cmd.delete('all')
	
	values=extract_values_pandas_list(avg_path, 0, 10)
	max_val_res=values.index(max(values))+1
	min_res=max_val_res-2  
	max_res=max_val_res+2
	
#	print(min_res)
#	print(max_res)
	for path in pdb_data:
		cmd.bg_color("white")
		cmd.set("ray_opaque_background", 1)

		rep_name = path.split('/')[-3]
		forcefield = path.split('/')[-2]
		selected = color_map.get(forcefield, 'black')

		md = glob.glob(path + 'md*ns_noPBC.xtc')[0]
		temp = glob.glob(path + 'temp*gro')[0]

		if len(md) > 0 and len(temp) > 0:
			cmd.load(temp, rep_name + forcefield)
			cmd.load_traj(md, rep_name + forcefield, state=1, interval=20000)
			cmd.color(selected, rep_name + forcefield)
			cmd.set('all_states', 'on')
			#cmd.color("yellow", f"{rep_name + forcefield} and resi {min_res}-{max_res}")
			#cmd.intra_fit(f"{rep_name + forcefield} and resi {min_res}-{max_res}")
			cmd.ray(300, 300)
			cmd.intra_fit(f"{rep_name + forcefield}")
		cmd.center()
		cmd.zoom()

		cmd.png(path + 'Ensemble_' + rep_name + '_' +  forcefield + '_aligned_fig.png', dpi=300)	
		cmd.delete('all')

	
	
run_pymol_operations()

plot_ensembles_images(sorted(glob.glob(SIM_DIR + "model*/*/Contact_map.png")), "Contact_map_combined")
plot_ensembles_images(sorted(glob.glob(SIM_DIR + "model*/*/Distance_map.png")), "Distance_map_combined")
plot_ensembles_images(sorted(glob.glob(SIM_DIR + "model*/*/*correlation*.png")), "Correlation_combined")
plot_ensembles_images(sorted(glob.glob(SIM_DIR + "model*/*/Ensemble*model*aligned_fig.png")), "Ensembles_aligned_combined")

subprocess.run(["python3", avg_script])


cmd.quit()



























