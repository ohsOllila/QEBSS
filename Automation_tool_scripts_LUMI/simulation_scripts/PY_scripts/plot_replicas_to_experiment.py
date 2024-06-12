#!/usr/bin/python3


import matplotlib.pyplot as plt
import os
import shutil
import fileinput
import glob
import math
import csv
import subprocess
import re
import numpy as np
from collections import Counter
from matplotlib.image import imread
import pandas as pd
import statistics
import random
import MDAnalysis as mda
from pymol import cmd
from matplotlib.patches import Rectangle
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
    'model_01': 'red',
    'model_02': 'blue',
    'model_03': 'green',
    'model_04': 'purple',
    'model_05': 'orange'
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



xticks=[1]
for i in range(10, res_nr, 20):
	xticks.append(i)


Names=['/'.join(i.split("/")[-3:-1]) for i in relax_data]




def extract_columns(csv_folder):
	result = {}
	for filename in csv_folder:
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
    
	return result

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

def extract_data(files, sim_case_nr, parameter):
	spec_key=list(files.keys())[sim_case_nr]
	parameter = files.get(spec_key, {}).get(parameter, [])

	return parameter



def use_list(files, sim_case_nr, parameter):
	parameter=extract_data(files, sim_case_nr, parameter)
	cleaned_list = [float(x) for x in parameter if x != 'n']
	
	return cleaned_list


def ranking_value(files, parameter):
	RMSD_lists=[]
	for i in range(len(list(files.keys()))):
		RMSD_lists.append((calculate_rmsd(use_list(files, i, parameter))))  
	rank=min(float(i) for i in RMSD_lists)

	return rank


def find_largest_value(data_path, column_number, include_header=False):
	ymax_list=[]
	ymin_list=[]
	for i in range(len(Names)):
		ymax_list.append(np.nanmax(extract_values_pandas(data_path, i, column_number, include_header=False)))
		ymax_list.append(np.nanmax(extract_values_pandas(data_path, i, column_number+1, include_header=False)))
		ymin_list.append(np.nanmin(extract_values_pandas(data_path, i, column_number, include_header=False)))
		ymin_list.append(np.nanmin(extract_values_pandas(data_path, i, column_number+1, include_header=False)))
	ymax=max(ymax_list)
	ymin=min(ymin_list)

	return ymin, ymax




def calculate_rmsd(values_diff):
	if values_diff:
		return math.sqrt(sum(float(x) ** 2 for x in values_diff) / len(values_diff))
	else:
		return None

def plot_all(file_path, output_file_name, execution):
	fig, axs = plt.subplots(5, 5, figsize=(15, 15))
	for ff_idx, ff in enumerate(FORCEFIELDS):
		for idx in range(len(file_path)):
			rep_name=list(extract_columns(relax_data).keys())[idx].split("/")[0]
			rep_idx=int(rep_name.split("/")[0][-1])-1
			if ff in file_path[idx]:
				execution(file_path, idx, axs[ff_idx, rep_idx])
	plt.tight_layout()
	plt.savefig(relax_folder + output_file_name)
	plt.close()


def avg_csv(input):
	dataframes = [pd.read_csv(file) for file in input]
	stacked_array = np.array([df.values for df in dataframes])
	mean_array = np.mean(stacked_array, axis=0)
	matrix = pd.DataFrame(mean_array)

	return matrix

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


def axs_plot(data_path, data_idx, column_nr, include_header=True, axs=None):
	x_header, x_list = extract_values_pandas(data_path, data_idx, 0, include_header)
	y_header, y_list = extract_values_pandas(data_path, data_idx, column_nr, include_header)

	name=Names[data_idx]
	rep_name=name.split('/')[0]	
	forcefield=name.split('/')[1]

	parameter=y_header.split('_')[0]

	df = pd.DataFrame({x_header: x_list, y_header: y_list})	
	df.dropna(inplace=True)


	if "hetNOE" in y_header:
		ax_label = parameter + ' values'
	elif "R1" in y_header or "R2" in y_header:
		ax_label = parameter + ' relaxation rates (1/s)'
	else:
		ax_label = 'Effective correlation times (ns)'

	if "exp" in y_header:
		selected = "black"
		tag = y_header.split('_')[0] + " experimental data"	
	else:
		selected = color_map.get(rep_name, "black")
		tag = name

	if "exp" in y_header:
		label = y_header.split('_')[0] + " experimental data"
	elif "results" in data_path[data_idx] and "sim" in y_header:
		label = y_header.split('_')[0] + " simulation average"
	else:
		label = tag


	df.plot(x=x_header, y=y_header, ax=axs, label=label, marker='o', linestyle="-", lw=1.0, markersize=2, color=selected)

	if "hetNOE_exp" in y_header and "results" not in data_path[data_idx]:
		y_min, y_max = find_largest_value(data_path, column_nr, include_header=False)
		axs.set_yticks(np.arange(y_min, y_max, 0.5))
		axs.set_ylim(y_min, y_max)
	elif "R1_exp" in y_header and "results" not in data_path[data_idx]:
		y_min, y_max = find_largest_value(data_path, column_nr, include_header=False)
		axs.set_yticks(np.arange(0, y_max, 0.5))
		axs.set_ylim(0, y_max)
	elif "R2_exp" in y_header and "results" not in data_path[data_idx]:
		y_min, y_max = find_largest_value(data_path, column_nr, include_header=False)
		axs.set_yticks(np.arange(0, y_max, 20))
		axs.set_ylim(0, y_max)


	if axs is not None:
		axs.set_xlim(0, res_nr)
		axs.set_xticks(xticks)
		axs.set_xlabel('Residue number')
		axs.set_ylabel(ax_label)
		axs.legend()
	else:
		plt.xticks(xticks)
		plt.xlabel('Residue number')
		plt.ylabel(ax_label)


RMSD_R1=[]
RMSD_R2=[]
RMSD_hetNOE=[]
for i in range(len(list(extract_columns(relax_data).keys()))):	
	RMSD_R1.append(calculate_rmsd(use_list(extract_columns(relax_data), i, "R1_diff")))
	RMSD_R2.append(calculate_rmsd(use_list(extract_columns(relax_data), i, "R2_diff")))
	RMSD_hetNOE.append(calculate_rmsd(use_list(extract_columns(relax_data), i, "hetNOE_diff")))




Best_cases=[]

ranking_file=relax_folder + 'ranking_table.csv'
with open(ranking_file, 'w', newline="") as csvfile:
	csvwriter = csv.writer(csvfile)
	csvwriter.writerow(["Force field", "Replica", "R1 RMSD", "R1 (%)", "R2 RMSD", "R2 (%)", "hetNOE RMSD", "hetNOE (%)", "Sum (%)"])
	for i in range(len(list(extract_columns(relax_data).keys()))):
			case=list(extract_columns(relax_data).keys())[i]
			rep_name=case.split("/")[0]
			ff_name=case.split("/")[1]
			RMSD_R1=calculate_rmsd(use_list(extract_columns(relax_data), i, "R1_diff"))
			RMSD_R2=calculate_rmsd(use_list(extract_columns(relax_data), i, "R2_diff"))
			RMSD_hetNOE=calculate_rmsd(use_list(extract_columns(relax_data), i, "hetNOE_diff"))

			R1_rank_value=(float(RMSD_R1)/ranking_value(extract_columns(relax_data), "R1_diff"))*100
			R2_rank_value=(float(RMSD_R2)/ranking_value(extract_columns(relax_data), "R2_diff"))*100
			hetNOE_rank_value=(float(RMSD_hetNOE)/ranking_value(extract_columns(relax_data), "hetNOE_diff"))*100
			Ranking_sum=R1_rank_value+R2_rank_value+hetNOE_rank_value
			csvwriter.writerow([ff_name, rep_name, RMSD_R1, R1_rank_value, RMSD_R2, R2_rank_value, RMSD_hetNOE, hetNOE_rank_value, Ranking_sum])
			if R1_rank_value/100 < 1.5 and R2_rank_value/100 < 1.5 and hetNOE_rank_value/100 < 1.5:
				Best_cases.append(case)
				print(case)
			else:
				print(case + " do not meet criteria")

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
avg_data=extract_columns(avg_path)

os.chdir(SIM_DIR)

for item in FORCEFIELDS:
	fig, axs = plt.subplots(1, 3, figsize=(15, 6))
	for i, name in enumerate(Names):
		if name.split("/")[1]==item:
			axs_plot(relax_data, i, 2, include_header=True, axs=axs[0])
			axs_plot(relax_data, i, 5, include_header=True, axs=axs[1])
			axs_plot(relax_data, i, 8, include_header=True, axs=axs[2])
	else:
		axs_plot(relax_data, i, 1, include_header=True, axs=axs[0])
		axs_plot(relax_data, i, 4, include_header=True, axs=axs[1])
		axs_plot(relax_data, i, 7, include_header=True, axs=axs[2])
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



fig, axs = plt.subplots(3, 1, figsize=(11, 8))
axs_plot(avg_path, 0, 2, include_header=True, axs=axs[0])
axs_plot(avg_path, 0, 5, include_header=True, axs=axs[1])
axs_plot(avg_path, 0, 8, include_header=True, axs=axs[2])
axs_plot(avg_path, 0, 1, include_header=True, axs=axs[0])
axs_plot(avg_path, 0, 4, include_header=True, axs=axs[1])
axs_plot(avg_path, 0, 7, include_header=True, axs=axs[2])
fig.suptitle(protein_sequence, fontsize=7)
plt.tight_layout()
plt.savefig(relax_folder + 'Accepted_cases/average_relaxation_compaired_plot.png')
plt.close('all')


axs_plot(avg_path, 0, 10, include_header=True, axs=None)
plt.savefig(relax_folder + 'Accepted_cases/Tau_effective_area_avg.png')
plt.close('all')


def relaxation_combined(sim, exp, output, execution):
	fig, axs = plt.subplots(5, 5, figsize=(30, 15))
	for j, ff in enumerate(FORCEFIELDS):
		for idx, item in enumerate(Names):
			if ff in item:
				rep_name=item.split("/")[0]
				i = int(rep_name[-1]) - 1
				execution(relax_data, idx, sim, include_header=True, axs=axs[j, i])
				if exp is not None:
					execution(relax_data, idx, exp, include_header=True, axs=axs[j, i])
				if rep_name + "/" + ff in Best_cases:
					rect = Rectangle((0, 0), 1, 1, transform=axs[j, i].transAxes,
						linewidth=5, edgecolor='black', facecolor='None')
					axs[j, i].add_patch(rect)
	plt.tight_layout()
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
		fig, axs = plt.subplots(5, 5, figsize=(30, 15))

	for path in data_path:
		name='/'.join(path.split("/")[-3:-1])
		rep_name=name.split("/")[0]
		forcefield=name.split("/")[1]

		rog_values, counts, values = read_rog_values(path)
		count_max.append(max(counts))
		if same==False:

			i = int(rep_name[-1]) - 1
			j = FORCEFIELDS.index(forcefield)

			axs[j, i].set_yticks(np.arange(0, 70000, 10000))
			axs[j, i].set_xticks(np.arange(0.5, 7, 1))
			axs[j, i].set_xlim(0.5, 7)
			axs[j, i].set_ylim(0, 70000)	
			axs[j, i].plot(values, counts, label=name)
			axs[j, i].axvline(x=statistics.mean(rog_values), color="black", linestyle='--', label="Average RoG")
			if rep_name + "/" + forcefield in Best_cases:
				rect = Rectangle((0, 0), 1, 1, transform=axs[j, i].transAxes,
					linewidth=5, edgecolor='black', facecolor='none')
				axs[j, i].add_patch(rect) 
			axs[j, i].set_xlabel('Count')
			axs[j, i].set_ylabel('Radius of Gyration (nm)')
			axs[j, i].legend()
		else:
			rog_values_all.append(rog_values)
	
	if same == True:
		flattened_list=[item for sublist in rog_values_all for item in sublist]
		rounded_values = [round(value, 1) for value in flattened_list]
		rounded_counts = {value: rounded_values.count(value) for value in set(rounded_values)}
		sorted_items = sorted(rounded_counts.items())
		values = [value[0] for value in sorted_items]
		counts = [value[1] for value in sorted_items]

		print(flattened_list)
		plt.plot(values, counts)
		plt.axvline(x=statistics.mean(flattened_list), color="black", linestyle='--', label="Average RoG")
		plt.xlabel('Count')
		plt.ylabel('Radius of Gyration (nm)')
		plt.title('Radius of gyration landscape for selected simulations')
		plt.legend()
	else:
		for ax in axs.flat:
			ax.set_yticks(np.arange(0, max(count_max), 10000))
			ax.set_ylim(0, max(count_max))

	plt.tight_layout()
	plt.savefig(relax_folder + output +'.png')
	plt.close()


plot_rog_density_landscape_all(glob.glob(SIM_DIR + "*/*/md*gyrate.xvg"), "density_landscape_plot", same=False)
#plot_rog_density_landscape_all([glob.glob(SIM_DIR + i + "/md*gyrate.xvg")[0] for i in Best_cases], "Accepted_cases/best_rog_density_landscape_plot", same=True)


def dif_to_exp(input, output):

	used_labels = [[], [], []]
	fig, axs = plt.subplots(1, 3, figsize=(15, 6))
	for i, ff in enumerate(FORCEFIELDS):
		list=[[], [], []]
		for num, item in enumerate(input):
			if ff in item and item in input:
				R1_diff_unfiltered=extract_values_pandas(relax_data, num, 3, include_header=False)
				R1_diff = [value for value in R1_diff_unfiltered if not math.isnan(value)]
				model=item.replace("/" + ff, "")
				try:
					if model in used_labels[0]:
						axs[0].scatter(i, statistics.mean(R1_diff), zorder=5, color=color_map.get(model, 'black'))
					else:
						axs[0].scatter(i, statistics.mean(R1_diff), zorder=5, label=model, color=color_map.get(model, 'black'))
						used_labels[0].append(model)
				except:
					pass
				list[0].extend(R1_diff)	

				R2_diff_unfiltered=extract_values_pandas(relax_data, num, 6, include_header=False)
				R2_diff = [value for value in R2_diff_unfiltered if not math.isnan(value)]
				try:
					if model in used_labels[1]:
						axs[1].scatter(i, statistics.mean(R2_diff), zorder=5, color=color_map.get(model, 'black'))
					else:
						axs[1].scatter(i, statistics.mean(R2_diff), zorder=5, label=model, color=color_map.get(model, 'black'))
						used_labels[1].append(model)
				except:
					pass
				list[1].extend(R2_diff)
			
				hetNOE_diff_unfiltered=extract_values_pandas(relax_data, num, 9, include_header=False)
				hetNOE_diff = [value for value in hetNOE_diff_unfiltered if not math.isnan(value)]
				try:
					if model in used_labels[2]:
						axs[2].scatter(i, statistics.mean(hetNOE_diff), zorder=5, color=color_map.get(model, 'black'))
					else:
						axs[2].scatter(i, statistics.mean(hetNOE_diff), zorder=5, label=model, color=color_map.get(model, 'black'))
						used_labels[2].append(model)
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

	axs[0].set_title('R1 averege differences')
	axs[0].legend()
	axs[1].set_title('R2 average differences')
	axs[1].legend()
	axs[2].set_title('hetNOE average differences')
	axs[2].legend()


	fig.suptitle("Average differences between simulations and experimental data")
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
		plt.xlabel('Residue number')
		plt.ylabel('Timescales (ns)')
		plt.title('Timescale Scatter Plot')
	else:
		axs.scatter(x_vals, y_vals, s=weights, zorder=5, color='red')
		axs.set_xlabel('Residue number')
		axs.set_ylabel('Timescales (ns)')
		if name in Best_cases:
			rect = Rectangle((0, 0), 1, 1, transform=axs.transAxes,
				linewidth=5, edgecolor='black', facecolor='None')
			axs.add_patch(rect)
		axs.set_yticks(np.arange(0, 81, 20)) 
		axs.set_ylim(0, 80)
		axs.set_xticks(xticks) 
		axs.set_xlim(0, res_nr)
		axs.set_title('/'.join(path[idx].split("/")[-3:-1]))


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
						plt.scatter(i, statistics.mean(rog_values), zorder=5, color=color_map.get(rep_name, 'black'))
					else:
						plt.scatter(i, statistics.mean(rog_values), zorder=5, label=rep_name, color=color_map.get(rep_name, 'black'))
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
	plt.title(SIM_DIR.split("/")[-1] )
	plt.ylabel('Average Gyration Radius (nm)')
	plt.tight_layout()
	plt.savefig(relax_folder + output + '.png')
	plt.close()


rog_values_best_all=[]
rog_values_best_round=[]
rog_counts_best=[]

for i in [glob.glob(SIM_DIR + i + "/md*gyrate.xvg")[0] for i in Best_cases]:
	rog_values, counts, values = read_rog_values(i)
	rog_values_best_all.extend(rog_values)
	rog_values_best_round.extend(values)
	rog_counts_best.extend(counts)



with open(Unst_folder + "Best_rog_landscape.txt", 'w') as f:
	f.write(f"Rog_avg:\t{statistics.mean(rog_values_best_all):.5f}\n")
	for x, y in zip(rog_values_best_round, rog_counts_best):
		f.write(f"{x:.3f}\t{y:.5f}\n")
def plot_ensembles_images(input, output):
	fig, axs = plt.subplots(5, 5, figsize=(15, 15))
	for j, ff in enumerate(FORCEFIELDS): # Use '*.csv' for CSV
		for k, item in enumerate(input):	
			if ff in item:
				ensemble = input[k]
				rep_name = ensemble.split("/")[-3]		
				i = int(rep_name[-1]) - 1
				img = imread(ensemble)
				axs[j, i].imshow(img)
				axs[j, i].set_title(f'{rep_name}/{ff}')
				axs[j, i].axis('off')

				if rep_name + "/" + ff in Best_cases:
					rect = Rectangle((0, 0), 1, 1, transform=axs[j, i].transAxes,
						linewidth=5, edgecolor='black', facecolor='none')
					axs[j, i].add_patch(rect) 

	plt.tight_layout()
	plt.savefig(f'{relax_folder}/{output}.png')
	plt.close()


pdb_data = sorted(glob.glob(SIM_DIR + "model*/*/"))

cmd.set("ray_opaque_background", 1)

for i in pdb_data:
	if i.split('/')[-3]+ "/" + i.split('/')[-2] in Best_cases:
		name = i.split('/')[-4]

		md = glob.glob(i + 'md*smooth*xtc')
		temp = glob.glob(i + 'temp*gro')

		if len(md) > 0 and len(temp) > 0:
			cmd.load(temp[0], name)
			cmd.load_traj(md[0], name, state=1, interval=2000)
			cmd.color('green', name)
			cmd.set('all_states', 'on')
			cmd.ray(300, 300)
obj = cmd.get_object_list('all')
for i in obj[1:]:
	cmd.align(obj[0], i, object='aln', transform=0)

cmd.png(best_cases_folder + 'Ensemble_' + SIM_DIR.split('/')[-2] + '_aligned_fig.png')	
cmd.delete('all')
cmd.quit()



for path in pdb_data:
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

plot_ensembles_images(sorted(glob.glob(SIM_DIR + "model*/*/*mdmat*.png")), "Contact_map_combined")
plot_ensembles_images(sorted(glob.glob(SIM_DIR + "model*/*/*correlation*.png")), "Correlation_combined")
plot_ensembles_images(sorted(glob.glob(SIM_DIR + "model*/*/Ensemble*model*aligned_fig.png")), "Ensembles_aligned_combined")

def contact_map_avg(input, output):
	matrix = avg_csv(input)
	fig, ax = plt.subplots()
	cmap = plt.cm.viridis
	plot = ax.imshow(matrix, cmap=cmap, origin="lower")
	
	#ax.set_xlabel("Residue",fontsize=15)
	#ax.set_ylabel("Residue",fontsize=15)

	#plt.xticks(np.arange(0, res_nr+1, 20), fontsize=15)
	#plt.yticks(np.arange(0, res_nr+1, 20), fontsize=15)

	cbar = plt.colorbar(plot, ax=ax)
	cbar.set_label("Distance (nm)",fontsize=15)
	plt.savefig(output + "Avg_contact.png", dpi=600)
	plt.close()

contact_map_avg([glob.glob(SIM_DIR + i + "/md*mdmat.csv")[0] for i in Best_cases], best_cases_folder)


def corr_map_avg(input, output):
	matrix = avg_csv(input)

	color_map = plt.imshow(matrix,vmin=-1, vmax=1, origin='lower')
	color_map.set_cmap("seismic")
	plt.colorbar()
	plt.savefig(output + "Avg_corr.png", dpi=600)
	plt.close()

corr_map_avg([glob.glob(SIM_DIR + i + "/*correlation*.csv")[0] for i in Best_cases], best_cases_folder)


subprocess.run(["python3", avg_script])
#subprocess.run(["python3", pymol_analysis])





