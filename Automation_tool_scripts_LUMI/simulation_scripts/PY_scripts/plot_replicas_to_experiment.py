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
import pymol


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

def extract_values_pandas(file_path, nr, column_number):
	df = pd.read_csv(file_path[nr])    
	column = df.iloc[:, column_number]
	return column


def filter_list(unfiltered_list):
        return [float(x) for x in unfiltered_list if x != 'n']

def use_list(files, sim_case_nr, parameter):
	spec_key=list(files.keys())[sim_case_nr]
	parameter = files.get(spec_key, {}).get(parameter, [])
	cleaned_list = filter_list(parameter)
	
	return cleaned_list


def adjust_existing_res(files, i, j):
    return [x for x, y in zip(files.get(list(files.keys())[i], {}).get("Residue_nr", []), files.get(list(files.keys())[i], {}).get(j, [])) if y != 'n']


def ranking_value(files, parameter):
	RMSD_lists=[]
	for i in range(len(list(files.keys()))):
		RMSD_lists.append((calculate_rmsd(use_list(files, i, parameter))))  
	rank=min(float(i) for i in RMSD_lists)

	return rank

def extreme_value(files, extreme, sim, exp):
	lists=[]
	for i in range(len(list(files.keys()))):
		lists.append(use_list(files, i, sim))
		lists.append(use_list(files, i, exp))
	ref=extreme(extreme(sublist) for sublist in lists)
	return ref
	


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


	


RMSD_R1=[]
RMSD_R2=[]
RMSD_NOE=[]
for i in range(len(list(extract_columns(relax_data).keys()))):	
	RMSD_R1.append(calculate_rmsd(use_list(extract_columns(relax_data), i, "R1_diff")))
	RMSD_R2.append(calculate_rmsd(use_list(extract_columns(relax_data), i, "R2_diff")))
	RMSD_NOE.append(calculate_rmsd(use_list(extract_columns(relax_data), i, "NOE_diff")))




Rog_data_list=[]
rog_values_best=[]
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
			RMSD_NOE=calculate_rmsd(use_list(extract_columns(relax_data), i, "NOE_diff"))

			R1_rank_value=(float(RMSD_R1)/ranking_value(extract_columns(relax_data), "R1_diff"))*100
			R2_rank_value=(float(RMSD_R2)/ranking_value(extract_columns(relax_data), "R2_diff"))*100
			NOE_rank_value=(float(RMSD_NOE)/ranking_value(extract_columns(relax_data), "NOE_diff"))*100
			Ranking_sum=R1_rank_value+R2_rank_value+NOE_rank_value
			csvwriter.writerow([ff_name, rep_name, RMSD_R1, R1_rank_value, RMSD_R2, R2_rank_value, RMSD_NOE, NOE_rank_value, Ranking_sum])
			if R1_rank_value/100 < 1.5 and R2_rank_value/100 < 1.5 and NOE_rank_value/100 < 1.5:
				Best_cases.append(case)
				print(Best_cases)
				with open(glob.glob(SIM_DIR+ case +'/md*gyrate.xvg')[0], 'r') as file:
					lines = file.readlines()
					for line in range(27, len(lines)):
						parts =(lines)[line].split()
						rog_values_best.append(round(float(parts[1]), 1))



if rog_values_best: 
	value_counts = Counter(rog_values_best)
	values, counts = zip(*sorted(value_counts.items(), key=lambda x: x[0]))
	plt.plot(values, counts)
	plt.xlabel('Radius of gyration')	
	plt.ylabel('Count')
	plt.title('Count of Radius of gyration values from best simulation extract_columns(relax_data)')
	plt.tight_layout()	
	plt.savefig(best_cases_folder + 'best_rog_density_landscape_plot.png')
	plt.close()
else:
	print("No simulation meets ranking criteria")


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
avg_data_path=glob.glob(Unst_folder + "relaxation_times.csv")
avg_data=extract_columns(avg_data_path)


os.chdir(SIM_DIR)

print(adjust_existing_res(extract_columns(relax_data), 0, "R1_sim"))
for item in FORCEFIELDS:
	fig, axs = plt.subplots(1, 3, figsize=(15, 6))
	for i, name in enumerate(list(extract_columns(relax_data).keys())):
		if name.split("/")[1]==item:
			rep_name=name.split("/")[0]
			selected = color_map.get(rep_name, 'black')
			axs[0].plot(adjust_existing_res(extract_columns(relax_data), i, "R1_sim"), use_list(extract_columns(relax_data), i, "R1_sim"), label='R1 Data_' + rep_name, marker='o', linestyle='-', lw=1.0, markersize=2, color=selected)
			axs[0].set_xlabel('Residue number')
			axs[0].set_ylabel('R1_values (1/s)')
			axs[0].set_title('R1_data')
			axs[0].legend()
			axs[1].plot(adjust_existing_res(extract_columns(relax_data), i, "R2_sim"), use_list(extract_columns(relax_data), i, "R2_sim"), label='R2 Data_' + rep_name, marker='o', linestyle='-', lw=1.0, markersize=2, color=selected)
			axs[1].set_xlabel('Residue number')
			axs[1].set_ylabel('R2_values (1/s)')
			axs[1].set_title('R2_data')
			axs[1].legend()
			axs[2].plot(adjust_existing_res(extract_columns(relax_data), i, "NOE_sim"), use_list(extract_columns(relax_data), i, "NOE_sim"), label='NOE Data_' + rep_name, marker='o', linestyle='-', lw=1.0, markersize=2, color=selected)
			axs[2].set_xlabel('Residue number')				
			axs[2].set_ylabel('NOE_values 1/s')
			axs[2].set_title('NOE_data')
			axs[2].legend()
	else:
		axs[0].plot(adjust_existing_res(extract_columns(relax_data), i, "R1_exp"), use_list(extract_columns(relax_data), i, "R1_exp"), label='R1 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
		axs[0].set_xlabel('Residue number')
		axs[0].set_ylabel('R1_values (1/s)')
		axs[0].set_xticks(np.arange(9, res_nr, 10))
		axs[0].set_xlim(0, res_nr)
		axs[0].set_title('R1_data')
		axs[0].legend()
		
		axs[1].plot(adjust_existing_res(extract_columns(relax_data), i, "R2_exp"), use_list(extract_columns(relax_data), i, "R2_exp"), label='R2 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
		axs[1].set_xlabel('Residue number')
		axs[1].set_ylabel('R2_values (1/s)')
		axs[1].set_xticks(np.arange(9, res_nr, 10))
		axs[1].set_xlim(0, res_nr)
		axs[1].set_title('R2_data')
		axs[1].legend()
	
		axs[2].plot(adjust_existing_res(extract_columns(relax_data), i, "NOE_exp"), use_list(extract_columns(relax_data), i, "NOE_exp"), label='NOE Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
		axs[2].set_xlabel('Residue number')
		axs[2].set_ylabel('NOE_values (1/s)')
		axs[2].set_xticks(np.arange(9, res_nr, 10))
		axs[2].set_xlim(0, res_nr)
		axs[2].set_title('NOE_data')
		axs[2].legend()
	
	plt.tight_layout()
	plt.savefig(relax_folder + item + '_plot.png')
	plt.close('all')

fig, axs = plt.subplots(1, 3, figsize=(15, 6))
for j in Best_cases:
	rep_name=j.split("/")[0]
	ff_name=j.split("/")[1]
	i=list(extract_columns(relax_data).keys()).index(j)

	selected = color_map.get(rep_name, 'black')
	axs[0].plot(adjust_existing_res(extract_columns(relax_data), i, "R1_sim"), use_list(extract_columns(relax_data), i, "R1_sim"), label='R1 Data_' + rep_name + '_' + ff_name, marker='o', linestyle='-', lw=1.0, markersize=2, color=selected)
	axs[0].set_xlabel('Residue number')
	axs[0].set_ylabel('R1_values')
	axs[0].set_title('R1_data')
	axs[0].legend()
	axs[1].plot(adjust_existing_res(extract_columns(relax_data), i, "R2_sim"), use_list(extract_columns(relax_data), i, "R2_sim"), label='R2 Data_' + rep_name + '_' + ff_name, marker='o', linestyle='-', lw=1.0, markersize=2, color=selected)
	axs[1].set_xlabel('Residue number')
	axs[1].set_ylabel('R2_values')
	axs[1].set_title('R2_data')
	axs[1].legend()
	axs[2].plot(adjust_existing_res(extract_columns(relax_data), i, "NOE_sim"), use_list(extract_columns(relax_data), i, "NOE_sim"), label='NOE Data_' + rep_name + '_' + ff_name, marker='o', linestyle='-', lw=1.0, markersize=2, color=selected)
	axs[2].set_xlabel('Residue number')				
	axs[2].set_ylabel('NOE_values')
	axs[2].set_title('NOE_data')
	axs[2].legend()
		
axs[0].plot(adjust_existing_res(extract_columns(relax_data), 0, "R1_exp"), use_list(extract_columns(relax_data), 0, "R1_exp"), label='R1 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
axs[0].set_xlabel('Residue number')
axs[0].set_ylabel('R1_values')
axs[0].set_xticks(np.arange(9, res_nr, 10))
axs[0].set_xlim(0, res_nr-20)
axs[0].set_title('R1_data')
axs[0].legend()

axs[1].plot(adjust_existing_res(extract_columns(relax_data), 0, "R2_exp"), use_list(extract_columns(relax_data), 0, "R2_exp"), label='R2 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
axs[1].set_xlabel('Residue number')
axs[1].set_ylabel('R2_values')
axs[1].set_xticks(np.arange(9, res_nr, 10))
axs[1].set_xlim(0, res_nr-20)
axs[1].set_title('R2_data')
axs[1].legend()

axs[2].plot(adjust_existing_res(extract_columns(relax_data), 0, "NOE_exp"), use_list(extract_columns(relax_data), 0, "NOE_exp"), label='NOE Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
axs[2].set_xlabel('Residue number')
axs[2].set_ylabel('NOE_values')
axs[2].set_xticks(np.arange(9, res_nr, 10))
axs[2].set_xlim(0, res_nr-20)
axs[2].set_title('NOE_data')
axs[2].legend()
	
plt.tight_layout()
plt.savefig(best_cases_folder + 'Best_cases_plot.png')
plt.close('all')

fig, axs = plt.subplots(3, 1, figsize=(11, 8))
axs[0].plot(adjust_existing_res(avg_data, 0, "R1_sim"), use_list(avg_data, 0, "R1_sim"), label='Simulation average', marker='o', linestyle='-', lw=1.0, markersize=2, color="blue")
axs[0].set_ylabel('R1 relaxation rate (1/s)')
axs[0].set_title(SIM_DIR.split("/")[-2])
axs[0].legend()
axs[1].plot(adjust_existing_res(avg_data, 0, "R2_sim"), use_list(avg_data, 0, "R2_sim"), label='Simulation average', marker='o', linestyle='-', lw=1.0, markersize=2, color="blue")
axs[1].set_ylabel('R2 relaxation rate (1/s)')
axs[1].legend()
axs[2].plot(adjust_existing_res(avg_data, 0, "NOE_sim"), use_list(avg_data, 0, "NOE_sim"), label='Simulation average', marker='o', linestyle='-', lw=1.0, markersize=2, color="blue")
axs[2].set_xlabel('Residue number')
axs[2].set_ylabel('hetNOE value (1/s)')
axs[2].legend()

axs[0].plot(adjust_existing_res(avg_data, 0, "R1_exp"), use_list(avg_data, 0, "R1_exp"), label='NMR experiment', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
axs[0].set_xticks(np.arange(9, res_nr, 10))
axs[0].set_xlim(0, res_nr-20)
axs[0].legend()
axs[1].plot(adjust_existing_res(avg_data, 0, "R2_exp"), use_list(avg_data, 0, "R2_exp"), label='NMR experiment', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
axs[1].set_xticks(np.arange(9, res_nr, 10))
axs[1].set_xlim(0, res_nr-20)
axs[1].legend()
axs[2].plot(adjust_existing_res(avg_data, 0, "NOE_exp"), use_list(avg_data, 0, "NOE_exp"), label='NMR experiment', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
axs[2].set_xticks(np.arange(9, res_nr, 10))
axs[2].set_xlim(0, res_nr-20)
axs[2].legend()

plt.tight_layout()
plt.savefig(relax_folder + 'Accepted_cases/average_relaxation_compaired_plot.png')
plt.close('all')


def R1_relaxation_combined(input, output):
	if len(input) > 10:
		fig, axs = plt.subplots(5, 5, figsize=(30, 15))
		for j, ff in enumerate(FORCEFIELDS):
			for item in input:
				if ff in item:
					index=input.index(item) 
					i=int(item.split("/")[0][-1])-1
					rep_name=item.split("/")[0]
					axs[j, i].plot(adjust_existing_res(extract_columns(relax_data), index, "R1_sim"), use_list(extract_columns(relax_data), index, "R1_sim"), label='R1 Data_' + rep_name + '_' + ff, marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
					axs[j, i].plot(adjust_existing_res(extract_columns(relax_data), index, "R1_exp"), use_list(extract_columns(relax_data), index, "R1_exp"), label='R1 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
					axs[j, i].set_yticks(np.arange(0, extreme_value(extract_columns(relax_data), max, "R1_sim", "R1_exp")+1, 0.5))
					axs[j, i].set_xticks(np.arange(9, res_nr, 10))
					axs[j, i].set_xlim(0, res_nr-10)
					axs[j, i].set_ylim(0, extreme_value(extract_columns(relax_data), max, "R1_sim", "R1_exp")+1)
					axs[j, i].set_xlabel('Residue number')
					axs[j, i].set_ylabel('R1 rates (1/s)')
					axs[j, i].legend()

	elif len(input) != 1:	
		col=len(input)
		fig, axs = plt.subplots(1, col, figsize=(15, 6))
		for i in range(col):
			index=list(extract_columns(relax_data).keys()).index(Best_cases[i])
			name=Best_cases[i]
			axs[i].plot(adjust_existing_res(extract_columns(relax_data), index, "R1_sim"), use_list(extract_columns(relax_data), index, "R1_sim"), label='R1 Data_' + name, marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
			axs[i].plot(adjust_existing_res(extract_columns(relax_data), index, "R1_exp"), use_list(extract_columns(relax_data), index, "R1_exp"), label='R1 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
			axs[i].set_yticks(np.arange(0, extreme_value(extract_columns(relax_data), max, "R1_sim", "R1_exp")+1, 0.5))
			axs[i].set_xticks(np.arange(9, res_nr, 10))
			axs[i].set_xlim(1, res_nr-10)
			axs[i].set_ylim(0, extreme_value(extract_columns(relax_data), max, "R1_sim", "R1_exp"))
			axs[i].set_xlabel('Residue number')
			axs[i].set_ylabel('R1 rates (1/s)')
			axs[i].legend()

	plt.tight_layout()
	plt.savefig(relax_folder + output +'.png')
	plt.close()

R1_relaxation_combined(list(extract_columns(relax_data).keys()), 'R1_relaxation_combined_plot')

def R2_relaxation_combined(input, output):
	if len(input) > 10:
		fig, axs = plt.subplots(5, 5, figsize=(30, 15))
		for j, ff in enumerate(FORCEFIELDS):
			for item in input:
				if ff in item:
					index=input.index(item) 
					i=int(item.split("/")[0][-1])-1
					rep_name=item.split("/")[0]
					axs[j, i].plot(adjust_existing_res(extract_columns(relax_data), index, "R2_sim"), use_list(extract_columns(relax_data), index, "R2_sim"), label='R2 Data_' + rep_name + '_' + ff, marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
					axs[j, i].plot(adjust_existing_res(extract_columns(relax_data), index, "R2_exp"), use_list(extract_columns(relax_data), index, "R2_exp"), label='R2 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
					axs[j, i].set_yticks(np.arange(0, extreme_value(extract_columns(relax_data), max, "R2_sim", "R2_exp"), 5))
					axs[j, i].set_xticks(np.arange(9, res_nr, 10))
					axs[j, i].set_xlim(1, res_nr-10)
					axs[j, i].set_ylim(0, extreme_value(extract_columns(relax_data), max, "R2_sim", "R2_exp"))
					axs[j, i].set_xlabel('Residue number')
					axs[j, i].set_ylabel('R2 rates (1/s)')
					axs[j, i].legend()
	elif len(input) != 1:	
		col=len(input)
		fig, axs = plt.subplots(1, col, figsize=(15, 6))
		for i in range(col):
			index=list(extract_columns(relax_data).keys()).index(Best_cases[i])
			name=Best_cases[i]
			axs[i].plot(adjust_existing_res(extract_columns(relax_data), index, "R2_sim"), use_list(extract_columns(relax_data), index, "R2_sim"), label='R2 Data_' + name, marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
			axs[i].plot(adjust_existing_res(extract_columns(relax_data), index, "R2_exp"), use_list(extract_columns(relax_data), index, "R2_exp"), label='R2 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
			axs[i].set_yticks(np.arange(0, extreme_value(extract_columns(relax_data), max, "R2_sim", "R2_exp"), 20))
			axs[i].set_xticks(np.arange(9, res_nr, 10))
			axs[i].set_xlim(1, res_nr-10)
			axs[i].set_ylim(0, extreme_value(extract_columns(relax_data), max, "R2_sim", "R2_exp"))
			axs[i].set_xlabel('Residue number')
			axs[i].set_ylabel('R2 rates (1/s)')
			axs[i].legend()
		


	plt.tight_layout()
	plt.savefig(relax_folder + output +'.png')
	plt.close()

R2_relaxation_combined(list(extract_columns(relax_data).keys()), 'R2_relaxation_combined_plot')

def NOE_relaxation_combined(input, output):
	if len(input) > 10:
		fig, axs = plt.subplots(5, 5, figsize=(30, 15))
		for j, ff in enumerate(FORCEFIELDS):
			for item in input:
				if ff in item:
					index=input.index(item) 
					i=int(item.split("/")[0][-1])-1
					rep_name=item.split("/")[0]
					axs[j, i].plot(adjust_existing_res(extract_columns(relax_data), index, "NOE_sim"), use_list(extract_columns(relax_data), index, "NOE_sim"), label='NOE Data_' + rep_name + '_' + ff, marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
					axs[j, i].plot(adjust_existing_res(extract_columns(relax_data), index, "NOE_exp"), use_list(extract_columns(relax_data), index, "NOE_exp"), label='NOE Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
					axs[j, i].set_yticks(np.arange(extreme_value(extract_columns(relax_data), min, "NOE_sim", "NOE_exp"), extreme_value(extract_columns(relax_data), max, "NOE_sim", "NOE_exp"), 0.5))
					axs[j, i].set_xticks(np.arange(9, res_nr, 10))
					axs[j, i].set_xlim(1, res_nr-10)
					axs[j, i].set_ylim(extreme_value(extract_columns(relax_data), min, "NOE_sim", "NOE_exp"), extreme_value(extract_columns(relax_data), max, "NOE_sim", "NOE_exp"))
					axs[j, i].set_xlabel('Residue number')
					axs[j, i].set_ylabel('NOE_values (1/s)')
					axs[j, i].set_title('NOE_data')
					axs[j, i].legend()
	elif len(input) != 1:
		col=len(input)
		fig, axs = plt.subplots(1, col, figsize=(15, 6))
		for i in range(col):
			index=list(extract_columns(relax_data).keys()).index(Best_cases[i])
			name=Best_cases[i]
			axs[i].plot(adjust_existing_res(extract_columns(relax_data), index, "NOE_sim"), use_list(extract_columns(relax_data), index, "NOE_sim"), label='NOE Data_' + name, marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
			axs[i].plot(adjust_existing_res(extract_columns(relax_data), index, "NOE_exp"), use_list(extract_columns(relax_data), index, "NOE_exp"), label='NOE Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
			axs[i].set_yticks(np.arange(0, extreme_value(extract_columns(relax_data), max, "NOE_sim", "NOE_exp"), 0.5))
			axs[i].set_xticks(np.arange(9, res_nr, 10))
			axs[i].set_xlim(1, res_nr-10)
			axs[i].set_ylim(extreme_value(extract_columns(relax_data), min, "NOE_sim", "NOE_exp"), extreme_value(extract_columns(relax_data), max, "NOE_sim", "NOE_exp"))
			axs[i].set_xlabel('Residue number')
			axs[i].set_ylabel('NOE_values (1/s)')
			axs[i].set_title('NOE_data')
			axs[i].legend()
		

	plt.tight_layout()
	plt.savefig(relax_folder + output +'.png')
	plt.close()

NOE_relaxation_combined(list(extract_columns(relax_data).keys()), 'NOE_relaxation_combined_plot')


rog_values_best_all=[]


def plot_rog_density_landscape_all(input, output):
	if len(input) > 10:
		fig, axs = plt.subplots(5, 5, figsize=(30, 15))
		for j, ff in enumerate(FORCEFIELDS):
			for item in input:
				if ff in item:
					rog_values=[]
					rep_name = item.split("/")[0]
					rog_data=glob.glob(SIM_DIR + rep_name + "/" + ff + "/md*gyrate.xvg")[0]
					with open(rog_data, 'r') as file:
						lines = file.readlines()
						for line in range(27, len(lines)):
							parts =lines[line].split()
							try:
								rog_values.append(float(parts[1]))
							except IndexError:
								pass
					rounded_values = [round(value, 1) for value in rog_values]
					rounded_counts = {value: rounded_values.count(value) for value in set(rounded_values)}
					sorted_items = sorted(rounded_counts.items())
					values = [value[0] for value in sorted_items]
					counts = [value[1] for value in sorted_items]

					Rog_data_list.append((rep_name + "/" + ff, rog_values))
	
					i=int(rep_name[-1])-1
#					axs[j, i].set_yticks(np.arange(0, 70000, 10000))
					axs[j, i].set_xticks(np.arange(0.5, 7, 1))
					axs[j, i].set_xlim(0.5, 7)
#					axs[j, i].set_ylim(0, 70000)
	
					axs[j, i].plot(values, counts)
					axs[j, i].axvline(x=statistics.mean(rog_values), color='black', linestyle='--', label='Mean Rog value')
					axs[j, i].set_xlabel('Count')
					axs[j, i].set_ylabel('RoG (nm)')
					axs[j, i].set_title(rep_name + "/" + ff)
	else:
		col=len(input)
		for i in range(col):
			rog_values=[]
			name=input[i]
			ff = name.split("/")[1]
			rep_name = name.split("/")[0]
			rog_data=glob.glob(SIM_DIR + name + "/md*gyrate.xvg")[0]
			with open(rog_data, 'r') as file:
				lines = file.readlines()
				for line in range(27, len(lines)):
					parts = lines[line].split()
					try:
						rog_values.append(float(parts[1]))
						rog_values_best_all.append(float(parts[1]))
					except IndexError:
						pass
			rounded_values = [round(value, 1) for value in rog_values]
			rounded_counts = {value: rounded_values.count(value) for value in set(rounded_values)}
			sorted_items = sorted(rounded_counts.items())
			values = [value[0] for value in sorted_items]
			counts = [value[1] for value in sorted_items]

			Rog_data_list.append((rep_name + "/" + ff, rog_values))
	
			
			plt.plot(values, counts, color=color_list[i], label=rep_name + '_' + ff)
			plt.xlabel('Count')
			plt.ylabel('RoG (nm)')
			plt.axvline(x=statistics.mean(rog_values), color=color_list[i], linestyle='--', label='Mean Rog value')
			plt.legend()
		plt.title("Rog landscape best cases")
	

	plt.tight_layout()
	plt.savefig(relax_folder + output + '.png')
	plt.close()

plot_rog_density_landscape_all(list(extract_columns(relax_data).keys()), 'density_landscape_plot')
plot_rog_density_landscape_all(Best_cases, 'Accepted_cases/density_landscape_plot_best')

def dif_to_exp(input, output):

	used_labels = [[], [], []]
	fig, axs = plt.subplots(1, 3, figsize=(15, 6))
	for i, ff in enumerate(FORCEFIELDS):
		list=[[], [], []]
		for num, item in enumerate(input):
			if ff in item and item in input:
				R1_diff=use_list(extract_columns(relax_data), num, "R1_diff")
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

				R2_diff=use_list(extract_columns(relax_data), num, "R2_diff")
				try:
					if model in used_labels[1]:
						axs[1].scatter(i, statistics.mean(R2_diff), zorder=5, color=color_map.get(model, 'black'))
					else:
						axs[1].scatter(i, statistics.mean(R2_diff), zorder=5, label=model, color=color_map.get(model, 'black'))
						used_labels[1].append(model)
				except:
					pass
				list[1].extend(R2_diff)
			
				NOE_diff=use_list(extract_columns(relax_data), num, "NOE_diff")
				try:
					if model in used_labels[2]:
						axs[2].scatter(i, statistics.mean(NOE_diff), zorder=5, color=color_map.get(model, 'black'))
					else:
						axs[2].scatter(i, statistics.mean(NOE_diff), zorder=5, label=model, color=color_map.get(model, 'black'))
						used_labels[2].append(model)
				except:
					pass
				list[2].extend(NOE_diff)
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
			NOE_avg=statistics.mean(list[2])
			axs[2].bar(i, NOE_avg, label=ff)
		except: 
			pass

	axs[0].set_title('R1 (1/s)')
	axs[0].legend()
	axs[1].set_title('R2 (1/s)')
	axs[1].legend()
	axs[2].set_title('hetNOE_difference')
	axs[2].legend()


	fig.suptitle("Comparison of difference to experimental extract_columns(relax_data)")
	plt.tight_layout()
	plt.savefig(relax_folder + output + '.png')
	plt.close()

dif_to_exp(list(extract_columns(relax_data).keys()), 'Difference_to_experiment_plot')



fig, axs = plt.subplots(1, 3, figsize=(15, 6))
R1_avg=statistics.mean(use_list(avg_data, 0, "R1_diff"))
axs[0].bar(0, R1_avg)
R2_avg=statistics.mean(use_list(avg_data, 0, "R2_diff"))
axs[1].bar(0, R2_avg)
NOE_avg=statistics.mean(use_list(avg_data, 0, "NOE_diff"))
axs[2].bar(0, NOE_avg)

axs[0].set_title('R1 (1/s)')
axs[1].set_title('R2 (1/s)')
axs[2].set_title('hetNOE_difference')

	
fig.suptitle("Difference to experimental extract_columns(relax_data) for average structure")
plt.tight_layout()
plt.savefig(relax_folder + 'Accepted_cases/Difference_to_experiment_plot_avg.png')
plt.close()


timescale_data=sorted(glob.glob(SIM_DIR + 'model*/*/Ctimes_Coeffs.csv'))
timescale_data_avg=glob.glob(Unst_folder + "Ctimes_Coeffs.csv")

def create_timescale_scatter_plot(path, idx, axs=None):
	C_times_ns=extract_values_pandas(path, idx, 0)
	Coeffs=extract_values_pandas(path, idx, 1)

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
		axs.set_yticks(np.arange(0, 81, 20)) 
		axs.set_ylim(0, 80)
		axs.set_title('/'.join(path[idx].split("/")[-3:-1]))


create_timescale_scatter_plot(timescale_data_avg, 0, axs=None)
plt.savefig(best_cases_folder + 'Timescale_plot_avg.png')
plt.close()

plot_all(timescale_data, "Timescale_plot_all.png", create_timescale_scatter_plot)



def tau_eff_area(path, idx, axs=None):
	Tau_values=extract_columns(path)
	Tau_values_adjust=[x*10**10 for x in use_list(Tau_values, idx, "Tau_eff_area")]
	if axs is None:
		plt.plot(adjust_existing_res(Tau_values, idx, "Tau_eff_area"), Tau_values_adjust, label='Average effective correlation time', marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
		plt.xlabel('Residue number')
		plt.ylabel('Effective correlation time ($10^{10}$ s)')
	else:
		axs.plot(adjust_existing_res(Tau_values, idx, "Tau_eff_area"), Tau_values_adjust, label='/'.join(path[idx].split("/")[-3:-1]), marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
		#axs.set_yticks(np.arange(0, extreme_value(Tau_values, max, "Tau_eff_area", None)*10**10, 0.5))
		axs.set_xticks(np.arange(9, res_nr, 10))
		axs.set_xlim(1, res_nr-10)
		#axs.set_ylim(0, extreme_value(Tau_values, max, "Tau_eff_area", None)*10**10)
		axs.set_xlabel('Residue number')
		axs.set_ylabel('Effective correlation time ($10^{10}$ s)')
		axs.legend()


tau_eff_area(avg_data_path, 0, axs=None)
plt.savefig(best_cases_folder + "Tau_effective_area_avg.png")
plt.close()

plot_all(relax_data, "Tau_effective_area_all.png", tau_eff_area)



def plot_avg_rog_bar(input, output):
	used_labels_rog = []
	for i, ff in enumerate(FORCEFIELDS):
		list=[]
		for item in input:
			if ff in item[0]:
				model=item[0].replace("/" + ff, "")
				try:
					if model in used_labels_rog:
						plt.scatter(i, statistics.mean(item[1]), zorder=5, color=color_map.get(model, 'black'))
					else:
						plt.scatter(i, statistics.mean(item[1]), zorder=5, label=model, color=color_map.get(model, 'black'))
						used_labels_rog.append(model)
				except:
					pass
				list.extend(item[1])
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

plot_avg_rog_bar(Rog_data_list, "Average_rog_plot")


def plot_ensembles_images(input, output):
	num_images = len(input)
	if num_images > 10:
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
	elif num_images == 1:
		ensemble = input[0]
		rep_name = ensemble.split("/")[-3]
		forcefield = ensemble.split("/")[-2]
		img = imread(ensemble)
		plt.imshow(img)
		plt.title(f'{rep_name}/{forcefield}')
		plt.axis('off')
	else:
		col = num_images
		fig, axs = plt.subplots(1, col, figsize=(15, 6)) 
		for i in range(col): 
			ensemble = input[i] 
			rep_name = ensemble.split("/")[-3] 
			forcefield = ensemble.split("/")[-2] 
			img = imread(ensemble) 
			axs[i].imshow(img) 
			axs[i].set_title(f'{rep_name}/{forcefield}') 
			axs[i].axis('off')

	plt.tight_layout()
	plt.savefig(f'{relax_folder}/{output}.png')
	plt.close()

plot_ensembles_images(sorted(glob.glob(SIM_DIR + "model*/*/*mdmat*.png")), "Contact_map_combined")
plot_ensembles_images([glob.glob(SIM_DIR + i + "/md*mdmat.png") for i in Best_cases][0], "Accepted_cases/Contact_map_combined_best")

plot_ensembles_images(sorted(glob.glob(SIM_DIR + "model*/*/*correlation*.png")), "Correlation_combined")
plot_ensembles_images([glob.glob(SIM_DIR + i + "/*correlation*.png") for i in Best_cases][0], "Accepted_cases/Correlation_combined_best")
'''
plot_ensembles_images(sorted(glob.glob(SIM_DIR + "model*/*/Ensemble*model*aligned_fig.png")), "Ensembles_aligned_combined")
plot_ensembles_images([glob.glob(SIM_DIR + i + "/Ensemble*model*aligned_fig.png") for i in Best_cases][0], "Accepted_cases/Ensembles_aligned_best")
'''


def contact_map_avg(input, output):
	matrix = avg_csv(input)
	fig, ax = plt.subplots()
	cmap = plt.cm.viridis
	plot = ax.imshow(matrix, cmap=cmap, origin="lower")
	
	ax.set_xlabel("Residue",fontsize=15)
	ax.set_ylabel("Residue",fontsize=15)

	plt.xticks(np.arange(1, res_nr+1, 10), fontsize=15)
	plt.yticks(np.arange(1, res_nr+1, 10), fontsize=15)

	cbar = plt.colorbar(plot, ax=ax)
	cbar.set_label("Distance (nm)",fontsize=15)
	plt.savefig(output + "Avg_contact.png", dpi=600)
	plt.close()

contact_map_avg([glob.glob(SIM_DIR + i + "/md*mdmat.csv") for i in Best_cases][0], best_cases_folder)

'''
def corr_map_avg(input, output):
	matrix = avg_csv(input)

	color_map = plt.imshow(matrix,vmin=-1, vmax=1, origin='lower')
	color_map.set_cmap("seismic")
	plt.colorbar()
	plt.savefig(output + "Avg_contact.png", dpi=600)
	plt.close()

contact_map_avg([glob.glob(SIM_DIR + i + "/*correlation*.png") for i in Best_cases][0], best_cases_folder)



rounded_values = [round(value, 1) for value in rog_values_best_all]
rounded_counts = {value: rounded_values.count(value) for value in set(rounded_values)}
sorted_items = sorted(rounded_counts.items())
values = [value[0] for value in sorted_items]
counts = [value[1] for value in sorted_items]


with open(Unst_folder + "Best_rog_landscape.txt", 'w') as f:
	f.write(f"Rog_avg:\t{statistics.mean(rog_values_best_all):.5f}\n")
	#f.write(f"{statistics.mean(rog_values_best_all):.3f}\n")
	for x, y in zip(values, counts):
		f.write(f"{x:.3f}\t{y:.5f}\n")


pdb_data = sorted(glob.glob(SIM_DIR + "/model*/*/"))
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

cmd.png(best_cases_folder + 'Ensemble_' + SIM_DIR.split('/')[-1] + '_aligned_fig.png')	
cmd.delete('all')
cmd.quit()

subprocess.run(["python3", avg_script])
subprocess.run(["python3", pymol_analysis])






'''
