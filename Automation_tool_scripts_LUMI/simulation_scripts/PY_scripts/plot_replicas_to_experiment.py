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


#SIM_DIR='/scratch/project_462000285/cmcajsa/systems/forcefield_compare/Unst_alphasynuclein/'
SIM_DIR = os.getcwd()
BASE_DIR=os.path.dirname(SIM_DIR)
exp_data = BASE_DIR + '/' + SIM_DIR.split("/")[-1] + '_exp_data.txt'
Unst_folder = SIM_DIR.replace(SIM_DIR.split("/")[-1], '') + "results/" + SIM_DIR.split("/")[-1] + '/'
relax_folder=SIM_DIR.replace(SIM_DIR.split("/")[-1], '') + "results/" + SIM_DIR.split("/")[-1] + '/rep_to_exp_data/'
avg_script=BASE_DIR + "/simulation_scripts/PY_scripts/average_contact.py"
py_script=BASE_DIR + "/simulation_scripts/PY_scripts/Old_Relaxations_for_Samuli.py"


best_cases_folder = os.path.join(relax_folder, "Accepted_cases/")

os.makedirs(relax_folder, exist_ok=True)
os.makedirs(best_cases_folder, exist_ok=True)

with open(exp_data, 'r') as file:
	lines = file.readlines()
	first_row = lines[0].split()
	magn_field=first_row[5]
	res_nr=len(lines)-1


R1_values_exp = [[], []]
R2_values_exp = [[], []]
NOE_values_exp = [[], []]

# Read data from the text file
with open(exp_data, 'r') as file:
	lines = file.readlines()
	for line in range(1, len(lines)):
                # Split the line into T1, T2, and NOE values
		parts = lines[line].split()
		try:
			R1_values_exp[1].append(float(parts[1]))
			R1_values_exp[0].append(line)
		except ValueError:
			R1_values_exp[1].append("n")
			R1_values_exp[0].append("n")
		try:
			R2_values_exp[1].append(float(parts[3]))
			R2_values_exp[0].append(line)
		except ValueError:
                        R2_values_exp[1].append("n")
                        R2_values_exp[0].append("n")
		try:
			NOE_values_exp[1].append(float(parts[5]))
			NOE_values_exp[0].append(line)
		except ValueError:
			NOE_values_exp[1].append("n")
			NOE_values_exp[0].append("n")

def filter_list(data):
        return [[x for x, y in zip(data[0], data[1]) if y != 'n'],
                [y for y in data[1] if y != 'n']]

R1_values_exp_filtered = filter_list(R1_values_exp)
R2_values_exp_filtered = filter_list(R2_values_exp)
NOE_values_exp_filtered = filter_list(NOE_values_exp)



combined_list_of_relax_times=[]

R1_diff=[]
R2_diff=[]
NOE_diff=[]

for item in FORCEFIELDS:
	relax_data=glob.glob(SIM_DIR+"/model*/"+item+'/relaxation_times.txt')
	for data in sorted(relax_data):
		rep_name = data.split("/")[-3]
		combined_list=[["R1"], ["R2"], ["NOE"]]
		existing_res=[[], [], []]
		R1_values_diff=[]
		R2_values_diff=[]
		NOE_values_diff=[]
		with open(data, 'r') as file:
			lines = file.readlines()
			for line in range(0, len(lines)):
				parts = lines[line].split()
				try:
					combined_list[0].append(1 / float(parts[1]))
					existing_res[0].append(line+1)
					if R1_values_exp[1][line] != "n":
						R1_values_diff.append((1 / float(parts[1])-R1_values_exp[1][line]))
				except:
					pass
				try:
					combined_list[1].append(1 / float(parts[3]))
					existing_res[1].append(line+1)
					if R2_values_exp[1][line] != "n":
						R2_values_diff.append((1 / float(parts[3])-R2_values_exp[1][line]))
				except:	
					pass
				try:
					combined_list[2].append(float(parts[5]))
					existing_res[2].append(line+1)
					if NOE_values_exp[1][line] != "n":
						NOE_values_diff.append((float(parts[5])-NOE_values_exp[1][line]))
				except:
					pass
		if len(R1_values_diff) != 0 and sum(R1_values_diff) != 0:
			RMSD_R1=f"R1: {math.sqrt(sum(x**2 for x in R1_values_diff) / len(R1_values_diff))}"
			R1_diff.append((rep_name + "/" + item, RMSD_R1, R1_values_diff))
		if len(R2_values_diff) != 0 and sum(R2_values_diff) != 0:
			RMSD_R2=f"R2: {math.sqrt(sum(x**2 for x in R2_values_diff) / len(R2_values_diff))}"
			R2_diff.append((rep_name + "/" + item, RMSD_R2, R2_values_diff))
		if len(NOE_values_diff) != 0 and sum(NOE_values_diff) != 0:
			RMSD_NOE=f"NOE: {math.sqrt(sum(x**2 for x in NOE_values_diff) / len(NOE_values_diff))}"
			NOE_diff.append((rep_name + "/" + item, RMSD_NOE, NOE_values_diff))
		if len(combined_list[0])>1:
			combined_list_of_relax_times.append((rep_name + "/" + item, combined_list, existing_res))





Names=[item[0] for item in combined_list_of_relax_times]
Value_lists=[item[1] for item in combined_list_of_relax_times]
Existing_res=[item[2] for item in combined_list_of_relax_times]


R1_lists=[item[0] for item in Value_lists]
R2_lists=[item[1] for item in Value_lists]
NOE_lists=[item[2] for item in Value_lists]

R1_existing=[item[0] for item in Existing_res]
R2_existing=[item[1] for item in Existing_res]
NOE_existing=[item[2] for item in Existing_res]

R1_RMSD_values = [item[1].replace('R1: ', '') for item in R1_diff]
R2_RMSD_values = [item[1].replace('R2: ', '') for item in R2_diff]
NOE_RMSD_values = [item[1].replace('NOE: ', '') for item in NOE_diff]

R1_RMSD_name = [item[0] for item in R1_diff]
R2_RMSD_name = [item[0] for item in R2_diff]
NOE_RMSD_name = [item[0] for item in NOE_diff]



Rog_data_list=[]

rog_values_best=[]

Best_cases=[]

data=relax_folder + 'ranking_table.csv'
with open(data, 'w', newline="") as csvfile:
	csvwriter = csv.writer(csvfile)
	csvwriter.writerow(["Force field", "Replica", "R1 RMSD", "R1 (%)", "R2 RMSD", "R2 (%)", "hetNOE RMSD", "hetNOE (%)", "Sum (%)"])
	for i in range(0, 25):
		try:
			rep_name=R1_RMSD_name[i].split("/")[0]
			ff_name=R1_RMSD_name[i].split("/")[1]
			R1_value=float(R1_RMSD_values[i])
			R2_value=float(R2_RMSD_values[i])
			NOE_value=float(NOE_RMSD_values[i])
			R1_rank_value=float(R1_RMSD_values[i])/min(float(i) for i in R1_RMSD_values)*100
			R2_rank_value=float(R2_RMSD_values[i])/min(float(i) for i in R2_RMSD_values)*100
			NOE_rank_value=float(NOE_RMSD_values[i])/min(float(i) for i in NOE_RMSD_values)*100
			Ranking_sum=R1_rank_value+R2_rank_value+NOE_rank_value
			csvwriter.writerow([ff_name, rep_name, R1_value, R1_rank_value, R2_value, R2_rank_value, NOE_value, NOE_rank_value, Ranking_sum])

			if R1_rank_value/100 < 1.5 and R2_rank_value/100 < 1.5 and NOE_rank_value/100 < 1.5:
				print(R1_diff[i][0], R1_diff[i][1])
				Best_cases.append((R1_diff[i][0], R1_diff[i][1]))
				with open(SIM_DIR+"/"+ R1_diff[i][0] +'/md_2000ns_gyrate.xvg', 'r') as file:
					lines = file.readlines()
					for line in range(27, len(lines)):
						parts = lines[line].split()
						rog_values_best.append(round(float(parts[1]), 1))
		except:
			pass


value_counts = Counter(rog_values_best)
values, counts = zip(*sorted(value_counts.items(), key=lambda x: x[0]))
plt.plot(values, counts)
plt.xlabel('Radius of gyration')
plt.ylabel('Count')
plt.title('Count of Radius of gyration values from best simulation data')
plt.tight_layout()	
plt.savefig(best_cases_folder + 'best_rog_density_landscape_plot.png')
plt.close()

Best_cases_names=[item[0] for item in Best_cases]

os.makedirs(Unst_folder + "correlation_functions/", exist_ok=True)

i = 1
while i <= res_nr :
	list_of_lists=[]
	for j in Best_cases_names:
		try:
			list=[]
			file_name = SIM_DIR + "/" + j + "/correlation_functions/NHrotaCF_"+ str(i) + ".xvg"
			print(file_name)
			with open(file_name, 'r') as file:
				lines = file.readlines()
				for line in range(0, len(lines)-1):
					parts = lines[line].split()
					list.append(float(parts[1]))
			list_of_lists.append(list)
		except:
			pass
	try:
		average_list=[sum(x) / len(x) for x in zip(*list_of_lists)]
		times = [i * 10.000 for i in range(0, len(list_of_lists[0]))]
		new_file=Unst_folder + 'correlation_functions/NHrotaCF_' + str(i) + '.xvg'
		with open(new_file, 'w') as f:
			for x, y in zip(times, average_list):
				f.write(f"{x:.3f}\t{y:.5f}\n")
	
	except:
		pass
	list_of_lists=[]			
	i += 1

shutil.copy(py_script, Unst_folder)

with fileinput.FileInput(Unst_folder + "Old_Relaxations_for_Samuli.py", inplace=True) as file:
	for line in file:
		line = line.replace('magn_field=magn_field', 'magn_field=' + str(magn_field))
		print(line, end="")

os.chdir(Unst_folder)
subprocess.run(["python3", Unst_folder + "Old_Relaxations_for_Samuli.py"], stdout=open(Unst_folder + "relaxation_times.txt", 'w'))
os.chdir(SIM_DIR)

combined_list_average=[[], [], []]
existing_res_avg=[[], [], []]
R1_values_diff_avg=[]
R2_values_diff_avg=[]
NOE_values_diff_avg=[]
with open(Unst_folder + "relaxation_times.txt", "r") as file:
	lines = file.readlines()
	for line in range(0, len(lines)):
		parts = lines[line].split()
		try:
			combined_list_average[0].append(1 / float(parts[1]))
			existing_res_avg[0].append(line+1)
			if R1_values_exp[1][line] != "n":
				R1_values_diff_avg.append((1 / float(parts[1])-R1_values_exp[1][line]))
		except:
			pass
		try:
			combined_list_average[1].append(1 / float(parts[3]))
			existing_res_avg[1].append(line+1)
			if R2_values_exp[1][line] != "n":
				R2_values_diff_avg.append((1 / float(parts[1])-R2_values_exp[1][line]))
		except:
			pass
		try:
			combined_list_average[2].append(float(parts[5]))
			existing_res_avg[2].append(line+1)
			if NOE_values_exp[1][line] != "n":
				NOE_values_diff_avg.append((float(parts[5])-NOE_values_exp[1][line]))
		except:
			pass


i=0

for item in FORCEFIELDS:
	fig, axs = plt.subplots(1, 3, figsize=(15, 6))
	for i, name in enumerate(Names):
		if name.split("/")[1]==item:
			rep_name=name.split("/")[0]
			selected = color_map.get(rep_name, 'black')
			list=Value_lists[i]
			axs[0].plot(R1_existing[i], list[0][1:], label='R1 Data_' + rep_name, marker='o', linestyle='-', lw=1.0, markersize=2, color=selected)
			axs[0].set_xlabel('Residue number')
			axs[0].set_ylabel('R1_values (1/s)')
			axs[0].set_title('R1_data')
			axs[0].legend()
			axs[1].plot(R2_existing[i], list[1][1:], label='R2 Data_' + rep_name, marker='o', linestyle='-', lw=1.0, markersize=2, color=selected)
			axs[1].set_xlabel('Residue number')
			axs[1].set_ylabel('R2_values (1/s)')
			axs[1].set_title('R2_data')
			axs[1].legend()
			axs[2].plot(NOE_existing[i], list[2][1:], label='NOE Data_' + rep_name, marker='o', linestyle='-', lw=1.0, markersize=2, color=selected)
			axs[2].set_xlabel('Residue number')				
			axs[2].set_ylabel('NOE_values 1/s')
			axs[2].set_title('NOE_data')
			axs[2].legend()
	else:
		axs[0].plot(R1_values_exp_filtered[0], R1_values_exp_filtered[1], label='R1 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
		axs[0].set_xlabel('Residue number')
		axs[0].set_ylabel('R1_values (1/s)')
		axs[0].set_title('R1_data')
		axs[0].legend()
		
		axs[1].plot(R2_values_exp_filtered[0], R2_values_exp_filtered[1], label='R2 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
		axs[1].set_xlabel('Residue number')
		axs[1].set_ylabel('R2_values (1/s)')
		axs[1].set_title('R2_data')
		axs[1].legend()
	
		axs[2].plot(NOE_values_exp_filtered[0], NOE_values_exp_filtered[1], label='NOE Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
		axs[2].set_xlabel('Residue number')
		axs[2].set_ylabel('NOE_values (1/s)')
		axs[2].set_title('NOE_data')
		axs[2].legend()
	
	plt.tight_layout()
	plt.savefig(relax_folder + item + '_plot.png')
	plt.close('all')


fig, axs = plt.subplots(1, 3, figsize=(15, 6))
axs[0].plot(existing_res_avg[0], combined_list_average[0], label='R1 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color="blue")
axs[0].set_xlabel('Residue number')
axs[0].set_ylabel('R1_values (1/s)')
axs[0].set_title('R1_data')
axs[0].legend()
axs[1].plot(existing_res_avg[1], combined_list_average[1], label='R2 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color="blue")
axs[1].set_xlabel('Residue number')
axs[1].set_ylabel('R2_values (1/s)')
axs[1].set_title('R2_data')
axs[1].legend()
axs[2].plot(existing_res_avg[2], combined_list_average[2], label='NOE Data', marker='o', linestyle='-', lw=1.0, markersize=2, color="blue")
axs[2].set_xlabel('Residue number')
axs[2].set_ylabel('NOE_values (1/s)')
axs[2].set_title('NOE_data')
axs[2].legend()

axs[0].plot(R1_values_exp_filtered[0], R1_values_exp_filtered[1], label='R1 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
axs[0].set_xlabel('Residue number')
axs[0].set_ylabel('R1_values (1/s)')
axs[0].set_title('R1_data')
axs[0].legend()
axs[1].plot(R2_values_exp_filtered[0], R2_values_exp_filtered[1], label='R2 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
axs[1].set_xlabel('Residue number')
axs[1].set_ylabel('R2_values (1/s)')
axs[1].set_title('R2_data')
axs[1].legend()
axs[2].plot(NOE_values_exp_filtered[0], NOE_values_exp_filtered[1], label='NOE Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
axs[2].set_xlabel('Residue number')
axs[2].set_ylabel('NOE_values (1/s)')
axs[2].set_title('NOE_data')
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
					axs[j, i].plot(R1_existing[index], R1_lists[index][1:], label='R1 Data_' + rep_name + '_' + ff, marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
					axs[j, i].plot(R1_values_exp_filtered[0], R1_values_exp_filtered[1], label='R1 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
					#axs[j, i].set_yticks(np.arange(0, int(max(max(sublist[1:]) for sublist in R1_lists))+0.5, 0.5))
					axs[j, i].set_yticks(np.arange(0, 2.5, 0.5))
					axs[j, i].set_xticks(np.arange(1, len(R1_lists[index]), 20))
					axs[j, i].set_xlim(1, len(R1_lists[index]))
					#axs[j, i].set_ylim(0, int(max(max(sublist[1:]) for sublist in R1_lists))+0.5)
					axs[j, i].set_ylim(0, 2.5)
					axs[j, i].set_xlabel('Residue number')
					axs[j, i].set_ylabel('R1_values (1/s)')
					axs[j, i].set_title('R1_data')
					axs[j, i].legend()
	elif len(input) != 1:	
		col=len(input)
		fig, axs = plt.subplots(1, col, figsize=(15, 6))
		for i in range(col):
			print(i)
			data = Names.index(Best_cases_names[i])
			axs[i].plot(R1_existing[data], R1_lists[data][1:], label='R1 Data_' + Names[data], marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
			axs[i].plot(R1_values_exp_filtered[0], R1_values_exp_filtered[1], label='R1 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
			#axs[i].set_yticks(np.arange(0, int(max(max(sublist[1:]) for sublist in R1_lists))+0.5, 0.5))
			axs[i].set_xticks(np.arange(1, len(R1_lists[i]), 20))
			axs[i].set_xlim(1, len(R1_lists[data]))
			axs[i].set_ylim(0, int(max(max(sublist[1:]) for sublist in R1_lists))+0.5)
			axs[i].set_xlabel('Residue number')
			axs[i].set_ylabel('R1_values (1/s)')
			axs[i].set_title('R1_data')
			axs[i].legend()
	else:
		data = Names.index(Best_cases_names[0])
		plt.plot(R1_existing[data], R1_lists[data][1:], label='R1 Data_' + Names[data], marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
		plt.plot(R1_values_exp_filtered[0], R1_values_exp_filtered[1], label='R1 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
		plt.xlabel('Residue number')
		plt.ylabel('R1_values (1/s)')
		plt.title('R1_data')
		plt.legend()

	plt.tight_layout()
	plt.savefig(relax_folder + output +'.png')
	plt.close()

R1_relaxation_combined(Names, 'R1_relaxation_combined_plot')
R1_relaxation_combined(Best_cases_names, 'Accepted_cases/R1_relaxation_combined_plot_best')

def R2_relaxation_combined(input, output):
	if len(input) > 10:
		fig, axs = plt.subplots(5, 5, figsize=(30, 15))
		for j, ff in enumerate(FORCEFIELDS):
			for item in input:
				if ff in item:
					index=input.index(item) 
					i=int(item.split("/")[0][-1])-1
					rep_name=item.split("/")[0]
					axs[j, i].plot(R2_existing[index], R2_lists[index][1:], label='R2 Data_' + rep_name + '_' + ff, marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
					axs[j, i].plot(R2_values_exp_filtered[0], R2_values_exp_filtered[1], label='R2 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
					#axs[j, i].set_yticks(np.arange(0, int(max(max(sublist[1:]) for sublist in R2_lists))+0.5, 0.5))
					axs[j, i].set_yticks(np.arange(0, 15, 5))
					axs[j, i].set_xticks(np.arange(1, len(R2_lists[index]), 20))
					axs[j, i].set_xlim(1, len(R2_lists[index]))
					#axs[j, i].set_ylim(0, int(max(max(sublist[1:]) for sublist in R2_lists))+0.5)
					axs[j, i].set_ylim(0, 15)
					axs[j, i].set_xlabel('Residue number')
					axs[j, i].set_ylabel('R2_values (1/s)')
					axs[j, i].set_title('R2_data')
					axs[j, i].legend()
	elif len(input) != 1:	
		col=len(input)
		fig, axs = plt.subplots(1, col, figsize=(15, 6))
		for i in range(col):
			data = Names.index(Best_cases_names[i])
			axs[i].plot(R2_existing[data], R2_lists[data][1:], label='R2 Data_' + Names[data], marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
			axs[i].plot(R2_values_exp_filtered[0], R2_values_exp_filtered[1], label='R2 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
			#axs[i].set_yticks(np.arange(0, int(max(max(sublist[1:]) for sublist in R2_lists))+0.5, 0.5))
			axs[i].set_xticks(np.arange(1, len(R2_lists[i]), 20))
			axs[i].set_xlim(1, len(R2_lists[data]))
			axs[i].set_ylim(0, int(max(max(sublist[1:]) for sublist in R2_lists))+0.5)
			axs[i].set_xlabel('Residue number')
			axs[i].set_ylabel('R2_values (1/s)')
			axs[i].set_title('R2_data')
			axs[i].legend()
		
	else:
		data = Names.index(Best_cases_names[0])
		plt.plot(R1_existing[data], R1_lists[data][1:], label='R1 Data_' + Names[data], marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
		plt.plot(R1_values_exp_filtered[0], R1_values_exp_filtered[1], label='R1 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
		plt.xlabel('Residue number')
		plt.ylabel('R1_values (1/s)')
		plt.title('R1_data')
		plt.legend()

	plt.tight_layout()
	plt.savefig(relax_folder + output +'.png')
	plt.close()

R2_relaxation_combined(Names, 'R2_relaxation_combined_plot')
R2_relaxation_combined(Best_cases_names, 'Accepted_cases/R2_relaxation_combined_plot_best')

def NOE_relaxation_combined(input, output):
	if len(input) > 10:
		fig, axs = plt.subplots(5, 5, figsize=(30, 15))
		for j, ff in enumerate(FORCEFIELDS):
			for item in input:
				if ff in item:
					index=input.index(item) 
					i=int(item.split("/")[0][-1])-1
					rep_name=item.split("/")[0]
					axs[j, i].plot(NOE_existing[index], NOE_lists[index][1:], label='NOE Data_' + rep_name + '_' + ff, marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
					axs[j, i].plot(NOE_values_exp_filtered[0], NOE_values_exp_filtered[1], label='NOE Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
					#axs[j, i].set_yticks(np.arange(0, int(max(max(sublist[1:]) for sublist in NOE_lists))+1, 1))
					axs[j, i].set_yticks(np.arange(0, 1, 0.2))
					axs[j, i].set_xticks(np.arange(1, len(NOE_lists[index]), 20))
					axs[j, i].set_xlim(1, len(NOE_lists[index]))
					axs[j, i].set_ylim(0, 1)
					#axs[j, i].set_ylim(0, int(max(max(sublist[1:]) for sublist in NOE_lists))+1)
					axs[j, i].set_xlabel('Residue number')
					axs[j, i].set_ylabel('NOE_values (1/s)')
					axs[j, i].set_title('NOE_data')
					axs[j, i].legend()
	elif len(input) != 1:
		col=len(input)
		fig, axs = plt.subplots(1, col, figsize=(15, 6))
		for i in range(col):
			data = Names.index(Best_cases_names[i])
			axs[i].plot(NOE_existing[data], NOE_lists[data][1:], label='NOE Data_' + Names[data], marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
			axs[i].plot(NOE_values_exp_filtered[0], NOE_values_exp_filtered[1], label='NOE Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
			#axs[i].set_yticks(np.arange(0, int(max(max(sublist[1:]) for sublist in NOE_lists))+0.5, 0.5))
			axs[i].set_xticks(np.arange(1, len(NOE_lists[i]), 20))
			axs[i].set_xlim(1, len(NOE_lists[data]))
			axs[i].set_ylim(0, int(max(max(sublist[1:]) for sublist in NOE_lists))+0.5)
			axs[i].set_xlabel('Residue number')
			axs[i].set_ylabel('NOE_values (1/s)')
			axs[i].set_title('NOE_data')
			axs[i].legend()
		
	else:
		data = Names.index(Best_cases_names[0])
		plt.plot(R1_existing[data], R1_lists[data][1:], label='R1 Data_' + Names[data], marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
		plt.plot(R1_values_exp_filtered[0], R1_values_exp_filtered[1], label='R1 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
		plt.xlabel('Residue number')
		plt.ylabel('R1_values (1/s)')
		plt.title('R1_data')
		plt.legend()

	plt.tight_layout()
	plt.savefig(relax_folder + output +'.png')
	plt.close()

NOE_relaxation_combined(Names, 'NOE_relaxation_combined_plot')
NOE_relaxation_combined(Best_cases_names, 'Accepted_cases/NOE_relaxation_combined_plot_best')


rog_values_best_all=[]


def plot_rog_density_landscape_all(input, output):
	if len(input) > 10:
		fig, axs = plt.subplots(5, 5, figsize=(30, 15))
		for j, ff in enumerate(FORCEFIELDS):
			for item in input:
				if ff in item:
					rog_values=[]
					rep_name = item.split("/")[0]
					data=SIM_DIR+"/" + rep_name + "/" + ff + "/md_2000ns_gyrate.xvg"
					with open(data, 'r') as file:
						lines = file.readlines()
						for line in range(27, len(lines)):
							parts = lines[line].split()
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
			data=SIM_DIR+"/" + name + "/md_2000ns_gyrate.xvg"
			with open(data, 'r') as file:
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
	
#			axs[i].set_yticks(np.arange(0, 70000, 10000))
#			axs[i].set_xticks(np.arange(0.5, 7, 1))
#			axs[i].set_xlim(0.5, 7)
#			axs[i].set_ylim(0, 70000)
			
			plt.plot(values, counts, color=color_list[i], label=rep_name + '_' + ff)
			plt.xlabel('Count')
			plt.ylabel('RoG (nm)')
			plt.axvline(x=statistics.mean(rog_values), color=color_list[i], linestyle='--', label='Mean Rog value')
			plt.legend()
		plt.title("Rog landscape best cases")
	

	plt.tight_layout()
	plt.savefig(relax_folder + output + '.png')
	plt.close()

plot_rog_density_landscape_all(Names, 'density_landscape_plot')
plot_rog_density_landscape_all(Best_cases_names, 'Accepted_cases/density_landscape_plot_best')

def dif_to_exp(input, output):

	used_labels = [[], [], []]
	fig, axs = plt.subplots(1, 3, figsize=(15, 6))
	for i, ff in enumerate(FORCEFIELDS):
		list=[[], [], []]
		for item in R1_diff:
			if ff in item[0] and item[0] in input:
				model=item[0].replace("/" + ff, "")
				try:
					if model in used_labels[0]:
						axs[0].scatter(i, statistics.mean(item[2]), zorder=5, color=color_map.get(model, 'black'))
					else:
						axs[0].scatter(i, statistics.mean(item[2]), zorder=5, label=model, color=color_map.get(model, 'black'))
						used_labels[0].append(model)
				except:
					pass
				list[0].extend(item[2])	
		for item in R2_diff:
			if ff in item[0] and item[0] in input:
				model=item[0].replace("/" + ff, "")
				try:
					if model in used_labels[1]:
						axs[1].scatter(i, statistics.mean(item[2]), zorder=5, color=color_map.get(model, 'black'))
					else:
						axs[1].scatter(i, statistics.mean(item[2]), zorder=5, label=model, color=color_map.get(model, 'black'))
						used_labels[1].append(model)
				except:
					pass
				list[1].extend(item[2])
		for item in NOE_diff:
			if ff in item[0] and item[0] in input:
				model=item[0].replace("/" + ff, "")
				try:
					if model in used_labels[2]:
						axs[2].scatter(i, statistics.mean(item[2]), zorder=5, color=color_map.get(model, 'black'))
					else:
						axs[2].scatter(i, statistics.mean(item[2]), zorder=5, label=model, color=color_map.get(model, 'black'))
						used_labels[2].append(model)
				except:
					pass
				list[2].extend(item[2])
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


	fig.suptitle("Comparison of difference to experimental data")
	plt.tight_layout()
	plt.savefig(relax_folder + output + '.png')
	plt.close()

dif_to_exp(Names, 'Difference_to_experiment_plot')


fig, axs = plt.subplots(1, 3, figsize=(15, 6))
R1_avg=statistics.mean(R1_values_diff_avg)
axs[0].bar(0, R1_avg)
R2_avg=statistics.mean(R2_values_diff_avg)
axs[1].bar(0, R2_avg)
NOE_avg=statistics.mean(NOE_values_diff_avg)
axs[2].bar(0, NOE_avg)

axs[0].set_title('R1 (1/s)')
axs[1].set_title('R2 (1/s)')
axs[2].set_title('hetNOE_difference')

	
fig.suptitle("Difference to experimental data for average structure")
plt.tight_layout()
plt.savefig(relax_folder + 'Accepted_cases/Difference_to_experiment_plot_avg.png')
plt.close()



def create_timescale_scatter_plot(data, axs=None):
	lines = open(data, 'r').readlines()
	rep_name = data.split("/")[-3]
	forcefield = data.split("/")[-2]
    
	x_vals, y_vals, weights = [], [], []
    
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
		axs.set_title(rep_name + "/" + forcefield)


coeff_data = sorted(glob.glob(SIM_DIR + "/model*/*/Ctimes_Coeffs.txt"))

for data in coeff_data:
	create_timescale_scatter_plot(data)
plt.tight_layout()
plt.savefig(relax_folder + 'Timescale_plot.png')
plt.close()


timescale_data_best=[]
for j in Best_cases_names:
        case_search=glob.glob(SIM_DIR + "/" + j + "/Ctimes_Coeffs.txt")
        timescale_data_best.append(case_search[0])
for data in timescale_data_best:
        create_timescale_scatter_plot(data)
plt.tight_layout()
plt.savefig(best_cases_folder + 'Timescale_plot_best.png')
plt.close()

fig, axs = plt.subplots(5, 5, figsize=(15, 15))
for i in range(5):
	for j in range(5):
		if 5*i+j < len(coeff_data):
			data = coeff_data[5*i+j]
			create_timescale_scatter_plot(data, axs[j, i])

plt.tight_layout()
plt.savefig(relax_folder + 'Timescale_plot_all.png')
plt.close()

fig, axs = plt.subplots(5, 5, figsize=(15, 15))
for i in range(5):
	for j in range(5):
		if 5*i+j < len(timescale_data_best):
			data = timescale_data_best[5*i+j]
			create_timescale_scatter_plot(data, axs[j, i])

plt.tight_layout()
plt.savefig(best_cases_folder + 'Timescale_plot_all_best.png')
plt.close()

def tau_eff_area(input, output):
	if len(input) > 10:
		fig, axs = plt.subplots(5, 5, figsize=(30, 15))
		for j, ff in enumerate(FORCEFIELDS):
			for item in input:
				if ff in item:
					rep_name=item.split("/")[0]
					data=SIM_DIR + '/' + rep_name + '/' + ff + '/relaxation_times.txt' 
					i=int(rep_name[-1])-1
					rep_name=item.split("/")[0]
					with open(data, 'r') as file:
						lines = file.readlines()
						y_values=[]
						for line in range(0, len(lines)):
							parts = lines[line].split()
							try:
								y_values.append(float(parts[7])*10**10)
							except:
								pass
						axs[j, i].plot(range(1, len(y_values)+1), y_values, label='R2 Data_' + rep_name + '/' + ff, marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
						#axs[j, i].set_yticks(np.arange(0, max(y_values), 5))
						axs[j, i].set_xticks(np.arange(1, len(y_values)+1, 10))
						#axs[j, i].set_xlim(1, len(y_values)+1)
						#axs[j, i].set_ylim(0, int(max(max(sublist[1:]) for sublist in R2_lists))+10)
						axs[j, i].set_xlabel('Residue number')
						axs[j, i].set_ylabel('Tau_area_value (10**10 s)')
						axs[j, i].set_title('Tau_effective_area')
						axs[j, i].legend()
	plt.tight_layout()
	plt.savefig(relax_folder + output + '.png')
	plt.close()

tau_eff_area(Names, 'Tau_effective_area')
#tau_eff_area(Best_cases_names, 'Accepted_cases/Tau_effective_area_best')


fig, axs = plt.subplots(1, 3, figsize=(15, 6))
for j in Best_cases_names:
	rep_name=j.split("/")[0]
	ff_name=j.split("/")[1]
	index = Names.index(j)
	list=Value_lists[index]

	selected = color_map.get(rep_name, 'black')
	axs[0].plot(range(1, len(list[0])), list[0][1:], label='R1 Data_' + rep_name + '_' + ff_name, marker='o', linestyle='-', lw=1.0, markersize=2, color=selected)
	axs[0].set_xlabel('Residue number')
	axs[0].set_ylabel('R1_values')
	axs[0].set_title('R1_data')
	axs[0].legend()
	axs[1].plot(range(1, len(list[1])), list[1][1:], label='R2 Data_' + rep_name + '_' + ff_name, marker='o', linestyle='-', lw=1.0, markersize=2, color=selected)
	axs[1].set_xlabel('Residue number')
	axs[1].set_ylabel('R2_values')
	axs[1].set_title('R2_data')
	axs[1].legend()
	axs[2].plot(range(1, len(list[2])), list[2][1:], label='NOE Data_' + rep_name + '_' + ff_name, marker='o', linestyle='-', lw=1.0, markersize=2, color=selected)
	axs[2].set_xlabel('Residue number')				
	axs[2].set_ylabel('NOE_values')
	axs[2].set_title('NOE_data')
	axs[2].legend()
		
axs[0].plot(R1_values_exp_filtered[0], R1_values_exp_filtered[1], label='R1 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
axs[0].set_xlabel('Residue number')
axs[0].set_ylabel('R1_values')
axs[0].set_title('R1_data')
axs[0].legend()

axs[1].plot(R2_values_exp_filtered[0], R2_values_exp_filtered[1], label='R2 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
axs[1].set_xlabel('Residue number')
axs[1].set_ylabel('R2_values')
axs[1].set_title('R2_data')
axs[1].legend()

axs[2].plot(NOE_values_exp_filtered[0], NOE_values_exp_filtered[1], label='NOE Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
axs[2].set_xlabel('Residue number')
axs[2].set_ylabel('NOE_values')
axs[2].set_title('NOE_data')
axs[2].legend()
	
plt.tight_layout()
plt.savefig(best_cases_folder + 'Best_cases_plot.png')
plt.close('all')


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
	plt.ylabel('Gyration radius (nm)')
	plt.tight_layout()
	plt.savefig(relax_folder + output + '.png')
	plt.close()

plot_avg_rog_bar(Rog_data_list, "Average_rog_plot")

'''
pdb_data = sorted(glob.glob(SIM_DIR + "/model*/*/"))
cmd.set("ray_opaque_background", 1)

for i in pdb_data[:25]:
	rep_name = i.split('/')[-3]
	forcefield = i.split('/')[-2]
	selected = color_map.get(rep_name, 'black')

	md = glob.glob(i + 'md*smooth*xtc')
	temp = glob.glob(i + 'temp*gro')

	if len(md) > 0 and len(temp) > 0:
		cmd.load(temp[0], 'structure')
		cmd.load_traj(md[0], 'structure', state=1, interval=1000)
		cmd.set('all_states', 'on')
		cmd.color(selected, 'structure')
		cmd.ray(300, 300)
		cmd.png(i + 'Ensemble_' + rep_name + '_' + forcefield + '.png')
		cmd.delete('all')

cmd.quit()

fig, axs = plt.subplots(1, 5, figsize=(15, 6))
current_ax = 0

for item in FORCEFIELDS:
	relax_data = glob.glob(SIM_DIR + "/model*/" + item)
	cmd.set("ray_opaque_background", 1)
	for data in sorted(relax_data):
		rep_name = data.split("/")[-2]

		md = glob.glob(data + '/md*smooth*xtc')
		temp = glob.glob(data + '/temp*gro')

		if len(md) > 0 and len(temp) > 0:
			cmd.load(temp[0], rep_name)
			cmd.load_traj(md[0], rep_name, state=1, interval=2000)
			cmd.color('green', rep_name)
			cmd.set('all_states', 'on')
			cmd.ray(300, 300)
	obj = cmd.get_object_list('all')
	for i in obj[1:]:
		cmd.align(obj[0], i, object='aln', transform=0)
	cmd.png(relax_folder + item + '_aligned_fig.png')
	axs[current_ax].imshow(plt.imread(relax_folder + item + '_aligned_fig.png'), aspect='auto')
	axs[current_ax].set_title(item)	
	current_ax += 1
	cmd.delete('all')
                
cmd.quit()

plt.tight_layout()
plt.savefig(relax_folder + 'aligned_fig.png')

pdb_data = sorted(glob.glob(SIM_DIR + "/model*/*/"))
cmd.set("ray_opaque_background", 1)

for i in pdb_data[:25]:
	rep_name = i.split('/')[-3]
	forcefield = i.split('/')[-2]
	selected = color_map.get(rep_name, 'black')

	md = glob.glob(i + 'md*smooth*xtc')
	temp = glob.glob(i + 'temp*gro')

	if len(md) > 0 and len(temp) > 0:
		cmd.load(temp[0], rep_name)
		cmd.load_traj(md[0], rep_name, state=1, interval=2000)
		cmd.color('green', rep_name)
		cmd.set('all_states', 'on')
		cmd.ray(300, 300)
	obj = cmd.get_object_list('all')
	for i in obj[1:]:
		cmd.align(obj[0], i, object='aln', transform=0)
	cmd.png(i + 'Ensemble_' + rep_name + '_' +  forcefield + '_aligned_fig.png')	
	cmd.delete('all')

cmd.quit()
'''

def plot_ensembles_images(input, output):
	if len(input) > 10:
		fig, axs = plt.subplots(5, 5, figsize=(15, 15))
		for j, ff in enumerate(FORCEFIELDS):
			for k, item in enumerate(input):
				if ff in item:
					data=input[k]
					rep_name = data.split("/")[-3]
					i=int(rep_name[-1])-1
					img = imread(data)
					axs[j, i].imshow(img)
					axs[j, i].set_title(rep_name + '/' + ff)
					axs[j, i].axis('off')

	elif len(input) != 1:
		col=len(input)
		fig, axs = plt.subplots(1, col, figsize=(15, 6))
		for i in range(col):
			data=input[i]
			rep_name = data.split("/")[-3]
			forcefield = data.split("/")[-2]
			img = imread(data)
			axs[i].imshow(img)
			axs[i].set_title(rep_name + '/' + forcefield)
			axs[i].axis('off')
	else:
		data = input[0]
		rep_name = data.split("/")[-3]
		forcefield = data.split("/")[-2]
		img = imread(data)
		plt.imshow(img)
		plt.title(rep_name + '/' + forcefield)
		plt.axis('off')

	
	plt.tight_layout()
	plt.savefig(relax_folder + output + '.png')
	plt.close()

contact_png = sorted(glob.glob(SIM_DIR + "/model*/*/*mdmat*.png"))
plot_ensembles_images(contact_png, "Contact_map_combined")

contact_png_best=[]
for j in Best_cases_names:
	case_search=sorted(glob.glob(SIM_DIR + "/" + j + "/*mdmat*.png"))
	contact_png_best.append(case_search[0])
plot_ensembles_images(contact_png_best, "Accepted_cases/Contact_map_combined_best")


correlation_png=sorted(glob.glob(SIM_DIR + "/model*/*/*correlation*.png"))
plot_ensembles_images(correlation_png, "Correlation_combined")

mdmat_png_best=[]
for j in Best_cases_names:
        case_search=sorted(glob.glob(SIM_DIR + "/" + j + "/*correlation*.png"))
        mdmat_png_best.append(case_search[0])
plot_ensembles_images(mdmat_png_best, "Accepted_cases/Correlation_combined_best")

ensemble_png = sorted(glob.glob(SIM_DIR + "/model*/*/Ensemble_model*.png"))
plot_ensembles_images(ensemble_png, "Ensembles_combined")

aligned_png = sorted(glob.glob(SIM_DIR + "/model*/*/Ensemble*model*aligned_fig.png"))
plot_ensembles_images(aligned_png, "Ensembles_aligned_combined")

aligned_best=[]
for j in Best_cases_names:
        contact_search=sorted(glob.glob(SIM_DIR + "/" + j + "/Ensemble*model*aligned_fig.png"))
        aligned_best.append(case_search[0])
plot_ensembles_images(aligned_best, "Accepted_cases/Ensembles_aligned_best")

'''
mdmat_best=[[], []]
for i in Best_cases_names:
	gro_files = sorted(glob.glob(SIM_DIR + '/' + i + "/temp*2000ns.gro"))
	xtc_files= sorted(glob.glob(SIM_DIR + '/' + i + "/md*2000ns.xtc"))
	mdmat_best[0].append(gro_files[0])
	mdmat_best[1].append(xtc_files[0])

if len(mdmat_best[0]) != len(mdmat_best[1]):
	raise ValueError("The number of .gro files does not match the number of .xtc files.")

trajectories = [mda.Universe(gro, xtc) for gro, xtc in zip(mdmat_best[0], mdmat_best[1])]

average_matrix = np.zeros((len(trajectories[0].select_atoms("name CA").ix) - 1, len(trajectories[0].select_atoms("name CA").ix) - 1))

for u in trajectories:
	CAatoms = u.select_atoms("name CA")

	vec=[]
	matrix = np.empty((len(CAatoms.ix) - 1, len(CAatoms.ix) - 1))

	for i in range(len(CAatoms.ix) - 1):
		vec.append(0)

	for frame in u.trajectory:
		for i in range(len(CAatoms.ix) - 1):
			vec[i] = CAatoms[i + 1].position - CAatoms[i].position
			vec[i] = vec[i] / np.sqrt(np.dot(vec[i], vec[i]))
		for i in range(len(vec)):
			for j in range(len(vec)):
				matrix[i, j] = matrix[i, j] + np.dot(vec[i], vec[j])
	matrix /= len(u.trajectory)

	average_matrix += matrix

average_matrix /= len(trajectories)

w = 10
h = 10
d = 600
plt.figure(figsize=(w, h), dpi=d)
#myFile=np.genfromtxt('result.csv', delimiter=',')

color_map = plt.imshow(average_matrix,vmin=-1, vmax=1)
color_map.set_cmap("seismic")
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=22)
plt.xticks(range(average_matrix.shape[1]), [str(i + 1) for i in range(average_matrix.shape[1])], fontsize=22)

plt.yticks(range(average_matrix.shape[0]), [str(i + 1) for i in range(average_matrix.shape[0])], fontsize=22)

plt.savefig(best_cases_folder + 'Avg_correlation.png', dpi=600) 
plt.close()
all_matrices = []

for i in Best_cases_names:
	csv_file = sorted(glob.glob(SIM_DIR + '/' + i + "/*mdmat*csv"))
	matrix = np.loadtxt(csv_file[0], delimiter=",")
	all_matrices.append(matrix)

average_matrix = sum(all_matrices)

average_matrix /= len(all_matrices)
np.savetxt(best_cases_folder + "Avg_contact.csv", average_matrix[::-1, :], delimiter=",")

fig, ax = plt.subplots()
cmap = plt.cm.viridis
plot = ax.imshow(average_matrix[::-1, :], cmap=cmap, origin="lower")
num_residues = np.shape(average_matrix)[0]
ax.set_xlabel("Residue",fontsize=18)
ax.set_ylabel("Residue",fontsize=18)
plt.xticks(np.arange(num_residues), np.arange(1, num_residues + 1), fontsize=18)
plt.yticks(np.arange(num_residues),np.arange(1, num_residues + 1), fontsize=18)

cbar = plt.colorbar(plot, ax=ax)
cbar.set_label("Distance (nm)",fontsize=18)
plt.savefig(best_cases_folder + "Avg_contact.png", dpi=600)
plt.close()


relax_avg_all=[]
for j in Best_cases_names:
	relax_avg=sorted(glob.glob(SIM_DIR + "/" + j + "/Ensemble*model*aligned_fig.png"))
	relax_avg_all.append(case_search[0])

# Define legend labels and colors for each file
legend_labels = {
    "relaxation_average_GGS.txt": "GGS",
    "relaxation_average_GPS.txt": "GPS",
    "relaxation_average_KAPK.txt": "KAPK",
    # Add more legend labels for other files as needed
}

colors = {
    "relaxation_average_GGS.txt": "orange",
    "relaxation_average_GPS.txt": "tab:blue",
    "relaxation_average_KAPK.txt": "green",
    # Add more colors for other files as needed
}

# Initialize a dictionary to store Tau_eff_area values for each file
tau_values_dict = {}

# Loop through all relaxation files in the folder
for file_name in os.listdir(relax_folder):
    if file_name.startswith("relaxation_average_"):
        file_path = os.path.join(relax_folder, file_name)
        with open(file_path, 'r') as file:
            tau_values = []
            x_values = []  # List to store x-axis values
            lines = file.readlines()
            line_number = 0  # Variable to track line numbers
            for line in lines:
                line_number += 1  # Increment line number
                if 'n' in line:
                    continue
                parts = line.split()
                if len(parts) >= 8:
                    try:
                        tau_value = float(parts[7])
                        tau_values.append(tau_value)
                        x_values.append(line_number)  # Append line number as x value
                    except ValueError:
                        pass  # Skip lines where Tau_eff_area cannot be converted to float
            # Store Tau_eff_area values in the dictionary with the file name as key
            tau_values_dict[file_name] = (x_values, tau_values)

# Plot the Tau_eff_area values
for file_name, (x_values, tau_values) in tau_values_dict.items():
    label = legend_labels.get(file_name, file_name)  # Get legend label or use file name if not defined
    color = colors.get(file_name, None)  # Get color or use None
    plt.plot(x_values, tau_values, marker='o', linestyle='-', lw=1.0, markersize=2, label=label, color=color)

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel('Residue', fontsize=16)
plt.ylabel('Effective correlation time (ns)', fontsize=16)

# Customize the y-axis scalar formatter
plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
plt.gca().yaxis.offsetText.set_fontsize(16)  # Set font size of the 1e-10 label

plt.legend(fontsize=12)
plt.tight_layout()

# Display the plot
plt.show()
'''

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


subprocess.run(["python3", avg_script])



