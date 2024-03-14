#!/usr/bin/python3


import matplotlib.pyplot as plt
import os
import glob
import math
import numpy as np
from collections import Counter
from matplotlib.image import imread
import statistics
from pymol import cmd
import pymol

R1_values_exp = [[], []]
R2_values_exp = [[], []]
NOE_values_exp = [[], []]

#snear_data='/scratch/project_2003809/cmcajsa/MD-stabilization/structures/forcefield_compare/simulation_scripts/MD_scripts/snear_exp_data.txt'
#SIM_DIR='/scratch/project_2003809/cmcajsa/MD-stabilization/structures/forcefield_compare/Unst_snear/'

#exp_data='/scratch/project_2003809/cmcajsa/MD-stabilization/structures/forcefield_compare/simulation_scripts/MD_scripts/hydrolase_exp_data.txt'
#SIM_DIR='/scratch/project_2003809/cmcajsa/MD-stabilization/structures/forcefield_compare/Unst_hydrolase/'
SIM_DIR = os.getcwd()
BASE_DIR=os.path.dirname(SIM_DIR)
exp_data = BASE_DIR + '/Unst_hydrolase_exp_data.txt'

relax_folder=SIM_DIR.replace('Unst_hydrolase', '') + "results/" + SIM_DIR.split("/")[-1] + '/rep_to_exp_data/'
print(relax_folder)

if not os.path.exists(relax_folder):
    os.makedirs(relax_folder)

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
			R2_values_exp[1].append(float(parts[2]))
			R2_values_exp[0].append(line)
		except ValueError:
                        R2_values_exp[1].append("n")
                        R2_values_exp[0].append("n")
		try:
			NOE_values_exp[1].append(float(parts[3]))
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


FORCEFIELDS=["AMBER03WS", "AMBER99SB-DISP", "AMBER99SBWS", "CHARMM36M", "DESAMBER"]

color_map = {
    'model_01': 'red',
    'model_02': 'blue',
    'model_03': 'green',
    'model_04': 'purple',
    'model_05': 'orange'
}

combined_list_of_lists=[]

R1_diff=[]
R2_diff=[]
NOE_diff=[]

for item in FORCEFIELDS:
	fig, axs = plt.subplots(1, 3, figsize=(15, 6))
	relax_data=glob.glob(SIM_DIR+"/model*/"+item+'/relaxation_data.txt')
	for data in sorted(relax_data):
		rep_name = data.split("/")[-3]
		combined_list=[["R1"], ["R2"], ["NOE"]]
		R1_values_diff=[]
		R2_values_diff=[]
		NOE_values_diff=[]
		with open(data, 'r') as file:
			lines = file.readlines()
			for line in range(0, len(lines)):
				parts = lines[line].split()
				try:
					combined_list[0].append(1 / float(parts[1]))
					if R1_values_exp[1][line] != "n":
						R1_values_diff.append((1 / float(parts[1])-R1_values_exp[1][line]))
				except:
					pass
				try:
					combined_list[1].append(1 / float(parts[3]))
					if R2_values_exp[1][line] != "n":
						R2_values_diff.append((1 / float(parts[3])-R2_values_exp[1][line]))
				except:	
					pass
				try:
					combined_list[2].append(float(parts[5]))
					if NOE_values_exp[1][line] != "n":
						NOE_values_diff.append((float(parts[5])-NOE_values_exp[1][line]))
				except:
					pass


		if sum(R1_values_diff) != 0:
			RMSD_R1=f"R1: {math.sqrt(sum(x**2 for x in R1_values_diff) / len(R1_values_diff))}"
			R1_diff.append((rep_name + "/" + item, RMSD_R1, R1_values_diff))
		if sum(R2_values_diff) != 0:
			RMSD_R2=f"R2: {math.sqrt(sum(x**2 for x in R2_values_diff) / len(R2_values_diff))}"
			R2_diff.append((rep_name + "/" + item, RMSD_R2, R2_values_diff))
		if sum(NOE_values_diff) != 0:
			RMSD_NOE=f"NOE: {math.sqrt(sum(x**2 for x in NOE_values_diff) / len(NOE_values_diff))}"
			NOE_diff.append((rep_name + "/" + item, RMSD_NOE, NOE_values_diff))
		try:
			combined_list_of_lists.append((rep_name + "/" + item, combined_list))
		except:
			pass
			
		selected = color_map.get(rep_name, 'black')			

		axs[0].plot(range(1, len(combined_list[0])), combined_list[0][1:], label='R1 Data_' + rep_name, marker='o', linestyle='-', lw=1.0, markersize=2, color=selected)
		axs[0].set_xlabel('Residue number')
		axs[0].set_ylabel('R1_values')
		axs[0].set_title('R1_data')
		axs[0].legend()
		axs[1].plot(range(1, len(combined_list[1])), combined_list[1][1:], label='R2 Data_' + rep_name, marker='o', linestyle='-', lw=1.0, markersize=2, color=selected)
		axs[1].set_xlabel('Residue number')
		axs[1].set_ylabel('R2_values')
		axs[1].set_title('R2_data')
		axs[1].legend()
		axs[2].plot(range(1, len(combined_list[2])), combined_list[2][1:], label='NOE Data_' + rep_name, marker='o', linestyle='-', lw=1.0, markersize=2, color=selected)
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
	plt.savefig(relax_folder + 'relaxation_' + item + '_plot.png')
plt.close('all')

Names=[item[0] for item in combined_list_of_lists]
Value_lists=[item[1] for item in combined_list_of_lists]

R1_lists=[item[0] for item in Value_lists]
R2_lists=[item[1] for item in Value_lists]
NOE_lists=[item[2] for item in Value_lists]

R1_max=max(max(sublist[1:]) for sublist in R1_lists)
R2_max=max(max(sublist[1:]) for sublist in R2_lists)
NOE_max=max(max(sublist[1:]) for sublist in NOE_lists)

NOE_min=min(min(sublist[1:]) for sublist in NOE_lists)

fig, axs = plt.subplots(5, 5, figsize=(30, 15))
for i in range(5):
	for j in range(5):
		data=Names[i+5*j]
		axs[j, i].plot(range(1, len(R1_lists[i+5*j])), R1_lists[i+5*j][1:], label='R1 Data_' + data, marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
		axs[j, i].plot(R1_values_exp_filtered[0], R1_values_exp_filtered[1], label='R1 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
		axs[j, i].set_yticks(np.arange(0, int(R1_max)+0.5, 0.5))
		axs[j, i].set_xticks(np.arange(1, len(R1_lists[i+5*j]), 20))
		axs[j, i].set_xlim(1, len(R1_lists[i+5*j]))
		axs[j, i].set_ylim(0, int(R1_max)+0.5)
		axs[j, i].set_xlabel('Residue number')
		axs[j, i].set_ylabel('R1_values')
		axs[j, i].set_title('R1_data')
		axs[j, i].legend()
plt.tight_layout()
plt.savefig(relax_folder + 'R1_relaxation_combined_plot.png')
plt.close()

fig, axs = plt.subplots(5, 5, figsize=(30, 15))
for i in range(5):
	for j in range(5):
		data=Names[i+5*j]
		axs[j, i].plot(range(1, len(R2_lists[i+5*j])), R2_lists[i+5*j][1:], label='R2 Data_' + data, marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
		axs[j, i].plot(R2_values_exp_filtered[0], R2_values_exp_filtered[1], label='R2 Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
		axs[j, i].set_yticks(np.arange(0, int(R2_max)+2, 2))
		axs[j, i].set_xticks(np.arange(1, len(R2_lists[i+5*j]), 20))
		axs[j, i].set_xlim(1, len(R2_lists[i+5*j]))
		axs[j, i].set_ylim(0, int(R2_max)+2)
		axs[j, i].set_xlabel('Residue number')
		axs[j, i].set_ylabel('R2_values')
		axs[j, i].set_title('R2_data')
		axs[j, i].legend()
plt.tight_layout()
plt.savefig(relax_folder + 'R2_relaxation_combined_plot.png')
plt.close()

fig, axs = plt.subplots(5, 5, figsize=(30, 15))
for i in range(5):
	for j in range(5):
		data=Names[i+5*j]
		axs[j, i].plot(range(1, len(NOE_lists[i+5*j])), NOE_lists[i+5*j][1:], label='NOE Data_' + data, marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
		axs[j, i].plot(NOE_values_exp_filtered[0], NOE_values_exp_filtered[1], label='NOE Data', marker='o', linestyle='-', lw=1.0, markersize=2, color='black')
		axs[j, i].set_yticks(np.arange(int(NOE_min)-0.5, int(NOE_max)+1, 0.5))
		axs[j, i].set_xticks(np.arange(1, len(NOE_lists[i+5*j]), 20))
		axs[j, i].set_xlim(1, len(NOE_lists[i+5*j]))
		axs[j, i].set_ylim(int(NOE_min)-0.5, int(NOE_max)+1)
		axs[j, i].set_xlabel('Residue number')
		axs[j, i].set_ylabel('NOE_values')
		axs[j, i].set_title('NOE_data')
		axs[j, i].legend()
plt.tight_layout()
plt.savefig(relax_folder + 'NOE_relaxation_combined_plot.png')
plt.close()

R1_RMSD_values = [item[1].replace('R1: ', '') for item in R1_diff]
R2_RMSD_values = [item[1].replace('R2: ', '') for item in R2_diff]
NOE_RMSD_values = [item[1].replace('NOE: ', '') for item in NOE_diff]


R1_RMSD_min=min(R1_RMSD_values)
R2_RMSD_min=min(R2_RMSD_values)
NOE_RMSD_min=min(NOE_RMSD_values)

rog_values_best=[]
Best_names=[]
for i in range(max(len(R1_RMSD_values), len(R2_RMSD_values), len(NOE_RMSD_values))):
	if float(R1_RMSD_values[i])/float(R1_RMSD_min) < 1.5 and float(R2_RMSD_values[i])/float(R2_RMSD_min) < 1.5 and float(NOE_RMSD_values[i])/float(NOE_RMSD_min) < 1.5:
		#print(R1_diff[i][0], R1_diff[i][1])
		Best_names.append((R1_diff[i][0], R1_diff[i][1]))

		with open(SIM_DIR + '/' + R1_diff[i][0] +'/md_2000ns_gyrate.xvg', 'r') as file:
			lines = file.readlines()
			for line in range(27, len(lines)):
				parts = lines[line].split()
				rog_values_best.append(round(float(parts[1]), 1))
value_counts = Counter(rog_values_best)
values, counts = zip(*sorted(value_counts.items(), key=lambda x: x[0]))
plt.plot(values, counts)
plt.xlabel('Radius of gyration')
plt.ylabel('Count')
plt.title('Count of Radius of gyration values from best simulation data')
plt.tight_layout()
plt.savefig(relax_folder + 'best_rog_density_landscape_plot.png')
plt.close()


rog_data=sorted(glob.glob(SIM_DIR+"/model*/*"+'/*gyrate*.xvg'))
fig, axs = plt.subplots(5, 5, figsize=(15, 15))

for i in range(5):
	for j in range(5):
		rog_values=[]
		data=rog_data[i+5*j]
		rep_name = data.split("/")[-3]
		forcefield = data.split("/")[-2]
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
#		axs[j, i].set_yticks(np.arange(0, 70000, 10000))
		axs[j, i].set_xticks(np.arange(1.5, 7, 1))
		axs[j, i].set_xlim(1.5, 7)
#		axs[j, i].set_ylim(0, 70000)

		axs[i, j].plot(values, counts)
		axs[i, j].axvline(x=statistics.mean(rog_values), color='black', linestyle='--', label='Mean Rog value')
		axs[i, j].set_title(rep_name + "/" + forcefield)


plt.tight_layout()
plt.savefig(relax_folder + 'density_landscape_plot.png')
plt.close()

used_labels = [[], [], []]
fig, axs = plt.subplots(1, 3, figsize=(15, 6))
for i, ff in enumerate(FORCEFIELDS):
	list=[[], [], []]
	for item in R1_diff:
		if ff in item[0]:
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
		if ff in item[0]:
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
		if ff in item[0]:
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
plt.savefig(relax_folder + 'Difference_to_experiment_plot.png')
plt.close()

'''
pdb_data = sorted(glob.glob(SIM_DIR + "/model*/*/"))
cmd.set("ray_opaque_background", 1)

for i in pdb_data[:25]:
	rep_name = i.split('/')[-3]
	forcefield = i.split('/')[-2]
	md = glob.glob(i + 'md*smooth*xtc')
	temp = glob.glob(i + 'temp*gro')

	if len(md) > 0 and len(temp) > 0:
		cmd.load(temp[0], 'structure')
		cmd.load_traj(md[0], 'structure', state=1, interval=1000)
		cmd.set('all_states', 'on')
		cmd.ray(300, 300)
		cmd.png(i + 'Ensemble_' + rep_name + '_' + forcefield + '.png')
		cmd.delete('all')

cmd.quit()

ensemble_images = sorted(glob.glob(SIM_DIR + "/model*/*/Ensemble_model*.png"))
fig, axs = plt.subplots(5, 5, figsize=(15, 15))
for i in range(5):
        for j in range(5):
                data=ensemble_images[i+5*j]
                rep_name = data.split("/")[-3]
                forcefield = data.split("/")[-2]
                try:
                        img = imread(data)
                        axs[i, j].imshow(img)
                        axs[i, j].set_title(rep_name + '/' + forcefield)
                        axs[i, j].axis('off')
                except:
                        pass

plt.tight_layout()
plt.savefig(relax_folder + 'Ensembles_combined.png')
plt.close()

'''

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
		axs.set_title(rep_name + "/" + forcefield)


coeff_data = sorted(glob.glob(SIM_DIR + "/model*/*/Ctimes_Coeffs.txt"))

for data in coeff_data:
	create_timescale_scatter_plot(data)
plt.tight_layout()
plt.savefig(relax_folder + 'Timescale_plot.png')
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

'''
ensemble_images=sorted(glob.glob(SIM_DIR+"/model*/*"+'/*correlation*.png'))
fig, axs = plt.subplots(5, 5, figsize=(15, 15))
for i in range(5):
	for j in range(5):
		data=ensemble_images[i+5*j]
		rep_name = data.split("/")[-3]
		forcefield = data.split("/")[-2]
		try:
			img = imread(data)
			axs[i, j].imshow(img)
			axs[i, j].set_title(rep_name + '/' + forcefield)
			axs[i, j].axis('off')
		except:
			pass
plt.tight_layout()
plt.savefig(relax_folder + 'Correlation_combined.png')
'''

ensemble_images = sorted(glob.glob(SIM_DIR + "/model*/*/*mdmat*.png"))
fig, axs = plt.subplots(5, 5, figsize=(15, 15))
for i in range(5):
	for j in range(5):
		data=ensemble_images[5*i+j]
		rep_name = data.split("/")[-3]
		forcefield = data.split("/")[-2]
		try:
			img = imread(data)
			axs[j, i].imshow(img)
			axs[j, i].set_title(rep_name + '/' + forcefield)
			axs[j, i].axis('off')
		except:
			pass

plt.tight_layout()
plt.savefig(relax_folder + 'Contact_map_combined.png')


relax_data=glob.glob(SIM_DIR+"/model*/*"+'/relaxation_data.txt')
fig, axs = plt.subplots(5, 5, figsize=(30, 15))
for i in range(5):
	for j in range(5):
		data=relax_data[5*i+j]
		rep_name = data.split("/")[-3]
		forcefield = data.split("/")[-2]
		with open(data, 'r') as file:
			lines = file.readlines()
			y_values=[]
			for line in range(0, len(lines)):
				parts = lines[line].split()
				try:
					y_values.append(float(parts[7])*10**10)
				except:
					pass
			axs[j, i].plot(range(1, len(y_values)+1), y_values, label='R2 Data_' + rep_name + '/' + forcefield, marker='o', linestyle='-', lw=1.0, markersize=2, color="green")
			#axs[j, i].set_yticks(np.arange(0, max(y_values), 5))
			axs[j, i].set_xticks(np.arange(1, len(y_values)+1, 10))
			#axs[j, i].set_xlim(1, len(y_values)+1)
			#axs[j, i].set_ylim(0, int(R2_max)+10)
			axs[j, i].set_xlabel('Residue number')
			axs[j, i].set_ylabel('Tau_area_value (Å)')
			axs[j, i].set_title('Tau_effective_area')
			axs[j, i].legend()
plt.tight_layout()
plt.savefig(relax_folder + 'Tau_effective_area.png')




print(relax_data)
