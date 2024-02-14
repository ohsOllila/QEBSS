#!/usr/bin/python3


import matplotlib.pyplot as plt
import os
import glob
import math
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

snear_data='/scratch/project_2003809/cmcajsa/MD-stabilization/structures/forcefield_compare/simulation_scripts/MD_scripts/snear_exp_data.txt'
SIM_DIR='/scratch/project_2003809/cmcajsa/MD-stabilization/structures/forcefield_compare/Unst_snear/'

relax_folder=SIM_DIR.replace('Unst_snear/', '') + 'results/' + SIM_DIR.split("/")[-2] + '/rep_to_exp_data/'
if not os.path.exists(relax_folder):
    os.makedirs(relax_folder)

'''
# Read data from the text file
with open(snear_data, 'r') as file:
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


FORCEFIELDS=["AMBER03WS", "AMBER99SB-DISP", "AMBER99SBWS", "CHARMM36M", "DESAMBER"]

color_map = {
    'model_01': 'red',
    'model_02': 'blue',
    'model_03': 'green',
    'model_04': 'purple',
    'model_05': 'orange'
}

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

R1_RMSD_values = [item[1].replace('R1: ', '') for item in R1_diff]
R2_RMSD_values = [item[1].replace('R2: ', '') for item in R2_diff]
NOE_RMSD_values = [item[1].replace('NOE: ', '') for item in NOE_diff]


R1_RMSD_min=min(R1_RMSD_values)
R2_RMSD_min=min(R2_RMSD_values)
NOE_RMSD_min=min(NOE_RMSD_values)

rog_values_best=[]
for i in range(max(len(R1_RMSD_values), len(R2_RMSD_values), len(NOE_RMSD_values))):
	if float(R1_RMSD_values[i])/float(R1_RMSD_min) < 1.5 and float(R2_RMSD_values[i])/float(R2_RMSD_min) < 1.5 and float(NOE_RMSD_values[i])/float(NOE_RMSD_min) < 1.5:
		#print(R1_diff[i][0], R1_diff[i][1])

		with open(SIM_DIR+ R1_diff[i][0] +'/md_1000ns_gyrate.xvg', 'r') as file:
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

rog_data=sorted(glob.glob(SIM_DIR+"model*/*"+'/*gyrate*.xvg'))
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
'''
coeff_data=sorted(glob.glob(SIM_DIR+"model*/AMBER03*/Ctimes_Coeffs.txt"))
for data in coeff_data:
	with open(data, 'r') as file:
		lines = file.readlines()
		for line_idx in range(len(lines)):
			parts = lines[line_idx].split()
			if "Res" in str(parts[0]):
				x = float(parts[0].replace('Res_nr_', ''))  # Extracting x-value
				y_num = 1
				while line_idx + y_num < len(lines):
					next_parts = lines[line_idx + y_num].split()
					if "Res" not in next_parts[0]:
						try:
							y = float(next_parts[2].rstrip(','))
							weight=float(next_parts[3])
							plt.scatter(x, y, weight*100, zorder=5, color='red')
							y_num += 1
						except (IndexError, ValueError):
							break
					else:
						break

# Modify the plot settings as needed
plt.xlabel('Residue number')
plt.ylabel('Timescales (ns)')
plt.title('Scatter Plot')

plt.savefig(relax_folder + 'Timescale_plot.png')
plt.close()

'''
pdb_data = sorted(glob.glob(SIM_DIR + "model*/*/"))
cmd.set("ray_opaque_background", 1)

ensemble_images = []  # To store file paths of rendered ensemble images

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
		cmd.png(relax_folder + 'Ensemble_' + rep_name + '_' + forcefield + '.png')
		ensemble_images.append(relax_folder + 'Ensemble_' + rep_name + '_' + forcefield + '.png')

		cmd.delete('all')

fig, axs = plt.subplots(5, 5, figsize=(15, 15))
for idx, img_path in enumerate(ensemble_images):
	img = imread(img_path)
	plt.subplot(5, 5, idx + 1)
	plt.imshow(img)
	plt.axis('off')

plt.tight_layout()
plt.savefig(relax_folder + 'Ensembles_combined.png')
cmd.quit()

'''
coeff_data=sorted(glob.glob(SIM_DIR+"model*/*/Ctimes_Coeffs.txt"))
fig, axs = plt.subplots(5, 5, figsize=(15, 15))
for i in range(5):
	for j in range(5):
		data=coeff_data[i+5*j]
		rep_name = data.split("/")[-3]
		forcefield = data.split("/")[-2]
		with open(data, 'r') as file:
			lines = file.readlines()
			for line_idx in range(len(lines)):
				parts = lines[line_idx].split()
				if "Res" in str(parts[0]):
					x = float(parts[0].replace('Res_nr_', ''))  # Extracting x-value
					y_num = 1
					while line_idx + y_num < len(lines):
						next_parts = lines[line_idx + y_num].split()
						if "Res" not in next_parts[0]:
							try:
								y = float(next_parts[2].rstrip(','))
								weight=float(next_parts[3])
								axs[i, j].scatter(x, y, weight*100, zorder=5, color='red')
								axs[i, j].set_title(rep_name + "/" + forcefield)
								y_num += 1
							except (IndexError, ValueError):
								break
						else:
							break	
plt.tight_layout()
plt.savefig(relax_folder + 'Timescale_plot_all.png')
plt.close()
'''






