#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=small
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --account=project_462000285
##SBATCH --account=project


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GMX_MAXBACKUP=-1#!/usr/bin/python3

export PATH="/scratch/project_462000285/cmcajsa/systems/forcefield_compare/env/bin:$PATH"

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
py_script=BASE_DIR + "/simulation_scripts/PY_scripts/Old_Relaxations_for_Samuli.py"


best_cases_folder = os.path.join(relax_folder, "Accepted_cases/")

os.makedirs(relax_folder, exist_ok=True)
os.makedirs(best_cases_folder, exist_ok=True)

Best_cases_names=["model_01/AMBER03WS", "model_02/AMBER03WS"]

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
