#!/bin/bash
: '
#SBATCH --time=02:00:00
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --account=Project_2003809
'

module purge
export PATH="/scratch/project_2003809/cmcajsa/env/bin:$PATH"

cd ..
SIM_SCRIPTS=${PWD}

cd ..
SIM_DIR=${PWD}

for dir in $SIM_DIR/Unst*/rep_to_exp_data; do
    mkdir -p $dir  
done


py_script=${SIM_SCRIPTS}/PY_scripts/Old_Relaxations_for_Samuli.py
plot_script=${SIM_SCRIPTS}/PY_scripts/plot_relaxation_data.py
plot_rep_to_exp_script=${SIM_SCRIPTS}/PY_scripts/plot_replicas_to_experiment.py

for i in $SIM_DIR/Unst_snear/*/*/
do
	cd ${i}
	cp $py_script ${i}
	sed -i "s|PATH_TO_CORR|${i}correlation_functions|" ${i}Old_Relaxations_for_Samuli.py
	python3 ${i}Old_Relaxations_for_Samuli.py > relaxation_data.txt 
	#python3 $plot_script
done



python3 $plot_rep_to_exp_script 
