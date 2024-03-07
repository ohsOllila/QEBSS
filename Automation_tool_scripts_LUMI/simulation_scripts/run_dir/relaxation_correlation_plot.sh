#!/bin/bash
: '
#SBATCH --time=04:00:00
#SBATCH --partition=small
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --account=Project_462000199
'
export PATH="/scratch/project_462000199/cmcajsa/modules/env/bin:$PATH"


cd ..
SIM_SCRIPTS=${PWD}

cd ..
SIM_DIR=${PWD}

for dir in $SIM_DIR/Unst_aphasynuclein/rep_to_exp_data; do
    mkdir -p $dir  
done


py_script=${SIM_SCRIPTS}/PY_scripts/Old_Relaxations_for_Samuli.py
plot_script=${SIM_SCRIPTS}/PY_scripts/plot_relaxation_data.py
plot_rep_to_exp_script=${SIM_SCRIPTS}/PY_scripts/plot_replicas_to_experiment.py

for i in $SIM_DIR/Unst_aphasynuclein/model*/*/
do
	cd ${i}
	cp $py_script ${i}
	sed -i "s|PATH_TO_CORR|${i}correlation_functions|" ${i}Old_Relaxations_for_Samuli.py
	python3 ${i}Old_Relaxations_for_Samuli.py > relaxation_data.txt 
	#python3 $plot_script
done



python3 $plot_rep_to_exp_script 
