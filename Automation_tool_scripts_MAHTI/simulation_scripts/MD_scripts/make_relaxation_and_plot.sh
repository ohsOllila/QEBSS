#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --array=0-24
#SBATCH --account=Project_2001058

module purge
module load pytorch

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GMX_MAXBACKUP=-1

PROT_PATH=${PWD}
PROT_DIR=$(basename $PROT_PATH)

SIM_DIR=$(cd .. && pwd)

py_script=${SIM_DIR}/simulation_scripts/MD_scripts/Old_Relaxations_for_Samuli.py
plot_script=${SIM_DIR}/simulation_scripts/MD_scripts/plot_relaxation_data.py

list=(${PROT_PATH}/*/*/)

cd ${list[${SLURM_ARRAY_TASK_ID}]}

path=${PWD}

cp $py_script ${path}
sed -i "s|PATH_TO_CORR|${path}correlation_functions|" ${path}Old_Relaxations_for_Samuli.py

python3 ${path}Old_Relaxations_for_Samuli.py > relaxation_data.txt 
python3 $plot_script
