#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --partition=medium
#SBATCH --mem-per-cpu=64
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --account=Project_2003809
##SBATCH --mail-type=END #uncomment to get mail

# this script runs a 256 core (2 full nodes, no hyperthreading) gromacs job, requesting 15 minutes time
# 64 tasks per node, each with 2 OpenMP threads

module purge
module load gcc/9.4.0 openmpi/4.1.2 gromacs/2021.5

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GMX_MAXBACKUP=-1


sim_time=sim_time

SCRIPTS=${PWD}
mkdir ${SCRIPTS}/rog_output
ROG_DIR=${SCRIPTS}/rog_output
mkdir ${ROG_DIR}/${sim_time}ns
DIR=${ROG_DIR}/${sim_time}ns

cd ..

SIM_DIR=${PWD}

PROTEIN=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${SIM_DIR}/MD_struct_directories.txt)
file_path=${PROTEIN}/md*${sim_time}ns*tpr
file=$(basename $file_path)
name=${file%.tpr}

cd $PROTEIN


echo 1 | srun gmx_mpi gyrate -s ${name}.tpr -f ${name}_noPBC.xtc -o ${name}_gyrate.xvg


cp ${name}_gyrate.xvg ${DIR}

