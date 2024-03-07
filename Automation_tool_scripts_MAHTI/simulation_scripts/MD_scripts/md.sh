#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=4
#SBATCH --account=project
##SBATCH --mail-type=END #uncomment to get mail

# this script runs a 256 core (2 full nodes, no hyperthreading) gromacs job, requesting 15 minutes time
# 64 tasks per node, each with 2 OpenMP threads

module purge
module load gromacs-env

export OMP_NUM_THREADS=1
export GMX_MAXBACKUP=0



sim_time=sim_time


srun gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -dlb yes



