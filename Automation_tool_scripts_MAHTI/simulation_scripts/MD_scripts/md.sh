#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=4
##SBATCH --account=project
#SBATCH --account=Project_2003809


export GMX_MAXBACKUP=-1

module load gromacs-env

sim_time=sim_time

srun gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -dlb yes


