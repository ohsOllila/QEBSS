#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=4
#SBATCH --account=Project_2001058
##SBATCH --mail-type=END #uncomment to get mail



module purge
module load gromacs-env

export OMP_NUM_THREADS=1
export GMX_MAXBACKUP=0


sim_time=1000


srun gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -dlb yes 

