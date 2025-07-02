#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=4
#SBATCH --account=project***


export GMX_MAXBACKUP=-1

module use /appl/local/csc/modulefiles
module load gromacs/2023.1-hipsycl



sim_time=sim_time


srun gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -dlb yes



