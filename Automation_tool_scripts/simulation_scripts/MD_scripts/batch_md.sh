#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=5
#SBATCH --account=Project_2001058
##SBATCH --mail-type=END #uncomment to get mail

# this script runs a 256 core (2 full nodes, no hyperthreading) gromacs job, requesting 15 minutes time
# 64 tasks per node, each with 2 OpenMP threads

module purge
module load gromacs-env

export OMP_NUM_THREADS=1
export GMX_MAXBACKUP=0



sim_time=1500



srun gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -dlb yes 

#cp md_1000ns.xtc md_1000ns.xtc.cpy
#cp md_1000ns.gro md_1000ns.gro.cpy

#gmx_mpi convert-tpr -s md_1000ns.tpr -o md_1500ns.tpr -until 1500000
#gmx_mpi mdrun -v -deffnm md_1500ns -cpi md_1500ns.cpt -noappend
#gmx trjcat -f *.xtc -o final_1500ns.xtc


