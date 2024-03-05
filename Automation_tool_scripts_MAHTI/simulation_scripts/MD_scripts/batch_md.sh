#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=4
#SBATCH --account=project_2005339
##SBATCH --mail-type=END #uncomment to get mail

# this script runs a 256 core (2 full nodes, no hyperthreading) gromacs job, requesting 15 minutes time
# 64 tasks per node, each with 2 OpenMP threads

module purge
module load gromacs-env

export OMP_NUM_THREADS=1
export GMX_MAXBACKUP=0



sim_time=2000


#gmx_mpi convert-tpr -s md_${sim_time}ns.tpr -extend 1000000 -o md_2000ns.tpr
#srun gmx_mpi mdrun -deffnm md_2000ns -cpi md_1000ns.cpt -noappend

#srun gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -dlb yes 
srun gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -dlb yes -noappend

#gmx_mpi mdrun -v -deffnm md_1500ns -cpi md_1500ns.cpt -noappend
#gmx trjcat -f *.xtc -o final_1500ns.xtc


