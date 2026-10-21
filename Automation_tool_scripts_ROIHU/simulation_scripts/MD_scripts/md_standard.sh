#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --nodes=nr_nodes
#SBATCH --ntasks-per-node=128
#SBATCH --account=project
##SBATCH --account=project_2003809

module load gromacs-env

export OMP_NUM_THREADS=1


sim_time=sim_time


#gmx_mpi convert-tpr -s md_${sim_time}ns.tpr -extend 1000000 -o md_2000ns.tpr
#srun gmx_mpi mdrun -deffnm md_2000ns -cpi md_1000ns.cpt -noappend

#srun gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -dlb yes
srun gmx_mpi mdrun -deffnm md -dlb yes

#srun gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -dlb yes -nsteps -1 
#srun gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -dlb yes -noappend

#gmx_mpi mdrun -v -deffnm md_1500ns -cpi md_1500ns.cpt -noappend
#mv md_2000ns.xtc md_2000ns.xtc.back
#gmx_mpi trjcat -f md*xtc -o md_2000ns.xtc

#echo 1 | gmx_mpi trjconv -f md_${sim_time}ns.xtc -s md_${sim_time}ns.tpr -b 0 -e 2000000 -o md_2000ns.xtc
#gmx_mpi check -f md_2000ns.xtc
