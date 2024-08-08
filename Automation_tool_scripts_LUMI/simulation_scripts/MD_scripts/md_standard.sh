#!/bin/bash
#SBATCH --partition=standard
##SBATCH --account=project
#SBATCH --account=project_462000540
#SBATCH --time=2-00:00:00
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=6

module use /appl/local/csc/modulefiles
module load gromacs/2023.3

export OMP_NUM_THREADS=1
export GMX_MAXBACKUP=0



sim_time=2000


#gmx_mpi convert-tpr -s md_${sim_time}ns.tpr -extend 1000000 -o md_2000ns.tpr
#srun gmx_mpi mdrun -deffnm md_2000ns -cpi md_1000ns.cpt -noappend

#srun gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -dlb yes
#srun gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -dlb yes -nsteps -1 
#srun gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -dlb yes -noappend

#gmx_mpi mdrun -v -deffnm md_1500ns -cpi md_1500ns.cpt -noappend
mv md_2000ns.xtc md_2000ns.xtc.back
gmx_mpi trjcat -f md*xtc -o md_2000ns.xtc

echo 1 | gmx_mpi trjconv -f md_${sim_time}ns.xtc -s md_${sim_time}ns.tpr -b 0 -e 2000000 -o md_2000ns.xtc
gmx_mpi check -f md_2000ns.xtc
