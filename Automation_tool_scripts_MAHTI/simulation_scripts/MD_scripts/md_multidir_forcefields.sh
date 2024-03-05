#!/bin/bash -l
#SBATCH --partition=standard
#SBATCH --account=project_462000199
#SBATCH --time=2-00:00:00
#SBATCH --nodes=80                             # change according to the number of simulations to be run, keeping in mind that each node can run 8 jobs
#SBATCH --ntasks-per-node=128

module use /appl/local/csc/modulefiles
module load gromacs/2023.1-hipsycl

export OMP_NUM_THREADS=1

list=(U*/*/*/)
list=("${list[@]%/}")


sim_time=sim_time

srun gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -multidir ${list[@]} -dlb yes  


#srun gmx_mpi mdrun -s md_${sim_time}ns.tpr -cpi md_100ns.cpt -append -multidir ${list[@]}

#srun --cpu-bind=$CPU_BIND ./select_gpu gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -nb gpu -bonded gpu -pme gpu -multidir ${list[@]}

