#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=1
#SBATCH --account=Project_2003809
##SBATCH --mail-type=END #uncomment to get mail

# this script runs a 256 core (2 full nodes, no hyperthreading) gromacs job, re$
# 64 tasks per node, each with 2 OpenMP threads

module purge
module load gcc/9.4.0 openmpi/4.1.2 gromacs/2021.5

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK



sim_time=sim_time





temp_name=(md_*100ns.edr)
name=${temp_name%.edr}


echo -e "Alpha\nAlpha" | gmx_mpi mdmat -f ${name}.xtc -s ${name}.tpr -mean ${name}_mdmat.xpm
gmx_mpi xpm2ps -f ${name}_mdmat.xpm -o ${name}_mdmat.eps  

