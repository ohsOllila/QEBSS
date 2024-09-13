#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --partition=small
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --account=Project_462000199


export EBU_USER_PREFIX=/project/project_462000154/EasyBuild
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GMX_MAXBACKUP=-1

module purge
module load LUMI/22.08
module load GROMACS/2021.6-cpeGNU-22.08-CPU


sim_time=sim_time





temp_name=(md_*100ns.edr)
name=${temp_name%.edr}


echo -e "Alpha\nAlpha" | gmx_mpi mdmat -f ${name}.xtc -s ${name}.tpr -mean ${name}_mdmat.xpm
gmx_mpi xpm2ps -f ${name}_mdmat.xpm -o ${name}_mdmat.eps  

