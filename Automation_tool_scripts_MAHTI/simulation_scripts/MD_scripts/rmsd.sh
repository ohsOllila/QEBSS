#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --partition=small
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --account=Project_462000199


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GMX_MAXBACKUP=-1


module use /appl/local/csc/modulefiles
module load gromacs/2023.1-hipsycl

sim_time=sim_time





temp_name=(md_*${sim_time}ns.edr)
name=${temp_name%.edr}

		
echo 4 4 | gmx_mpi rms -s ${name}.tpr -f ${name}_noPBC.xtc -o ${name}_rmsd.xvg -tu ns #compute RMSD
		

