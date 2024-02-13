#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --account=Project_2003809
##SBATCH --mail-type=END #uncomment to get mail

# this script runs a 256 core (2 full nodes, no hyperthreading) gromacs job, requesting 15 minutes time
# 64 tasks per node, each with 2 OpenMP threads

module purge
module load gcc/9.4.0 openmpi/4.1.2 gromacs/2021.5

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GMX_MAXBACKUP=-1


sim_time=1000





temp_name=(md_*${sim_time}ns.edr)
name=${temp_name%.edr}

echo 1 1 | gmx_mpi trjconv -s ${name}.tpr -f ${name}.xtc -o ${name}_noPBC.xtc -pbc mol -center #Correct for periodic boundary condition (PBC)
echo 4 4 | gmx_mpi rms -s ${name}.tpr -f ${name}_noPBC.xtc -o ${name}_rmsd.xvg -tu ns #compute RMSD


echo -e "Protein\nProtein" | gmx_mpi trjconv -f ${name}.xtc -s ${name}.tpr -pbc mol -center -dump 0 -o temp_${name}.gro
gmx_mpi filter -f ${name}_noPBC.xtc -s temp_${name}.gro -nf 20 -all -ol md_smooth_${name}.xtc


echo 1 | srun gmx_mpi gyrate -s ${name}.tpr -f ${name}_noPBC.xtc -o ${name}_gyrate.xvg		
		

