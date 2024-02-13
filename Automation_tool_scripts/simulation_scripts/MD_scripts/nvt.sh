#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --account=Project_2001058 
##SBATCH --mail-type=END #uncomment to get mail

# this script runs a 256 core (2 full nodes, no hyperthreading) gromacs job, requesting 15 minutes time
# 64 tasks per node, each with 2 OpenMP threads

module purge
module load gcc/9.4.0 openmpi/4.1.2 gromacs/2021.5

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GMX_MAXBACKUP=-1

FORCEFIELD=$(basename $PWD)
export GMXLIB=/scratch/project_2003809/cmcajsa/MD-stabilization/MD_parameter_files/$FORCEFIELD


MD_BASEDIR=/scratch/project_2003809/cmcajsa/MD-stabilization
STRUCTURE_DIR=${MD_BASEDIR}/structures
PARAM_DIR=${MD_BASEDIR}/MD_parameter_files


srun -n 1 gmx_mpi grompp -f ${PARAM_DIR}/${FORCEFIELD}/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
srun gmx_mpi mdrun -deffnm nvt


