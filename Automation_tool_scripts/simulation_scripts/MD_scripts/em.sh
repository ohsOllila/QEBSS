#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=64
#SBATCH --mem-per-cpu=1000
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

FORCEFIELD=$(basename $PWD)
export GMXLIB=/scratch/project_2003809/cmcajsa/MD-stabilization/MD_parameter_files/$FORCEFIELD


temp_name=(*.pdb)
PROTEIN=${temp_name%.pdb}


MD_BASEDIR=/scratch/project_2003809/cmcajsa/MD-stabilization
STRUCTURE_DIR=${MD_BASEDIR}/structures
PARAM_DIR=${MD_BASEDIR}/MD_parameter_files

srun -n 1 gmx_mpi grompp -f ${MD_BASEDIR}/MD_parameter_files/${FORCEFIELD}/em.mdp -c ${PROTEIN}_solv_ions.gro -p topol.top -o em.tpr
srun gmx_mpi mdrun -deffnm em
