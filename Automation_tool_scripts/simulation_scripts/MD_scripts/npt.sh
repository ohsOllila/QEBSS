#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --account=Project_2003809
##SBATCH --mail-type=END #uncomment to get mail

module purge
module load gcc/9.4.0 openmpi/4.1.2 gromacs/2021.5

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GMX_MAXBACKUP=-1

FORCEFIELD=$(basename $PWD)
export GMXLIB=/scratch/project_2003809/cmcajsa/MD-stabilization/MD_parameter_files/$FORCEFIELD

MD_BASEDIR=/scratch/project_2003809/cmcajsa/MD-stabilization
SCRIPT_DIR=${MD_BASEDIR}/scripts
STRUCTURE_DIR=${MD_BASEDIR}/structures
PARAM_DIR=${MD_BASEDIR}/MD_parameter_files




srun -n 1 gmx_mpi grompp -f ${PARAM_DIR}/${FORCEFIELD}/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr

srun gmx_mpi mdrun -deffnm npt

