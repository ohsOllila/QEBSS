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

FORCEFIELD=CHARMM27
SUFFIX=c27


MD_BASEDIR=/scratch/project_2003809/cmcajsa/MD-stabilization_puhti
STRUCTURE_DIR=${MD_BASEDIR}/structures
PARAM_DIR=${MD_BASEDIR}/MD_parameter_files


srun -n 1 gmx_mpi grompp -f ${PARAM_DIR}/${FORCEFIELD}/nvt_${SUFFIX}.mdp -c em_${SUFFIX}.gro -r em_${SUFFIX}.gro -p topol.top -o nvt_${SUFFIX}.tpr
srun gmx_mpi mdrun -deffnm nvt_${SUFFIX}


