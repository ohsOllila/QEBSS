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
water_model=tip3p
forcefield=charmm27

temp_name=(*.pdb)
PROTNAME=${temp_name%.edr}



MD_BASEDIR=/scratch/project_2003809/cmcajsa/MD-stabilization_puhti
STRUCTURE_DIR=${MD_BASEDIR}/structures
PARAM_DIR=${MD_BASEDIR}/MD_parameter_files



gmx_mpi pdb2gmx -f ${PROTNAME}.pdb -o ${PROTNAME}.gro -water ${water_model} -ff ${forcefield} -ignh
gmx_mpi editconf -f ${PROTNAME}.gro -o ${PROTNAME}_newbox.gro -c -d 1.0 -bt cubic
gmx_mpi solvate -cp ${PROTNAME}_newbox.gro -cs spc216.gro -o ${PROTNAME}_solv.gro -p topol.top
gmx_mpi grompp -f ${MD_BASEDIR}/MD_parameter_files/${FORCEFIELD}/ions_${SUFFIX}.mdp -c ${PROTNAME}_solv.gro -p topol.top -o ions_${SUFFIX}.tpr
echo SOL | gmx_mpi -quiet genion -s ions_${SUFFIX}.tpr -o ${PROTNAME}_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

