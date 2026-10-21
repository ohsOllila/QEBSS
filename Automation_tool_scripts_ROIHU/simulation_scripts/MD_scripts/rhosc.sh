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
sim_time=sim_time

MD_BASEDIR=/scratch/project_2003809/cmcajsa/MD-stabilization_puhti
STRUCTURE_DIR=${MD_BASEDIR}/structures
PARAM_DIR=${MD_BASEDIR}/MD_parameter_files


pwd
var=$(pwd)
basename $(pwd)
protein="$(basename $PWD)"


temp_name=(md_*Un*${sim_time}ns.edr)
name=${temp_name%.edr}

echo -e "keep 0\na CA\nname 1 Alpha\nq" | gmx_mpi make_ndx -f ${name}.tpr -o alpha.ndx
echo -e "Alpha\nAlpha" | gmx_mpi rms -f ${name}.xtc -s ${name}.tpr -what rhosc -n alpha.ndx -o ${name}_rhosc.xvg -tu ns

		

PROT_PATH=${PWD}

cd ..
SIM_DIR=${PWD}
RHOSC_DIR=${SIM_DIR}/simulation_scripts/rhosc_output
RHOSC=${RHOSC_DIR}/${sim_time}ns

cd ${PROT_PATH}

cp ${name}_rhosc.xvg ${RHOSC}

