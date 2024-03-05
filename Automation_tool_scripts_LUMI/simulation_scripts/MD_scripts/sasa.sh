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

MD_BASEDIR=/scratch/project_462000040/cmcajsa/MD-stabilization_puhti
STRUCTURE_DIR=${MD_BASEDIR}/structures
PARAM_DIR=${MD_BASEDIR}/MD_parameter_files


pwd
var=$(pwd)
basename $(pwd)
protein="$(basename $PWD)"


temp_name=(md_*seq*${sim_time}ns.edr)
name=${temp_name%.edr}

		
echo -e "\"non-Water\"" | gmx sasa -f ${name}.xtc -s ${name}.tpr -b ${sim_time}000 -o ${name}_sasa.xvg		

PROT_PATH=${PWD}

cd ..
SIM_DIR=${PWD}
SASA_DIR=${SIM_DIR}/simulation_scripts_${SUFFIX}/sasa_output

cd ${PROT_PATH}

cp ${name}_sasa.xvg ${SASA_DIR}/${sim_time}ns
