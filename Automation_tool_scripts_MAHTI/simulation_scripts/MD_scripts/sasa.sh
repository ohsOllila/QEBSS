#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --partition=medium
#SBATCH --mem-per-cpu=64
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
