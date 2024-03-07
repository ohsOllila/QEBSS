#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=1
#SBATCH --account=Project_2003809
##SBATCH --mail-type=END #uncomment to get mail

# this script runs a 256 core (2 full nodes, no hyperthreading) gromacs job, re$
# 64 tasks per node, each with 2 OpenMP threads

module purge
module load gcc/9.4.0 openmpi/4.1.2 gromacs/2021.5

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export GMX_MAXBACKUP=-1



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

