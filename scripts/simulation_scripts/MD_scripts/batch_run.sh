#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=6
#SBATCH --account=Project_2001058
##SBATCH --mail-type=END #uncomment to get mail

# this script runs a 256 core (2 full nodes, no hyperthreading) gromacs job, requesting 15 minutes time
# 64 tasks per node, each with 2 OpenMP threads

module purge
module load gromacs-env

export OMP_NUM_THREADS=1
export GMX_MAXBACKUP=0


temp_name=(*.pdb)
protein=${temp_name%.pdb}


sim_time=1000
FORCEFIELD=$(basename $PWD)
export GMXLIB=/scratch/project_2003809/cmcajsa/MD-stabilization/MD_parameter_files/$FORCEFIELD

MD_BASEDIR=/scratch/project_2003809/cmcajsa/MD-stabilization
SCRIPT_DIR=${MD_BASEDIR}/scripts
STRUCTURE_DIR=${MD_BASEDIR}/structures
PARAM_DIR=${MD_BASEDIR}/MD_parameter_files


#srun gmx_mpi grompp -f ${PARAM_DIR}/${FORCEFIELD}/md_diff_sim_time/md_${sim_time}ns.mdp -c npt.gro -t npt.cpt -p topol.top -o md_run_${sim_time}ns.tpr

srun gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -dlb yes 

