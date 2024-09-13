#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --partition=small
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=1
#SBATCH --account=Project_462000199


export OMP_NUM_THREADS=1
export GMX_MAXBACKUP=-1

module use /appl/local/csc/modulefiles
module load gromacs/2023.1-hipsycl

temp_name=(*.pdb)
PROTEIN=${temp_name%.pdb}
NAME=$(echo $protein | sed 's/_1000ns//g')


sim_time=1000
FORCEFIELD=$(basename $PWD)
export GMXLIB=/scratch/project_462000199/cmcajsa/$FORCEFIELD


echo -e "Protein\nProtein" | gmx_mpi trjconv -f ${PROTEIN}.xtc -s ${PROTEIN}.tpr -pbc mol -center -o ${PROTEIN}_noPBC.xtc
echo 1 | srun gmx_mpi gyrate -s ${PROTEIN}.tpr -f ${PROTEIN}_noPBC.xtc -o ${PROTEIN}_gyrate.xvg
