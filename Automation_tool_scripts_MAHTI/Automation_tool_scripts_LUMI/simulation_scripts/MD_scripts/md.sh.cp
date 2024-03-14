#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=4
#SBATCH --account=Project_462000199


export OMP_NUM_THREADS=1
export GMX_MAXBACKUP=-1

module use /appl/local/csc/modulefiles
module load gromacs/2023.1-hipsycl

temp_name=(*.pdb)
protein=${temp_name%.pdb}
NAME=$(echo $protein | sed 's/_1000ns//g')


sim_time=sim_time
FORCEFIELD=$(basename $PWD)
export GMXLIB=/scratch/project_462000199/cmcajsa/$FORCEFIELD


srun gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -dlb yes

: '
TEMP_NAME=(md_*$sim_time*.edr)
PROTEIN=${TEMP_NAME%.edr}


echo -e "Protein\nProtein" | gmx_mpi trjconv -f ${PROTEIN}.xtc -s ${PROTEIN}.tpr -pbc mol -center -dump 0 -o temp_${NAME}.gro 
echo -e "Protein\nProtein" | gmx_mpi trjconv -f ${PROTEIN}.xtc -s ${PROTEIN}.tpr -pbc mol -center -o ${PROTEIN}_noPBC.xtc
gmx_mpi filter -f ${PROTEIN}_noPBC.xtc -s temp_${NAME}.gro -nf 20 -all -ol ${PROTEIN}_smooth.xtc


'
