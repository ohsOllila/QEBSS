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

temp_name=(md_*100ns.xtc)
PROTEIN=${temp_name%.xtc}



echo 1 1 | gmx_mpi trjconv -s $PROTEIN.tpr -f $PROTEIN.xtc -o $PROTEIN_noPBC.xtc -pbc mol -center

echo 4 4 | gmx_mpi trjconv -s $PROTEIN.tpr -f $PROTEIN_noPBC.xtc -dump 100000 -o Frame_100ns_${PROTEIN}.pdb
