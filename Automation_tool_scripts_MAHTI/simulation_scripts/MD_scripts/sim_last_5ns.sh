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
svans=$(echo $PROTEIN | sed 's/100ns/5ns_tail/g')


echo 1 | gmx_mpi trjconv -s $PROTEIN.tpr -f $PROTEIN.xtc -b 95000 -o $svans.xtc
echo 1 1 | gmx_mpi trjconv -s $PROTEIN.tpr -f $PROTEIN.xtc -b 95000 -pbc mol -center -dump 0 -o $svans.gro

echo 1 1 | gmx_mpi trjconv -s $PROTEIN.tpr -f $PROTEIN.xtc -o $PROTEIN_noPBC.xtc -pbc mol -center
echo 4 4 | gmx_mpi rms -s $PROTEIN.tpr -f $PROTEIN_noPBC.xtc -b 95000 -m matrix.xpm

smallest_average=$(python3 ../simulation_scripts/PY_scripts/determine_smallest_frame.py)
echo $smallest_average

smallest_frame_time=$((95000 + 10*$smallest_average))
echo $smallest_frame_time

echo 4 4 | gmx_mpi trjconv -s $PROTEIN.tpr -f $PROTEIN_noPBC.xtc -dump $smallest_frame_time -o avg_$svans.pdb
