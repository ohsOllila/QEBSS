#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=64
#SBATCH --mem-per-cpu=64
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --account=Project_2003809
##SBATCH --mail-type=END #uncomment to get mail

# this script runs a 256 core (2 full nodes, no hyperthreading) gromacs job, requesting 15 minutes ti$
# 64 tasks per node, each with 2 OpenMP threads

module purge
module load gcc/9.4.0 openmpi/4.1.2 gromacs/2021.5
module load python-data

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export GMX_MAXBACKUP=-1


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
