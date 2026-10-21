#!/bin/bash
#SBATCH --partition=standard
##SBATCH --account=project
#SBATCH --account=project_462000540
#SBATCH --time=2-00:00:00
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=1

module use /appl/local/csc/modulefiles
module load gromacs/2023.3

export OMP_NUM_THREADS=1

: '

for i in md*2000*xtc; do
	echo 1 | gmx_mpi trjconv -f "$i" -s md_2000ns.tpr -o "$i"
done

echo 1 | gmx_mpi trjconv -f md_1000ns.xtc -s md_1000ns.tpr -o md_1000ns.xtc

#mv md_2000ns.xtc md_2000ns.xtc.cp

rm md_2000ns.xtc

gmx_mpi trjcat -f md*xtc -o md_2000ns.xtc
gmx_mpi check -f md_2000ns.xtc

echo 1 | gmx_mpi trjconv -f md_2000ns.xtc -s md_2000ns.tpr -b 0 -e 2000000 -o md_2000ns.xtc
'

gmx_mpi check -f md_1500ns.xtc
