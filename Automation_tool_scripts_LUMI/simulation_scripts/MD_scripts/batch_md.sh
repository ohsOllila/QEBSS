#!/bin/bash
#SBATCH --partition=standard-g
##SBATCH --account=project
#SBATCH --account=project_462000404
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=8
#SBATCH --ntasks-per-node=8

module use /appl/local/csc/modulefiles
module load gromacs/2023.3-gpu

export OMP_NUM_THREADS=7

export MPICH_GPU_SUPPORT_ENABLED=1
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_FORCE_GPU_AWARE_MPI=1

cat << EOF > select_gpu
#!/bin/bash

export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
exec \$*
EOF

chmod +x ./select_gpu

CPU_BIND="mask_cpu:fe000000000000,fe00000000000000"
CPU_BIND="${CPU_BIND},fe0000,fe000000"
CPU_BIND="${CPU_BIND},fe,fe00"
CPU_BIND="${CPU_BIND},fe00000000,fe0000000000"

temp_name=(*.pdb)
protein=${temp_name%.pdb}

PARAM_DIR=/scratch/project_462000199/cmcajsa
i=$(basename $PWD)


sim_time=2000
FORCEFIELD=$(basename $PWD)
export GMXLIB=/scratch/project_462000199/cmcajsa/$i

#srun gmx_mpi grompp -f ${PARAM_DIR}/${i}/md_diff_sim_time/md_${sim_time}ns.mdp -c npt_${i}.gro -t npt_${i}.cpt -p topol.top -o md_${sim_time}ns.tpr
srun --cpu-bind=$CPU_BIND ./select_gpu gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -nb gpu -bonded gpu -pme gpu -npme 1 -noappend
gmx_mpi trjcat -f md_1500ns*xtc -o md_2000ns.xtc

#gmx_mpi convert-tpr -s md_${sim_time}ns.tpr -extend 1000000 -o md_1500ns.tpr
#srun --cpu-bind=$CPU_BIND ./select_gpu gmx_mpi mdrun -deffnm md_1500ns -cpi md_1000ns.cpt -noappend -nb gpu -bonded gpu -pme gpu -npme 1
