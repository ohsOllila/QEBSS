#!/bin/bash
#SBATCH --partition=standard-g
#SBATCH --account=project
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

sim_time=sim_time
water_model=tip4p

SIM_DIR=${PWD}
PARAM_DIR=$(cd ../MD_parameter_files && pwd)
time_input=$((500000 * $sim_time))


i=$(basename $PWD)
export GMXLIB=$PARAM_DIR/$i

sed "s/time_input/${time_input}/" "${PARAM_DIR}/${i}/md_diff_sim_time/md_any_ns.mdp" > "${PARAM_DIR}/${i}/md_diff_sim_time/md_${sim_time}ns.mdp"
sed -i "s/sim_time/${sim_time}/" "${PARAM_DIR}/${i}/md_diff_sim_time/md_${sim_time}ns.mdp"



temp_name=(*.pdb)
PROTEIN=${temp_name%.pdb}

if [ -f md*gro ]; then
       exit 0
fi


gmx_mpi pdb2gmx -f ${PROTEIN}.pdb -o ${PROTEIN}.gro -water ${water_model} -ff ${i,,} -ignh
gmx_mpi editconf -f ${PROTEIN}.gro -o ${PROTEIN}_newbox.gro -c -d 1.5 -bt dodecahedron

gmx_mpi solvate -cp ${PROTEIN}_newbox.gro -cs ${PARAM_DIR}/${i}/${i,,}.ff/tip4p2005.gro -o ${PROTEIN}_solv.gro -p topol.top

sed -i 's/tip4p.itp/tip4p2005s.itp/' topol.top

gmx_mpi grompp -f ${PARAM_DIR}/${i}/ions.mdp -c ${PROTEIN}_solv.gro -p topol.top -o ions.tpr
echo SOL | gmx_mpi -quiet genion -s ions.tpr -o ${PROTEIN}_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

srun -n 1 gmx_mpi grompp -f ${PARAM_DIR}/${i}/em.mdp -c ${PROTEIN}_solv_ions.gro -p topol.top -o em.tpr  
srun gmx_mpi mdrun -deffnm em

srun -n 1 gmx_mpi grompp -f ${PARAM_DIR}/${i}/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
srun gmx_mpi mdrun -deffnm nvt

srun -n 1 gmx_mpi grompp -f ${PARAM_DIR}/${i}/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
srun gmx_mpi mdrun -deffnm npt

srun gmx_mpi grompp -f ${PARAM_DIR}/${i}/md_diff_sim_time/md_${sim_time}ns.mdp -c npt.gro -t npt.cpt -p topol.top -o md_${sim_time}ns.tpr


srun --cpu-bind=$CPU_BIND ./select_gpu gmx_mpi mdrun -deffnm md_${sim_time}ns -cpi md_${sim_time}ns.cpt -nb gpu -bonded gpu -pme gpu -npme 1


echo 1 1 | gmx_mpi trjconv -f md_${sim_time}ns.xtc -s md_${sim_time}ns.tpr -pbc mol -center -dump 0 -o temp_md_${sim_time}ns.gro
echo 1 1 | gmx_mpi trjconv -f md_${sim_time}ns.xtc -s md_${sim_time}ns.tpr -pbc mol -center -o md_${sim_time}ns_noPBC.xtc
gmx_mpi filter -f md_${sim_time}ns_noPBC.xtc -s temp_md_${sim_time}ns.gro -nf 20 -all -ol md_${sim_time}ns_smooth.xtc
echo 1 | gmx_mpi gyrate -s md_${sim_time}ns.tpr -f md_${sim_time}ns_noPBC.xtc -o md_${sim_time}ns_gyrate.xvg
