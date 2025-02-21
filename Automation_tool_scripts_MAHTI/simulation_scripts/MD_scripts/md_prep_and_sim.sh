#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --partition=medium
##SBATCH --mem-per-cpu=64
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --account=Project_2003809
##SBATCH --mail-type=END #uncomment to get mail

# this script runs a 256 core (2 full nodes, no hyperthreading) gromacs job, requesting 15 minutes time
# 64 tasks per node, each with 2 OpenMP threads

module purge
module load gcc/9.4.0 openmpi/4.1.2 gromacs/2021.5

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GMX_MAXBACKUP=-1


FORCEFIELD=CHARMM27
SUFFIX=c27
water_model=tip3p
forcefield=charmm27

pwd
var=$(pwd)
basename $(pwd)
PROTEIN="$(basename $PWD)"

MD_BASEDIR=/scratch/project_2003809/cmcajsa/MD-stabilization_puhti
STRUCTURE_DIR=${MD_BASEDIR}/structures
PARAM_DIR=${MD_BASEDIR}/MD_parameter_files

sim_time=sim_time


gmx_mpi pdb2gmx -f ${PROTEIN}.pdb -o ${PROTEIN}.gro -water ${water_model} -ff ${forcefield} -ignh
gmx_mpi editconf -f ${PROTEIN}.gro -o ${PROTEIN}_newbox.gro -c -d 1.0 -bt cubic
gmx_mpi solvate -cp ${PROTEIN}_newbox.gro -cs spc216.gro -o ${PROTEIN}_solv.gro -p topol.top
gmx_mpi grompp -f ${MD_BASEDIR}/MD_parameter_files/${FORCEFIELD}/ions_${SUFFIX}.mdp -c ${PROTEIN}_solv.gro -p topol.top -o ions_${SUFFIX}.tpr
echo SOL | gmx_mpi -quiet genion -s ions_${SUFFIX}.tpr -o ${PROTEIN}_solv_ions.gro -p topol.top -pname NA -nname CL -neutral


srun -n 1 gmx_mpi grompp -f ${MD_BASEDIR}/MD_parameter_files/${FORCEFIELD}/em_${SUFFIX}.mdp -c ${PROTEIN}_solv_ions.gro -p topol.top -o em_${SUFFIX}.tpr
srun gmx_mpi mdrun -deffnm em_${SUFFIX} 

srun -n 1 gmx_mpi grompp -f ${PARAM_DIR}/${FORCEFIELD}/nvt_${SUFFIX}.mdp -c em_${SUFFIX}.gro -r em_${SUFFIX}.gro -p topol.top -o nvt_${SUFFIX}.tpr
srun gmx_mpi mdrun -deffnm nvt_${SUFFIX}


srun -n 1 gmx_mpi grompp -f ${PARAM_DIR}/${FORCEFIELD}/npt_${SUFFIX}.mdp -c nvt_${SUFFIX}.gro -r nvt_${SUFFIX}.gro -t nvt_${SUFFIX}.cpt -p topol.top -o npt_${SUFFIX}.tpr
srun gmx_mpi mdrun -deffnm npt_${SUFFIX}

if [ -f md_run_${PROTEIN}_${SUFFIX}_${sim_time}ns.xtc ] && [ ! -f md_run_${PROTEIN}_${SUFFIX}_${sim_time}ns.gro ]; then
	srun gmx_mpi mdrun -deffnm md_run_${PROTEIN}_${SUFFIX}_${sim_time}ns -dlb yes -cpi md_run_${PROTEIN}_${SUFFIX}_${sim_time}ns.cpt
else
	srun gmx_mpi grompp -f ${PARAM_DIR}/${FORCEFIELD}/md_diff_sim_time/md_${SUFFIX}_${sim_time}ns.mdp -c npt_${SUFFIX}.gro -t npt_${SUFFIX}.cpt -p topol.top -o md_run_${PROTEIN}_$$
	srun gmx_mpi mdrun -deffnm md_run_${PROTEIN}_${SUFFIX}_${sim_time}ns -dlb yes
fi


temp_name=(md_${sim_time}ns.edr)
name=${temp_name%.edr}

echo 1 1 | gmx_mpi trjconv -s ${name}.tpr -f ${name}.xtc -o ${name}_noPBC.xtc -pbc mol -center #Correct for periodic boundary condition (PBC)
echo 4 4 | gmx_mpi rms -s ${name}.tpr -f ${name}_noPBC.xtc -o ${name}_rmsd.xvg -tu ns #compute RMSD


echo -e "Protein\nProtein" | gmx_mpi trjconv -f ${name}.xtc -s ${name}.tpr -pbc mol -center -dump 0 -o temp_${name}.gro
gmx_mpi filter -f md_clean_${name}.xtc -s temp_${name}.gro -nf 20 -all -ol md_smooth_${name}.xtc

echo -e "keep 0\na CA\nname 1 Alpha\nq" | gmx_mpi make_ndx -f ${name}.tpr -o alpha.ndx
echo -e "Alpha\nAlpha" | gmx_mpi rms -f ${name}.xtc -s ${name}.tpr -what rhosc -n alpha.ndx -o ${name}_rhosc.xvg -tu ns

echo 1 | srun gmx_mpi gyrate -s ${name}.tpr -f ${name}_noPBC.xtc -o ${name}_gyrate.xvg
