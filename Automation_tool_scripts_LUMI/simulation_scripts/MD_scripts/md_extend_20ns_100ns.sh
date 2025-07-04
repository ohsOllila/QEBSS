#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --partition=small
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --account=Project_462000199


export EBU_USER_PREFIX=/project/project_462000154/EasyBuild
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GMX_MAXBACKUP=-1

module purge
module load LUMI/22.08
module load GROMACS/2021.6-cpeGNU-22.08-CPU

sim_time=100
FORCEFIELD=CHARMM27
SUFFIX=c27

MD_BASEDIR=/scratch/project_2003809/cmcajsa/MD-stabilization_puhti
SCRIPT_DIR=${MD_BASEDIR}/scripts
STRUCTURE_DIR=${MD_BASEDIR}/structures
PARAM_DIR=${MD_BASEDIR}/MD_parameter_files

temp_name=(md_*H*20ns.edr)
PROTEIN=${temp_name%.edr}
EXTENDED=$(echo $PROTEIN | sed 's/20ns/100ns/g')
NAME=$(echo $PROTEIN | sed 's/_20ns//g')

if [ ! -f ${EXTENDED}.xtc ]; then
	srun gmx_mpi convert-tpr -s ${PROTEIN}.tpr -extend 80000 -o ${EXTENDED}.tpr
	srun gmx_mpi mdrun -deffnm ${EXTENDED} -dlb yes
elif [ -f ${EXTENDED}.xtc ] && [ ! -f ${EXTENDED}.gro ]; then
	srun gmx_mpi mdrun -deffnm ${EXTENDED} -cpi ${EXTENDED}.cpt
fi



echo -e "Protein\nProtein" | gmx_mpi trjconv -f ${PROTEIN}.xtc -s ${PROTEIN}.tpr -pbc mol -center -o temp_${NAME}.gro 
echo -e "Protein\nProtein" | gmx_mpi trjconv -f ${EXTENDED}.xtc -s ${EXTENDED}.tpr -pbc mol -center -o ${EXTENDED}_clean.xtc
gmx_mpi filter -f ${EXTENDED}_clean.xtc -s temp_${NAME}.gro -nf 20 -all -ol ${EXTENDED}_smooth.xtc




echo 4 4 | gmx_mpi rms -s ${EXTENDED}.tpr -f ${EXTENDED}_clean.xtc -o ${EXTENDED}_rmsd.xvg -tu ns 

echo -e "keep 0\na CA\nname 1 Alpha\nq" | gmx_mpi make_ndx -f ${EXTENDED}.tpr -o alpha.ndx
echo -e "Alpha\nAlpha" | gmx_mpi rms -f ${EXTENDED}.xtc -s ${EXTENDED}.tpr -what rhosc -n alpha.ndx -o ${EXTENDED}_rhosc.xvg -tu ns

echo -e "Alpha\nAlpha" | gmx_mpi mdmat -f ${EXTENDED}.xtc -s ${EXTENDED}.tpr -mean ${EXTENDED}_mdmat.xpm
gmx_mpi xpm2ps -f ${EXTENDED}_mdmat.xpm -o ${EXTENDED}_mdmat.eps

