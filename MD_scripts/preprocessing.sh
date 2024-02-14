#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --partition=medium
#SBATCH --mem-per-cpu=64
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --account=Project_2003809
##SBATCH --mail-type=END #uncomment to get mail

# this script runs a 256 core (2 full nodes, no hyperthreading) gromacs job, requesting 15 minutes time
# 64 tasks per node, each with 2 OpenMP threads

module purge
module load gromacs-env

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GMX_MAXBACKUP=-1


FORCEFIELD=$(basename $PWD)

export GMXLIB=/scratch/project_2003809/cmcajsa/MD-stabilization_puhti/MD_parameter_files/$FORCEFIELD

if [ "$FORCEFIELD" = "AMBER03WS" ]; then
    watermodel1="tip4p2005.itp"
    watermodel2="tip4p2005.gro"
elif [ "$FORCEFIELD" = "AMBER99SB-DISP" ]; then
    watermodel1="a99SBdisp_water.itp"
    watermodel2="a99SBdisp_water.gro"
elif [ "$FORCEFIELD" = "AMBER99SBWS" ]; then
    watermodel1="tip4p2005.itp"
    watermodel2="tip4p2005.gro"
elif [ "$FORCEFIELD" = "DESAMBER" ]; then
    watermodel1="tip4pd.itp" 
    watermodel2="tip4p.gro"
elif [ "$FORCEFIELD" = "CHARMM36M" ]; then
    watermodel2="spc216.gro"
fi

temp_name=(*.pdb)
PROTNAME=${temp_name%.pdb}

MD_BASEDIR=/scratch/project_2003809/cmcajsa/MD-stabilization_puhti
STRUCTURE_DIR=${MD_BASEDIR}/structures
PARAM_DIR=${MD_BASEDIR}/MD_parameter_files


if [ "$FORCEFIELD" = "CHARMM36M" ]; then
	gmx_mpi pdb2gmx -f ${PROTNAME}.pdb -o ${PROTNAME}.gro -water tip3p -ff ${FORCEFIELD,,} -ignh
else
	gmx_mpi pdb2gmx -f ${PROTNAME}.pdb -o ${PROTNAME}.gro -water tip4p -ff ${FORCEFIELD,,} -ignh
fi

gmx_mpi editconf -f ${PROTNAME}.gro -o ${PROTNAME}_newbox.gro -c -d 1.5 -bt dodecahedron


gmx_mpi solvate -cp ${PROTNAME}_newbox.gro -cs $PARAM_DIR/$FORCEFIELD/${FORCEFIELD,,}.ff/${watermodel2} -o ${PROTNAME}_solv.gro -p topol.top
sed -i s/tip4p.itp/$watermodel1/ topol.top

gmx_mpi grompp -f ${MD_BASEDIR}/MD_parameter_files/${FORCEFIELD}/ions.mdp -c ${PROTNAME}_solv.gro -p topol.top -o ions.tpr



if [ "$FORCEFIELD" = "CHARMM36M" ]; then
        echo SOL | gmx_mpi -quiet genion -s ions.tpr -o ${PROTNAME}_solv_ions.gro -p topol.top -pname SOD -nname CLA -neutral
else
    	echo SOL | gmx_mpi -quiet genion -s ions.tpr -o ${PROTNAME}_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
fi
