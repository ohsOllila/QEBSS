#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --array=0-num_jobs
#SBATCH --output=array_job_output_%A_%a.txt
#SBATCH --account=project


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GMX_MAXBACKUP=-1

module load gromacs-env



sim_time=sim_time
water_model=tip4p

SIM_DIR=${PWD}
PARAM_DIR=$(cd ../MD_parameter_files && pwd)
time_input=$((500000 * $sim_time))


FORCEFIELD=("AMBER03WS" "AMBER99SB-DISP" "AMBER99SBWS" "CHARMM36M" "DESAMBER")


PROT_FOLDERS=()

for pdb_file in $SIM_DIR/model*/*; do
		PROT_FOLDERS+=($pdb_file)
done


cd ${PROT_FOLDERS[${SLURM_ARRAY_TASK_ID}]}

i=$(basename $PWD)
export GMXLIB=$PARAM_DIR/$i

sed "s/time_input/${time_input}/" "${PARAM_DIR}/${i}/md_diff_sim_time/md_any_ns.mdp" > "${PARAM_DIR}/${i}/md_diff_sim_time/md_${sim_time}ns.mdp"
sed -i "s/sim_time/${sim_time}/" "${PARAM_DIR}/${i}/md_diff_sim_time/md_${sim_time}ns.mdp"


if [[ $i == "CHARMM36M" ]]; then
	pos=SOD
	neg=CLA
else
	pos=NA
	neg=CL
fi


temp_name=(*.pdb)
PROTEIN=${temp_name%.pdb}

: '
if [ -f md*tpr ]; then
       exit 0
fi
'

if [[ $i == "CHARMM36M" || $i == "AMBER99SB-DISP" ]]; then
	sed -i.bak 's/HIP/HIS/g' $temp_name
        echo 1 0 | gmx_mpi pdb2gmx -f ${PROTEIN}.pdb -o ${PROTEIN}.gro -ter -water ${water_model} -ff ${i,,} -ignh
else
        gmx_mpi pdb2gmx -f ${PROTEIN}.pdb -o ${PROTEIN}.gro -water ${water_model} -ff ${i,,} -ignh
fi



gmx_mpi editconf -f ${PROTEIN}.gro -o ${PROTEIN}_newbox.gro -c -d 1.5 -bt dodecahedron

if [[ $i == "AMBER03WS" || $i == AMBER99SBWS ]]; then
	gmx_mpi solvate -cp ${PROTEIN}_newbox.gro -cs ${PARAM_DIR}/${i}/${i,,}.ff/tip4p2005.gro -o ${PROTEIN}_solv.gro -p topol.top
elif [[ $i == "AMBER99SB-DISP" ]]; then
	gmx_mpi solvate -cp ${PROTEIN}_newbox.gro -cs ${PARAM_DIR}/${i}/${i,,}.ff/a99SBdisp_water.gro -o ${PROTEIN}_solv.gro -p topol.top
else
	gmx_mpi solvate -cp ${PROTEIN}_newbox.gro -cs tip4p.gro -o ${PROTEIN}_solv.gro -p topol.top
fi


if [[ $i == "DESAMBER" ]]; then
        sed -i 's/tip4p.itp/tip4pd.itp/' topol.top
elif [[ $i == "AMBER03WS" || $i == AMBER99SBWS ]]; then
        sed -i 's/tip4p.itp/tip4p2005s.itp/' topol.top
elif [[ $i == "AMBER99SB-DISP" ]]; then
        sed -i 's/tip4p.itp/a99SBdisp_water.itp/' topol.top
fi

gmx_mpi grompp -f ${PARAM_DIR}/${i}/ions.mdp -c ${PROTEIN}_solv.gro -p topol.top -o ions.tpr
echo SOL | gmx_mpi -quiet genion -s ions.tpr -o ${PROTEIN}_solv_ions.gro -p topol.top -pname $pos -nname $neg -neutral



srun -n 1 gmx_mpi grompp -f ${PARAM_DIR}/${i}/em.mdp -c ${PROTEIN}_solv_ions.gro -p topol.top -o em_${i}.tpr  
srun gmx_mpi mdrun -deffnm em_${i}


srun -n 1 gmx_mpi grompp -f ${PARAM_DIR}/${i}/nvt.mdp -c em_${i}.gro -r em_${i}.gro -p topol.top -o nvt_${i}.tpr
srun gmx_mpi mdrun -deffnm nvt_${i}



srun -n 1 gmx_mpi grompp -f ${PARAM_DIR}/${i}/npt.mdp -c nvt_${i}.gro -r nvt_${i}.gro -t nvt_${i}.cpt -p topol.top -o npt_${i}.tpr
srun gmx_mpi mdrun -deffnm npt_${i}

srun gmx_mpi grompp -f ${PARAM_DIR}/${i}/md_diff_sim_time/md_${sim_time}ns.mdp -c npt_${i}.gro -t npt_${i}.cpt -p topol.top -o md_${sim_time}ns.tpr



srun gmx_mpi grompp -f ${PARAM_DIR}/${i}/md_diff_sim_time/md_${sim_time}ns.mdp -c npt_${i}.gro -t npt_${i}.cpt -p topol.top -o md_${sim_time}ns.tpr

