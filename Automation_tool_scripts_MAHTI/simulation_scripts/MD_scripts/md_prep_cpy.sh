#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --ntasks=1
##SBATCH --mem-per-cpu=500
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=1
#SBATCH --array=0-num_jobs
#SBATCH --account=project_2003809
##SBATCH --account=project

#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GMX_MAXBACKUP=-1


module load gromacs-env

sim_time=sim_time
water_model=tip4p

SIM_DIR=${PWD}
PARAM_DIR=$(cd ../../MD_parameter_files && pwd)
time_input=$((500000 * $sim_time))


FORCEFIELD=("AMBER03WS" "AMBER99SB-DISP" "AMBER99SBWS" "CHARMM36M" "DESAMBER")
#FORCEFIELD=("CHARMM36M")

PROT_FOLDERS=()
for pdb_file in $SIM_DIR/replica*/*; do
		PROT_FOLDERS+=($pdb_file)
done

cd ${PROT_FOLDERS[${SLURM_ARRAY_TASK_ID}]}

JOB_FOLD=${PWD}
mdp=$(basename $PWD)
echo $JOB_FOLD
export GMXLIB=$PARAM_DIR/$mdp


sed "s/time_input/${time_input}/" "${PARAM_DIR}/$mdp/md_diff_sim_time/md_any_ns.mdp" > "${PARAM_DIR}/$mdp/md_diff_sim_time/md_${sim_time}ns.mdp"
sed -i "s/sim_time/${sim_time}/" "${PARAM_DIR}/$mdp/md_diff_sim_time/md_${sim_time}ns.mdp"


if [[ $mdp == "CHARMM36M" ]]; then
	pos=SOD
	neg=CLA
else
	pos=NA
	neg=CL
fi


temp_name=(rep*.pdb)
PROTEIN=${temp_name%.pdb}



# Define an array for box types and distances
box_types=("triclinic" "dodecahedron")
#box_types=("triclinic")
distances=("0.1" "0.2" "0.5" "0.8")  # You can adjust these as needed
#distances=("0.8")

i=$(basename $(dirname $PWD))_$(basename $PWD)
# Loop through each combination of box type and distance
for box_type in "${box_types[@]}"; do
    for dist in "${distances[@]}"; do
        mkdir $JOB_FOLD/md_${i}_${box_type}_${dist}
        cp $JOB_FOLD/$temp_name $JOB_FOLD/md_${i}_${box_type}_${dist}
        cd $JOB_FOLD/md_${i}_${box_type}_${dist} 
        # Check if the simulation already exists for the given box type and distance
        if [ -f md_${i}_${box_type}_${dist}*tpr ]; then
            echo "Simulation for box type $box_type and distance $dist already exists, skipping..."
            continue  # Skip this iteration if file exists
        fi

        # Prepare the system with pdb2gmx and editconf
        echo "Running simulation with box type $box_type and distance $dist..."

        if [[ $mdp == "CHARMM36M" || $mdp == "AMBER99SB-DISP" ]]; then
            sed -i.bak 's/HIP/HIS/g' $temp_name
            echo 0 1 | gmx_mpi pdb2gmx -f ${PROTEIN}.pdb -o ${PROTEIN}.gro -ter -water ${water_model} -ff ${mdp,,} -ignh
        else
            gmx_mpi pdb2gmx -f ${PROTEIN}.pdb -o ${PROTEIN}.gro -water ${water_model} -ff ${mdp,,} -ignh
        fi

        # Adjust the box type and distance
        gmx_mpi editconf -f ${PROTEIN}.gro -o ${PROTEIN}_newbox.gro -c -d $dist -bt $box_type


	if [[ $mdp == "AMBER03WS" || $mdp == AMBER99SBWS ]]; then
		gmx_mpi solvate -cp ${PROTEIN}_newbox.gro -cs ${PARAM_DIR}/${mdp}/${mdp,,}.ff/tip4p2005.gro -o ${PROTEIN}_solv_${box_type}_${dist}.gro -p topol.top
	elif [[ $mdp == "AMBER99SB-DISP" ]]; then
		gmx_mpi solvate -cp ${PROTEIN}_newbox.gro -cs ${PARAM_DIR}/${mdp}/${mdp,,}.ff/a99SBdisp_water.gro -o ${PROTEIN}_solv_${box_type}_${dist}.gro -p topol.top
	else
		gmx_mpi solvate -cp ${PROTEIN}_newbox.gro -cs tip4p.gro -o ${PROTEIN}_solv_${box_type}_${dist}.gro -p topol.top
	fi


	if [[ $mdp == "DESAMBER" ]]; then
        	sed -i 's/tip4p.itp/tip4pd.itp/' topol.top
	elif [[ $mdp == "AMBER03WS" || $mdp == AMBER99SBWS ]]; then
        	sed -i 's/tip4p.itp/tip4p2005s.itp/' topol.top
	elif [[ $mdp == "AMBER99SB-DISP" ]]; then
        	sed -i 's/tip4p.itp/a99SBdisp_water.itp/' topol.top
	fi


        # Add ions
        gmx_mpi grompp -f ${PARAM_DIR}/$mdp/ions.mdp -c ${PROTEIN}_solv_${box_type}_${dist}.gro -p topol.top -o ions_${box_type}_${dist}.tpr
        echo SOL | gmx_mpi -quiet genion -s ions_${box_type}_${dist}.tpr -o ${PROTEIN}_solv_ions_${box_type}_${dist}.gro -p topol.top -pname $pos -nname $neg -neutral

        # Energy minimization
        srun -n 1 gmx_mpi grompp -f ${PARAM_DIR}/$mdp/em.mdp -c ${PROTEIN}_solv_ions_${box_type}_${dist}.gro -p topol.top -o em_${i}_${box_type}_${dist}.tpr
        [ ! -f "em_${i}_${box_type}_${dist}.gro" ] || srun gmx_mpi mdrun -deffnm em_${i}_${box_type}_${dist} -cpi em_${i}_${box_type}_${dist}.cpt


	srun -n 1 gmx_mpi grompp -f ${PARAM_DIR}/$mdp/nvt.mdp -c em_${i}_${box_type}_${dist}.gro -r em_${i}_${box_type}_${dist}.gro -p topol.top -o nvt_${i}_${box_type}_${dist}.tpr
	[ ! -f "nvt_${i}_${box_type}_${dist}.gro" ] || srun gmx_mpi mdrun -deffnm nvt_${i}_${box_type}_${dist} -cpi nvt_${i}_${box_type}_${dist}.cpt



	srun -n 1 gmx_mpi grompp -f ${PARAM_DIR}/$mdp/npt.mdp -c nvt_${i}_${box_type}_${dist}.gro -r nvt_${i}_${box_type}_${dist}.gro -t nvt_${i}_${box_type}_${dist}.cpt -p topol.top -o npt_${i}_${box_type}_${dist}.tpr
	[ ! -f "npt_${i}_${box_type}_${dist}.gro" ] || srun gmx_mpi mdrun -deffnm npt_${i}_${box_type}_${dist} -cpi npt_${i}_${box_type}_${dist}.cpt



	srun gmx_mpi grompp -f ${PARAM_DIR}/$mdp/md_diff_sim_time/md_${sim_time}ns.mdp -c npt_${i}_${box_type}_${dist}.gro -t npt_${i}_${box_type}_${dist}.cpt -p topol.top -o md_${i}_${box_type}_${dist}.tpr
        cp ${PARAM_DIR}/$mdp/md_diff_sim_time/md_${sim_time}ns.mdp .
    done
done

# Wait for all background jobs to finish before continuing
wait

