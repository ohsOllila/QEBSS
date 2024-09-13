#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=128
#SBATCH --nodes=1
#SBATCH --array=0-num_jobs
#SBATCH --account=project_2003809
##SBATCH --account=project


#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GMX_MAXBACKUP=-1


module load gromacs-env

sim_time=sim_time

SIM_PATH=${PWD}
SIM_DIR=$(basename $SIM_PATH)

BASE_DIR=$(cd .. && pwd)

magn_field=$(awk 'NR==1 {print $6}' "${BASE_DIR}/${SIM_DIR}_exp_data.txt" 2>/dev/null)
make_index=${BASE_DIR}/simulation_scripts/MD_scripts/makeNHindex.awk
py_script=${BASE_DIR}/simulation_scripts/PY_scripts/Old_Relaxations_for_Samuli.py
contact_plot=${BASE_DIR}/simulation_scripts/PY_scripts/create_contact.py
dist_plot=${BASE_DIR}/simulation_scripts/PY_scripts/create_distance.py
#secondary=${BASE_DIR}/simulation_scripts/PY_scripts/pymol_structure_analysis.py
relax_plot=${BASE_DIR}/simulation_scripts/PY_scripts/plot_replicas_to_experiment.py
corr_plot=${BASE_DIR}/simulation_scripts/PY_scripts/correlationCALC.py


mkdir -p $BASE_DIR/results/${SIM_DIR}

list=(${SIM_PATH}/model*/*/)


cd ${list[${SLURM_ARRAY_TASK_ID}]}



path=${PWD}

base_dir=$(basename "$(dirname "$path")")
replicas=$(basename "$base_dir")
ff=$(basename $path)

cp $py_script ${path}
sed -i "s|magn_field=magn_field|magn_field=$magn_field|" ${path}/Old_Relaxations_for_Samuli.py

TEMP_NAME=(md_${sim_time}ns.tpr)
name=${TEMP_NAME%.tpr}

gmx_mpi check -f ${name}.xtc

mkdir correlation_functions


echo 1 1 | gmx_mpi trjconv -f ${name}.xtc -s ${name}.tpr -pbc mol -center -dump 0 -o temp_${name}.gro
echo 1 1 | gmx_mpi trjconv -f ${name}.xtc -s ${name}.tpr -pbc mol -center -o ${name}_noPBC.xtc
echo 1 | gmx_mpi gyrate -s ${name}.tpr -f ${name}_noPBC.xtc -o ${name}_gyrate.xvg

#echo -e "Alpha\nAlpha" | gmx_mpi mdmat -f ${name}_noPBC.xtc -s ${name}.tpr -mean ${name}_mdmat.xpm
#gmx_mpi xpm2ps -f ${name}_mdmat.xpm -o ${name}_mdmat.eps


GRO_FILE=(temp_md_${sim_time}ns.gro)
sed -i.bak 's/ H /HN /g' $GRO_FILE
sed -i.bak 's/H1/HN/g' $GRO_FILE
awk -f ${make_index} $GRO_FILE > HN.ndx

line_number=1  # Initialize the line number

numberOFfuncs=$(awk -v lines="$(wc -l < HN.ndx)" 'BEGIN {print int(lines / 2)}')
for ((i = 0; i <= $numberOFfuncs; i++)); do
	num=$(awk -v line="$line_number" 'NR==line {print $2}' HN.ndx)
	echo $i | gmx_mpi rotacf -f ${name}_noPBC.xtc -s ${name}.tpr -n HN.ndx -o correlation_functions/NHrotaCF_$num.xvg -P 2 -d -xvg none  #-nice 20 
	((line_number += 2)) 
done

module purge
export PATH="$(cd ../../../env/bin && pwd):$PATH"
#export PATH="/scratch/project_2003809/cmcajsa/forcefield/env/bin:$PATH"

python3 $contact_plot
python3 $dist_plot
#python3 $secondary
python3 $corr_plot
python3 ${path}/Old_Relaxations_for_Samuli.py


mkdir -p $BASE_DIR/results/${SIM_DIR}/$replicas/$ff
sim_results=$BASE_DIR/results/${SIM_DIR}/$replicas/$ff
cp $path/*gyrate*xvg $sim_results
cp $path/*mdmat* $sim_results
cp $path/*png $sim_results
cp $path/*eps  $sim_results
cp -r $path/correlation_functions $sim_results

cd $SIM_PATH


if [ "$SLURM_ARRAY_TASK_ID" -eq "$SLURM_ARRAY_TASK_MAX" ]; then
	python3 $relax_plot
fi
