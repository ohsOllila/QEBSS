#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --partition=medium
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --array=0-num_jobs
#SBATCH --account=Project_2001058	
##SBATCH --mail-type=END #uncomment to get mail



export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GMX_MAXBACKUP=-1


module purge
module load gromacs-env


sim_time=sim_time

SIM_PATH=${PWD}
SIM_DIR=$(basename $SIM_PATH)

BASE_DIR=$(cd .. && pwd)

make_index=${BASE_DIR}/simulation_scripts/MD_scripts/makeNHindex.awk
py_script=${BASE_DIR}/simulation_scripts/PY_scripts/Old_Relaxations_for_Samuli.py
mdmat_plot=${BASE_DIR}/simulation_scripts/PY_scripts/xpm_plot.py
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
sed -i "s|PATH_TO_CORR|${path}/correlation_functions|" ${path}/Old_Relaxations_for_Samuli.py

TEMP_NAME=(md_${sim_time}ns.tpr)
name=${TEMP_NAME%.tpr}


mkdir correlation_functions

#sed -i.bak 's/ H /HN /g' ${name}.gro
#awk -f ${make_index} ${name}.gro > HN.ndx

echo 1 1 | gmx_mpi trjconv -f ${name}.xtc -s ${name}.tpr -pbc mol -center -dump 0 -o temp_${name}.gro
echo 1 1 | gmx_mpi trjconv -f ${name}.xtc -s ${name}.tpr -pbc mol -center -o ${name}_noPBC.xtc
gmx_mpi filter -f ${name}_noPBC.xtc -s temp_${name}.gro -nf 20 -all -ol ${name}_smooth.xtc
echo 1 | gmx_mpi gyrate -s ${name}.tpr -f ${name}_noPBC.xtc -o ${name}_gyrate.xvg

echo -e "Alpha\nAlpha" | gmx_mpi mdmat -f ${name}.xtc -s ${name}.tpr -mean ${name}_mdmat.xpm
gmx_mpi xpm2ps -f ${name}_mdmat.xpm -o ${name}_mdmat.eps

GRO_FILE=(temp_md_1000ns.gro)
sed -i.bak 's/ H /HN /g' $GRO_FILE
sed -i.bak 's/H1/HN/g' $GRO_FILE
sed -i.bak 's/CD/HN/g' $GRO_FILE
awk -f ${make_index} $GRO_FILE > HN.ndx
numberOFfuncs=$(grep "\[" HN.ndx | tail -n 1 | awk '{print $2}')
for ((i = 0; i <= $numberOFfuncs; i++))
do
	if [ ! -s correlation_functions/NHrotaCF_$i.xvg ]; then

		echo $i | gmx_mpi rotacf -f ${name}_noPBC.xtc -s ${name}.tpr -n HN.ndx -o correlation_functions/NHrotaCF_$i.xvg -P 2 -d -xvg none  #-nice 20 &
	fi
done
	

module purge
export PATH="$(cd ../../../env/bin && pwd):$PATH"

python3 $mdmat_plot
python3 $corr_plot
python3 ${path}/Old_Relaxations_for_Samuli.py > relaxation_data.txt
python3 $relax_plot

module purge

mkdir -p $BASE_DIR/results/${SIM_DIR}/$replicas/$ff
sim_results=$BASE_DIR/results/${SIM_DIR}/$replicas/$ff
cp $path/*gyrate*xvg $sim_results
cp $path/*mdmat* $sim_results
cp $path/*png $sim_results
cp $path/*eps  $sim_results
cp -r $path/correlation_functions $sim_results




