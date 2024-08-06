#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=50000
#SBATCH --array=0-24
#SBATCH --account=project_462000540
##SBATCH --account=project


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GMX_MAXBACKUP=-1


module use /appl/local/csc/modulefiles
module load gromacs/2023.3-gpu
sim_time=sim_time

SIM_PATH=${PWD}
SIM_DIR=$(basename $SIM_PATH)

BASE_DIR=$(cd .. && pwd)

magn_field=$(awk 'NR==1 {print $6}' "${BASE_DIR}/${SIM_DIR}_exp_data.txt" 2>/dev/null)
make_index=${BASE_DIR}/simulation_scripts/MD_scripts/makeNHindex.awk
py_script=${BASE_DIR}/simulation_scripts/PY_scripts/Old_Relaxations_for_Samuli.py
py_script2=${BASE_DIR}/simulation_scripts/MD_scripts/relaxation_times/relaxation_times.py
mdmat_plot=${BASE_DIR}/simulation_scripts/PY_scripts/xpm_plot.py
secondary=${BASE_DIR}/simulation_scripts/PY_scripts/pymol_structure_analysis.py
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
cp $py_script2 ${path}

TEMP_NAME=(md_${sim_time}ns.tpr)
name=${TEMP_NAME%.tpr}

gmx_mpi check -f ${name}.xtc

mkdir correlation_functions

#sed -i.bak 's/ H /HN /g' ${name}.gro
#awk -f ${make_index} ${name}.gro > HN.ndx


echo 1 1 | gmx_mpi trjconv -f ${name}.xtc -s ${name}.tpr -pbc mol -center -dump 0 -o temp_${name}.gro
echo 1 1 | gmx_mpi trjconv -f ${name}.xtc -s ${name}.tpr -pbc mol -center -o ${name}_noPBC.xtc
echo 1 | gmx_mpi trjconv -f ${name}_noPBC.xtc -s temp_${name}.gro -skip 10 -o ${name}_skip.xtc

if [ -e ${name}_gyrate.xvg ]
then
    echo "Radius of gyration already calculated"
else
    echo 1 | gmx_mpi gyrate -s ${name}.tpr -f ${name}_noPBC.xtc -o ${name}_gyrate.xvg
fi
    
if [ -e ${name}_mdmat.eps ]
then
    echo "MDmat Already calculated"
else
    echo -e "Alpha\nAlpha" | gmx_mpi mdmat -f ${name}_skip.xtc -s ${name}.tpr -mean ${name}_mdmat.xpm
    gmx_mpi xpm2ps -f ${name}_mdmat.xpm -o ${name}_mdmat.eps
fi


GRO_FILE=(temp_md_${sim_time}ns.gro)
sed -i.bak 's/ H /HN /g' $GRO_FILE
sed -i.bak 's/H1/HN/g' $GRO_FILE
awk -f ${make_index} $GRO_FILE > HN.ndx

line_number=1  # Initialize the line number

numberOFfuncs=$(grep "\[" HN.ndx | tail -n 1 | awk '{print $2}')
for ((i = 0; i <= $numberOFfuncs; i++)); do
    num=$(awk -v line="$line_number" 'NR==line {print $2}' HN.ndx)
    if [ -e correlation_functions/NHrotaCF_$num.xvg ]
    then
	echo "Correlation function already calculated"
    else
	echo $i | gmx_mpi rotacf -f ${name}.xtc -s ${name}.tpr -n HN.ndx -o correlation_functions/NHrotaCF_$num.xvg -P 2 -d -xvg none  #-nice 20 
	((line_number += 2))
    fi
done

module purge
export PATH="$(cd ../../../env/bin && pwd):$PATH"
#export PATH="/scratch/project_462000285/cmcajsa/systems/forcefield_compare/env/bin:$PATH"

python3 $mdmat_plot
python3 $secondary
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
python3 $relax_plot

