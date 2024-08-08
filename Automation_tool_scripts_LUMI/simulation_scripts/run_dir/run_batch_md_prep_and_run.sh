#!/bin/bash

time_input=2000


cd ..
cd ..
BASE_DIR=${PWD}
SCRIPTS=${BASE_DIR}/simulation_scripts/MD_scripts
md_script=${SCRIPTS}/md_prep_and_run.sh


for pdb_file in "$BASE_DIR"/*.pdb; do
    model="$(basename "${pdb_file%.pdb}")"
    mkdir -p "$SIM_DIR/AMBER03WS/$model"
done


your_projects=$(csc-projects | grep -o "project_.*" | awk '{print $1}')
echo "Select the number of the project you want to use:"

num=1
list=()

for i in $your_projects; do
        list+=($i)
        echo "("$num")" ${i} 
        ((num++))
done

read choice
project=${list[choice-1]}




for i in $SIM_DIR/AMBER03WS/$model/
do
  	cd $i

	cp ${md_script} ${SCRIPTS}/batch_md.sh	
	JOB_SCRIPT=${SCRIPTS}/batch_md.sh

	sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"
	sed -i "s/project/${project}/" "${JOB_SCRIPT}"

	sbatch ${JOB_SCRIPT}
done


