#!/bin/bash

time_input=2000

cd ..
cd ..
BASE_DIR=$PWD

SIM_DIR=$BASE_DIR/Unst*

#your_projects=$(csc-projects | grep -o "project_.*" | awk '{print $1}')
#echo "Select the number of the project you want to use:"

num=1
list=()

#for i in $your_projects; do
#        list+=($i)
#        echo "("$num")" ${i} 
#        ((num++))
#done

#read choice
#project=${list[choice-1]}


SCRIPTS=$BASE_DIR/simulation_scripts/MD_scripts
md_script=${SCRIPTS}/md.sh
cp ${md_script} ${SCRIPTS}/batch_md.sh
JOB_SCRIPT=${SCRIPTS}/batch_md.sh

sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"
#sed -i "s/project/${project}/" "${JOB_SCRIPT}"

for i in $SIM_DIR/model*/*; do
	cd $i
	if [ ! -f md*$time_input*gro ]; then 
		sbatch ${JOB_SCRIPT}
	fi
done
