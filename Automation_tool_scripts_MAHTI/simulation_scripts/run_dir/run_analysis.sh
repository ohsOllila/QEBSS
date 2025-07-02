#!/bin/bash



cd ..
cd ..
BASE_DIR=${PWD}
SCRIPTS=${BASE_DIR}/simulation_scripts/MD_scripts
analysis_script=${SCRIPTS}/analysis.sh


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

echo "Input simulation time in ns:"
read time_input


for i in $BASE_DIR/Unst_prot/
do
  	cd $i
	jobs=$(( $(find "$i" -mindepth 2 -maxdepth 2 -type d | wc -l) - 1 ))

	cp ${analysis_script} ${SCRIPTS}/batch_analysis.sh
	JOB_SCRIPT=${SCRIPTS}/batch_analysis.sh
	sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"
	sed -i "s/num_jobs/${jobs}/" "${JOB_SCRIPT}"
	sed -i "s/project/${project}/" "${JOB_SCRIPT}"

	sbatch ${JOB_SCRIPT}
done
