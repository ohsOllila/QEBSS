#!/bin/bash

time_input=2000
PROT_NAME="EB1"

cd ..
cd ..
BASE_DIR=${PWD}
SCRIPTS=${BASE_DIR}/simulation_scripts/MD_scripts
md_script=${SCRIPTS}/md_prep.sh

: '
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
'


FORCEFIELD=(AMBER03WS AMBER99SB-DISP AMBER99SBWS CHARMM36M DESAMBER)
#FORCEFIELD=(AMBER03WS)

for pdb_file in $BASE_DIR/TEST/$PROT_NAME/rep*.pdb; do
	directory_path="${pdb_file%/*}"
	replicas=$(basename ${pdb_file%.pdb})
	for i in "${FORCEFIELD[@]}"; do
        	mkdir -p $directory_path/$replicas/${i}
        	cp -R -u -p $pdb_file $directory_path/$replicas/${i}
	done
done


for i in $BASE_DIR/TEST/$PROT_NAME/
do
  	cd $i
	jobs=$(( $(find $i -mindepth 2 -maxdepth 2 -type d | wc -l) - 1 ))

	#jobs=1
	cp ${md_script} ${SCRIPTS}/batch_md.sh	
	JOB_SCRIPT=${SCRIPTS}/batch_md.sh
	sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"
	sed -i "s/num_jobs/${jobs}/" "${JOB_SCRIPT}"
	#sed -i "s/project/${project}/" "${JOB_SCRIPT}"

	sbatch ${JOB_SCRIPT}
done


