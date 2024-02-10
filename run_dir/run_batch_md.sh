#!/bin/bash


time_input=1000

cd ..
cd ..
BASE_DIR=$PWD
SIM_DIR=$BASE_DIR/Unst*

SCRIPTS=$BASE_DIR/simulation_scripts/MD_scripts
md_script=${SCRIPTS}/md.sh
cp ${md_script} ${SCRIPTS}/batch_md.sh


JOB_SCRIPT=${SCRIPTS}/batch_md.sh

list=$SIM_DIR/*/*/
jobs=$(( $(ls -l $list | grep -c '^d') - 1 ))

sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"
sed -i "s/num_jobs/${jobs}/" "${JOB_SCRIPT}"

for i in $SIM_DIR/model*/*; do
	cd $i
	if [ ! -f md*1000ns*gro ]; then 
		sbatch ${JOB_SCRIPT}
	fi
done

