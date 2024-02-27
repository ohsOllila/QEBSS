#!/bin/bash


time_input=1500

cd ..
cd ..
BASE_DIR=$PWD
SIM_DIR=$BASE_DIR/Unst*

SCRIPTS=$BASE_DIR/simulation_scripts/MD_scripts
md_script=${SCRIPTS}/md.sh
cp ${md_script} ${SCRIPTS}/batch_md.sh
JOB_SCRIPT=${SCRIPTS}/batch_md.sh

sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"

for i in $SIM_DIR/model*/*; do
	cd $i
	if [ ! -f md*$time_input*gro ]; then 
		sbatch ${JOB_SCRIPT}
	fi
done

