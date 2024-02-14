#!/bin/bash



time_input=1000


cd ..
cd ..
BASE_DIR=${PWD}

SCRIPTS=${BASE_DIR}/simulation_scripts/MD_scripts
analysis_script=${SCRIPTS}/analysis.sh
cp ${analysis_script} ${SCRIPTS}/batch_analysis.sh

jobs=$(( $(find "$BASE_DIR/Unst_prot" -mindepth 2 -maxdepth 2 -type d | wc -l) - 1 ))

JOB_SCRIPT=${SCRIPTS}/batch_analysis.sh

sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"
sed -i "s/num_jobs/${jobs}/" "${JOB_SCRIPT}"

for i in $BASE_DIR/Unst_prot
do
  	cd $i
	sbatch ${JOB_SCRIPT}
done

