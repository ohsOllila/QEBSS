
#!/bin/bash



time_input=1500


cd ..
cd ..
BASE_DIR=${PWD}

SCRIPTS=${BASE_DIR}/simulation_scripts/MD_scripts
analysis_script=${SCRIPTS}/analysis.sh
cp ${analysis_script} ${SCRIPTS}/batch_analysis.sh


JOB_SCRIPT=${SCRIPTS}/batch_analysis.sh

sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"

for i in $BASE_DIR/Unst*
do
  	cd $i
	directory_path=$(basename $i)
	jobs=$(( $(find "$i" -mindepth 2 -maxdepth 2 -type d | wc -l) - 1 ))
	sed -i "s/num_jobs/${jobs}/" "${JOB_SCRIPT}"
	sbatch ${JOB_SCRIPT}
done

