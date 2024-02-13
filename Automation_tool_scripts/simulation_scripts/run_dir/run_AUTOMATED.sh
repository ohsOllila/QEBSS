#!/bin/bash


echo "Input simulation time"
read time_input

RUN_DIR=${PWD}

cd ..
SIM_SCRIPTS=${PWD}

cd ..
SIM_DIR=${PWD}
mkdir output_files

md_script1=${SIM_SCRIPTS}/MD_scripts/md_prep.sh
cp ${md_script1} ${SIM_SCRIPTS}/MD_scripts/batch_prep.sh

md_script2=${SIM_SCRIPTS}/MD_scripts/md.sh
cp ${md_script2} ${SIM_SCRIPTS}/MD_scripts/batch_run.sh

md_script3=${SIM_SCRIPTS}/MD_scripts/analysis.sh
cp ${md_script3} ${SIM_SCRIPTS}/MD_scripts/batch_analysis.sh

SCRIPT1=${SIM_SCRIPTS}/MD_scripts/batch_prep.sh
SCRIPT2=${SIM_SCRIPTS}/MD_scripts/batch_run.sh
SCRIPT3=${SIM_SCRIPTS}/MD_scripts/batch_analysis.sh

list=$SIM_DIR/Unst_hydrolase/*/*/
seq_dir=($(find $SIM_DIR -maxdepth 2 -type d -name "Unst_hydrolase"))

jobs=$(ls -l $list | grep -c '^d')

sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${SCRIPT1}" "${SCRIPT2}" "${SCRIPT3}"
sed -i "s/num_jobs/${jobs}/" "${SCRIPT1}" "${SCRIPT3}"


cd $seq_dir
sbatch ${SCRIPT1}

for i in $list
do
	cd $i
	sbatch --dependency=afterany:$(squeue -h -o %i -n batch_prep.sh) ${SCRIPT2}
	while [ ! -f md*$time_input*gro ] && [ ! squeue -h -o %i -n batch_run.sh ]; do
        	sbatch ${SCRIPT2}
		sleep 300
	done

done	

cd $seq_dir
while true; do
	if [ ! -f md*$time_input*gro ]; then
		sleep 300
	else            
		sbatch ${SCRIPT3}
		break
	fi
done
