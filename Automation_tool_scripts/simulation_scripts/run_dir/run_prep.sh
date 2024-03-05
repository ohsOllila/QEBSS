#!/bin/bash


time_input=1500


cd ..
cd ..
BASE_DIR=${PWD}

SCRIPTS=${BASE_DIR}/simulation_scripts/MD_scripts
md_script=${SCRIPTS}/md_prep.sh
cp ${md_script} ${SCRIPTS}/batch_md.sh


JOB_SCRIPT=${SCRIPTS}/batch_md.sh

FORCEFIELD=(AMBER03WS AMBER99SB-DISP AMBER99SBWS CHARMM36M DESAMBER)


for pdb_file in $BASE_DIR/Unst*/*.pdb; do
	directory_path="${pdb_file%/*}"
	replicas=$(basename ${pdb_file%.pdb})
	for i in "${FORCEFIELD[@]}"; do
        	mkdir -p $directory_path/$replicas/${i}
        	cp -R -u -p $pdb_file $directory_path/$replicas/${i}
	done
done




for i in $BASE_DIR/Unst*
do
  	cd $i
	jobs=$(( $(find $i -mindepth 2 -maxdepth 2 -type d | wc -l) - 1 ))
	sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"
	sed -i "s/num_jobs/${jobs}/" "${JOB_SCRIPT}"
	sbatch ${JOB_SCRIPT}
done



