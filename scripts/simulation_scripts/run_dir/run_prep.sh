#!/bin/bash


echo "Input simulation time"
read time_input


cd ..
cd ..
BASE_DIR=${PWD}

SCRIPTS=${BASE_DIR}/simulation_scripts/MD_scripts
md_script=${SCRIPTS}/md_prep.sh
cp ${md_script} ${SCRIPTS}/batch_md.sh


JOB_SCRIPT=${SCRIPTS}/batch_md.sh

FORCEFIELD=(AMBER03WS AMBER99SB-DISP AMBER99SBWS CHARMM36M DESAMBER)


for pdb_file in $BASE_DIR/Unst_prot/*.pdb; do
	replicas=$(basename ${pdb_file%.pdb})
	for i in "${FORCEFIELD[@]}"; do
        	mkdir -p $BASE_DIR/Unst_prot/$replicas/${i}
        	cp -R -u -p $pdb_file $BASE_DIR/Unst_prot/$replicas/${i}
	done
done

list=$BASE_DIR/Unst_prot/*/*/
jobs=$(( $(ls -l $list | grep -c '^d') - 1 ))


sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"
sed -i "s/num_jobs/${jobs}/" "${JOB_SCRIPT}"

for i in $BASE_DIR/Unst_prot
do
  	cd $i
	sh ${JOB_SCRIPT}
done


