#!/bin/bash

cd ..
SIM_SCRIPTS=${PWD}

cd ..
SIM_DIR=${PWD}

mkdir results > /dev/null 2>&1
cd results

FORCEFIELD=(AMBER03WS AMBER99SB-DISP AMBER99SBWS CHARMM36M DESAMBER)

for replica_path in $SIM_DIR/U*/*/; do
	protein=$(basename $(dirname ${replica_path}))
	replica=$(basename ${replica_path})
	
	
	mkdir -p $SIM_DIR/results/${protein}/${replica} > /dev/null 2>&1
	cd $SIM_DIR/results/${protein}/${replica}
	for f in "${FORCEFIELD[@]}"; do
		mkdir $f > /dev/null 2>&1
	done
	# Loop over the ions and copy all files from the corresponding subdirectory to the protein's directory

	for f in "${FORCEFIELD[@]}"; do
	#	cp ${replica_path}/${f}/temp*gro ${SIM_DIR}/results/${protein}/${replica}/${f}/ > /dev/null 2>&1
		cp ${replica_path}/${f}/*noPBC* ${SIM_DIR}/results/${protein}/${replica}/${f}/ > /dev/null 2>&1
    	#	cp ${replica_path}/${f}/*smooth* ${SIM_DIR}/results/${protein}/${replica}/${f}/ > /dev/null 2>&1
		cp ${replica_path}/${f}/md_1000ns.tpr ${SIM_DIR}/results/${protein}/${replica}/${f}/ > /dev/null 2>&1
		cp ${replica_path}/${f}/md_1000ns_noPBC.xtc ${SIM_DIR}/results/${protein}/${replica}/${f}/ > /dev/null 2>&1
		cp ${replica_path}/${f}/HN.ndx ${SIM_DIR}/results/${protein}/${replica}/${f}/ > /dev/null 2>&1
		cp ${replica_path}/${f}/*xvg ${SIM_DIR}/results/${protein}/${replica}/${f}/ > /dev/null 2>&1
		cp -r ${replica_path}/${f}/correlation_functions ${SIM_DIR}/results/${protein}/${replica}/${f}/ > /dev/null 2>&1
  	done

done
