#!/bin/bash

cd ..
SIM_SCRIPTS=${PWD}

cd ..
SIM_DIR=${PWD}

mkdir results

cd results
ions=(CA MG FE2P ZN CU2P MN2P NA K)


for i in $SIM_DIR/H*; do
	protein_path=${SIM_DIR}/${i}
	protein=$(basename ${protein_path})
	mkdir ${protein}
	cd ${protein}
	for ion in "${ions[@]}"; do
		mkdir $ion
	done
	# Loop over the ions and copy all files from the corresponding subdirectory to the protein's directory
	for ion in "${ions[@]}"; do
    		cp ${SIM_DIR}/${protein}/${ion}/md*xtc* ${SIM_DIR}/results/${protein}/${ion}/
		cp ${SIM_DIR}/${protein}/${ion}/temp*gro ${SIM_DIR}/results/${protein}/${ion}/
		cp ${SIM_DIR}/${protein}/${ion}/*rmsd* ${SIM_DIR}/results/${protein}/${ion}/
  	done
	cd ..
done

