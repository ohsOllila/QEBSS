#!/bin/bash

cd ..
cd ..
SIM_DIR=${PWD}
batch=${SIM_DIR}/simulation_scripts/MD_scripts/em.sh

FORCEFIELDS=(AMBER03WS AMBER99SB-DISP AMBER99SBWS CHARMM36 DESAMBER)


for pdb_file in *pdb; do
	pdb_filename=$(basename $pdb_file)
	pdb_name=${pdb_filename%.pdb}
	pdb_folder=$SIM_DIR/$pdb_name
	# Create the folder if it doesn't exist
    
	# Move the PDB file into the folder
	cd $pdb_folder
	for i in "${FORCEFIELDS[@]}"; do
    		cd $pdb_folder/$i
    		sbatch $batch
	done
done
