#!/bin/bash

cd ..
cd ..
SIM_DIR=${PWD}
mkdir results > /dev/null 2>&1
cd results


FORCEFIELDS=(AMBER03WS AMBER99SB-DISP AMBER99SBWS CHARMM36 DESAMBER)


for pdb_file in $SIM_DIR/*pdb; do
        pdb_filename=$(basename $pdb_file)
        pdb_name=${pdb_filename%.pdb}
        pdb_folder=$SIM_DIR/$pdb_name
        # Create the folder if it doesn't exist
	mkdir -p $SIM_DIR/results/$pdb_name > /dev/null 2>&1    
        # Move the PDB file into the folder
        cd $pdb_folder
        for i in "${FORCEFIELDS[@]}"; do
		mkdir $SIM_DIR/results/$pdb_name/$i > /dev/null 2>&1
                cp $pdb_folder/$i/*pdb $SIM_DIR/results/$pdb_name/$i
	#	cp $pdb_folder/$i/*xtc $SIM_DIR/results/$pdb_name/$i
      		cp $pdb_folder/$i/*xvg $SIM_DIR/results/$pdb_name/$i
	  done
done


