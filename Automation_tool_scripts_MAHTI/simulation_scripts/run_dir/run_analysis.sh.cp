#!/bin/bash

time_input=1000

cd ..
cd ..
SIM_DIR=${PWD}

SCRIPTS=${SIM_DIR}/simulation_scripts/MD_scripts
analysis_script=${SCRIPTS}/analysis.sh
cp ${analysis_script} ${SCRIPTS}/batch_md.sh


JOB_SCRIPT=${SCRIPTS}/batch_md.sh

sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"


FORCEFIELDS=(AMBER03WS AMBER99SB-DISP AMBER99SBWS CHARMM36M DESAMBER)


for pdb_file in Unst_snear/*pdb; do
        pdb_filename=$(basename $pdb_file)
        pdb_name=${pdb_filename%.pdb}
        pdb_folder=$SIM_DIR/Unst_snear/$pdb_name
        # Create the folder if it doesn't exist
    
        # Move the PDB file into the folder
        cd $pdb_folder
        for i in "${FORCEFIELDS[@]}"; do
                cd $pdb_folder/$i
		sbatch ${JOB_SCRIPT}
        done
done
