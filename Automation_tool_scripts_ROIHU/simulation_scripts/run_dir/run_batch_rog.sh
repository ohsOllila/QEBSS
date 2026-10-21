#!/bin/bash

echo "Input simulation time"
read time_input


cd ..
cd ..
SIM_DIR=${PWD}

SCRIPTS=${SIM_DIR}/simulation_scripts/MD_scripts
rog_script=${SCRIPTS}/rog.sh
cp ${rog_script} ${SCRIPTS}/batch_rog.sh


JOB_SCRIPT=${SCRIPTS}/batch_rog.sh

sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"


FORCEFIELDS=(AMBER03WS AMBER99SB-DISP AMBER99SBWS CHARMM36M DESAMBER)


for pdb_file in U*/*pdb; do
        pdb_filename=$(basename $pdb_file)
        pdb_name=${pdb_filename%.pdb}
        pdb_folder=$SIM_DIR/U*/$pdb_name
        # Create the folder if it doesn't exist
    
        # Move the PDB file into the folder
        cd $pdb_folder
        for i in "${FORCEFIELDS[@]}"; do
                cd $pdb_folder/$i
                sbatch ${JOB_SCRIPT}
        done
done
