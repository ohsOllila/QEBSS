#!/bin/bash


echo "Input simulation time"
read time_input

cd ..
SIM_SCRIPTS=${PWD}

cd ..
SIM_DIR=${PWD}


SCRIPTS=${SIM_DIR}/simulation_scripts/MD_scripts
md_script=${SCRIPTS}/md_multidir_forcefields.sh
cp ${md_script} ${SCRIPTS}/batch_md.sh


JOB_SCRIPT=${SCRIPTS}/batch_md.sh

sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"


sbatch ${JOB_SCRIPT}





