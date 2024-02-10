#!/bin/bash

echo "Choose simulation time (ns) to calculate rmsd, (10, 20, 100)"
read time_input

cd ..
SIM_SCRIPTS=${PWD}

cd ..
SIM_DIR=${PWD}

rog_script=${SIM_SCRIPTS}/MD_scripts/rog.sh
cp ${rog_script} ${SIM_SCRIPTS}/batch_rog.sh


JOB_SCRIPT=${SIM_SCRIPTS}/batch_rog.sh


sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"

cd ${SIM_SCRIPTS}

sbatch ${JOB_SCRIPT}

