#!/bin/bash

echo "Input simulation time"
read time_input

cd ..
cd ..

SIM_DIR=${PWD}

SCRIPTS=${SIM_DIR}/simulation_scripts/MD_scripts
rhosc_script=${SCRIPTS}/rhosc.sh
cp ${rhosc_script} ${SCRIPTS}/batch_rhosc.sh


JOB_SCRIPT=${SCRIPTS}/batch_rhosc.sh

sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"


for i in $SIM_DIR/*
do
  	cd $i
	sbatch ${JOB_SCRIPT}
done
