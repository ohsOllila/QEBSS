#!/bin/bash

echo "Choose simulation time (ns) to calculate rmsd, (10, 20, 100)"
read time_input

cd ..
cd ..
SIM_DIR=${PWD}

SCRIPTS=${SIM_DIR}/simulation_scripts/MD_scripts
sasa_script=${SCRIPTS}/sasa.sh
cp ${sasa_script} ${SCRIPTS}/batch_sasa.sh


JOB_SCRIPT=${SCRIPTS}/batch_sasa.sh

sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"


for i in $SIM_DIR/*
do
  	cd $i
	sbatch ${JOB_SCRIPT}
done
