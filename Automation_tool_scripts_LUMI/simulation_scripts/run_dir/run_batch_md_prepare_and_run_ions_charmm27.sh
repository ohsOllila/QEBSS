#!/bin/bash


echo "Input simulation time"
read time_input


cd ..
cd ..
SIM_DIR=${PWD}

SCRIPTS=${SIM_DIR}/simulation_scripts/MD_scripts
md_script=${SCRIPTS}/md_prepare_and_run_ions_charmm27.sh
cp ${md_script} ${SCRIPTS}/batch_md.sh


JOB_SCRIPT=${SCRIPTS}/batch_md.sh

sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"


for i in $SIM_DIR/H*
do
  	cd $i
	sbatch ${JOB_SCRIPT}
done
