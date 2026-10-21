#!/bin/bash


echo "Input simulation time"
read time_input


cd ..
cd ..
SIM_DIR=${PWD}

SCRIPTS=${SIM_DIR}/simulation_scripts/MD_scripts
md_script=${SCRIPTS}/md_prep_and_sim.sh
cp ${md_script} ${SCRIPTS}/batch_md.sh


JOB_SCRIPT=${SCRIPTS}/batch_md.sh

sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"


for i in $SIM_DIR/Helices_Phase2_set2_0007_unrelaxed_rank_1_model_5
do
  	cd $i
	sbatch ${JOB_SCRIPT}
done

