#!/bin/bash



cd ..
cd ..
SIM_DIR=${PWD}

SCRIPTS=${SIM_DIR}/simulation_scripts/MD_scripts
analysis_script=${SCRIPTS}/make_relaxation_and_plot.sh
cp ${analysis_script} ${SCRIPTS}/batch_analysis.sh


JOB_SCRIPT=${SCRIPTS}/batch_analysis.sh



for i in $SIM_DIR/U*
do
  	cd $i
	sbatch ${JOB_SCRIPT}
done

