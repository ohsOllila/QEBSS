#!/bin/bash


cd ..
cd ..
BASE_DIR=$PWD

SIM_DIR=$BASE_DIR/Unst*hyd*


SCRIPTS=$BASE_DIR/simulation_scripts/MD_scripts
md_script=${SCRIPTS}/shorten_traj.sh
cp ${md_script} ${SCRIPTS}/batch_md.sh
JOB_SCRIPT=${SCRIPTS}/batch_md.sh

for i in $SIM_DIR/model*/*; do
	cd $i
	sbatch ${JOB_SCRIPT}
done
