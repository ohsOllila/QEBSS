#!/bin/bash


cd ..
cd ..
SIM_DIR=${PWD}

SCRIPTS=${SIM_DIR}/simulation_scripts/MD_scripts
JOB_SCRIPT=${SCRIPTS}/mdmat.sh



for i in $SIM_DIR/*
do
  	cd $i
	sbatch ${JOB_SCRIPT}
done
