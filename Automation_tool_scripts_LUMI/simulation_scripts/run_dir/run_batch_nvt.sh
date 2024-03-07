#!/bin/bash
cd ..
cd ..
SIM_DIR=${PWD}
batch=${SIM_DIR}/simulation_scripts/MD_scripts/nvt.sh

for i in $SIM_DIR/*
do
  	cd $i
	if [ ! -f nvt*gro ]; then
		sbatch $batch
	fi
done
