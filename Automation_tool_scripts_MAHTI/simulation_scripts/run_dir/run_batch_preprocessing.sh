#!/bin/bash

cd ..
cd ..
SIM_DIR=${PWD}
batch=${SIM_DIR}/simulation_scripts/MD_scripts/preprocessing.sh

for i in $SIM_DIR/*
do
        cd $i
        sbatch $batch
done


