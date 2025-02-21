#!/bin/bash

cd ..
SIM_SCRIPTS=${PWD}

cd ..
SIM_DIR=${PWD}

script1=${SIM_DIR}/simulation_scripts/MD_scripts/mdmat.sh
script2=${SIM_DIR}/simulation_scripts/PY_scripts/xpm_plot.py

for i in  $SIM_DIR/*
do
  	protein_path=${i}
        protein=($(basename ${i}))
	cd $protein_path
	temp_name=(md*100ns.edr)
        name=${temp_name%.edr}
        python3 ${script1} 
	python3 ${script2}
done



