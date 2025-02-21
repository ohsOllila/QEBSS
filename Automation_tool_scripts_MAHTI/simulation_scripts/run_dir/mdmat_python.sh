#!/bin/bash

cd ..
SIM_SCRIPTS=${PWD}

cd ..
SIM_DIR=${PWD}

script=${SIM_DIR}/simulation_scripts/PY_scripts/get_mdmat.py


for i in  `seq -f "%03g" 177`
do
  	protein_path=${SIM_DIR}/*${i}*
        protein=($(basename ${protein_path}))
	cd $protein_path
	temp_name=(md*100ns.edr)
        name=${temp_name%.edr}

        python3 ${script} ${name}.xtc ${protein}.gro

done



