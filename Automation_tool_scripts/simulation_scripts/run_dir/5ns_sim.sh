#!/bin/bash

cd ..
SIM_SCRIPTS=${PWD}

script=$SIM_SCRIPTS/MD_scripts/sim_last_5ns.sh

cd ..
SIM_DIR=${PWD}

mkdir results/5ns
mkdir results/5ns/xtc
#mkdir results/5ns/gro
mkdir results/5ns/pdb

#xtc_dir=$SIM_DIR/results/5ns/xtc
#gro_dir=$SIM_DIR/results/5ns/gro
pdb_dir=$SIM_DIR/results/5ns/pdb

for i in $SIM_DIR/U*
do
    cd $i
	if [ ! -f *5ns_tail.pdb ]; then
		sbatch $script
	fi
done

#find $SIM_DIR -name '*5ns_tail.xtc' -exec cp {} $xtc_dir \;
#find $SIM_DIR -name '*5ns_tail.gro' -exec cp {} $gro_dir \;
find $SIM_DIR -name '*5ns_tail.pdb' -exec cp {} $pdb_dir \;
