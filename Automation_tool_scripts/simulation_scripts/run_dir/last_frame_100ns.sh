#!/bin/bash

cd ..
SIM_SCRIPTS=${PWD}

script=$SIM_SCRIPTS/MD_scripts/sim_lastframe_100ns.sh

cd ..
SIM_DIR=${PWD}

mkdir results/last_frame_100ns
#mkdir results/5ns/xtc
#mkdir results/5ns/gro
mkdir results/last_frame_100ns/pdb

#xtc_dir=$SIM_DIR/results/5ns/xtc
#gro_dir=$SIM_DIR/results/5ns/gro
pdb_dir=$SIM_DIR/results/last_frame_100ns/pdb

for i in $SIM_DIR/*
do
    cd $i
	if [ ! -f Frame_100ns* ]; then
		sbatch $script
	fi
done

#find $SIM_DIR -name '*5ns_tail.xtc' -exec cp {} $xtc_dir \;
#find $SIM_DIR -name '*5ns_tail.gro' -exec cp {} $gro_dir \;
find $SIM_DIR -name 'Frame_100ns*' -exec cp {} $pdb_dir \;
