#!/bin/bash

cd ..
SIM_SCRIPTS=${PWD}

cd ..
SIM_DIR=${PWD}


mkdir results
cd results



for i in  $SIM_DIR/*
do
        protein_path=$i
        protein=($(basename $i))
        mkdir ${protein}
	cp ${protein_path}/*pdb ${SIM_DIR}/results/${protein}/
	cp ${protein_path}/*gro ${SIM_DIR}/results/${protein}/
	cp ${protein_path}/*clean* ${SIM_DIR}/results/${protein}/
	cp ${protein_path}/md*ns.tpr ${SIM_DIR}/results/${protein}/
        cp ${protein_path}/md*ns.xtc ${SIM_DIR}/results/${protein}/
        cp ${protein_path}/*smooth* ${SIM_DIR}/results/${protein}/
	cp ${protein_path}/*mdmat.xpm ${SIM_DIR}/results/${protein}/
	cp ${protein_path}/*mdmat.eps ${SIM_DIR}/results/${protein}/
	cp ${protein_path}/*rhosc* ${SIM_DIR}/results/${protein}/
	cp ${protein_path}/*rmsd* ${SIM_DIR}/results/${protein}/
done

