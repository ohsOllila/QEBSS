#!/bin/bash


cd ..
cd ..
SIM_DIR=${PWD}


for i in $SIM_DIR/*
do
	cd $i
	pymol ~/distancetohelix.py
	cd ..
done
