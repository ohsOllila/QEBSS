#!/bin/bash


dir=${PWD}
cd IDPConformerGenerator
python setup.py develop --no-deps 
cd $dir


fasta_files=*.fasta

echo "Select the number that represent the fasta file you want to generate pdb:s for:"

num=1
list=()

for file in $fasta_files; do
	list+=($file)
	echo "("$num")" ${file} 
	((num++))
done

read choice

fasta_name=${list[choice-1]}


culled_lists=("culled_list1" "culled_list2" "culled_list3" "culled_list4" "culled_list5")
replicas=("replica_01.pdb" "replica_02.pdb" "replica_03.pdb" "replica_04.pdb" "replica_05.pdb")


for i in {0..4}; do
    idpconfgen pdbdl culled_lists/${culled_lists[i]} -u -n -d pdbs.tar 
    idpconfgen sscalc pdbs.tar -rd -n 
    idpconfgen torsions sscalc_splitted.tar -sc sscalc.json -o idpconfgen_database.json -n 
    idpconfgen build \
        -db idpconfgen_database.json \
        -seq $fasta_name \
        -nc 10 \
        -et 'pairs' \
        --dstrand \
        --dloop-off \
        -n 
    wait
    
    mv conformer_1.pdb ${replicas[i]}
    rm pdbs.tar
    rm idpconfgen_database.json
    rm idpconfgen.version
    rm energies.log
    rm sscalc*
done

mkdir ../Unst_prot
mv replica*pdb ../Unst_prot


