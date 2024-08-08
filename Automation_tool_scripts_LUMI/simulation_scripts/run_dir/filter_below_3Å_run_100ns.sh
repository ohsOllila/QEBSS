#!/bin/bash


cd ..
SIM_SCRIPTS=${PWD}

filter=$SIM_SCRIPTS/MD_scripts/md_extend_20ns_100ns.sh

cd ..
SIM_DIR=${PWD}


lim=0.3 #limit in ns	

for i in $SIM_DIR/H*
do
		protein_path=$i
		protein=($(basename ${protein_path}))
		cd ${protein_path}
		rmsd_file1=*20ns*rmsd*
		tail $rmsd_file1  >  md_20ns_run_${protein}_rmsd_last_10_points.xvg #Extract 10 last points
		avg1=$( awk 'BEGIN {total = 0}{total+=$2} END {print total/NR}' md_20ns_run_${protein}_rmsd_last_10_points.xvg)	
		fasta=$(grep "CA" topol.top | awk '{print $4}' | tr -d '\n')
                echo "$protein $fasta $avg1" >> ${SIM_DIR}/proteins_20ns_rmsd_avg_fasta.txt
		echo "$protein $avg1" >> ${SIM_DIR}/proteins_20ns_rmsd_avg.txt
		if (( $(echo "$avg1 < $lim" | bc -l) )) && [ ! -f *100ns*gro ]
        	then
			sbatch $filter
			rmsd_file2=*100ns*rmsd*
			tail $rmsd_file2  >  md_100ns_run_${protein}_rmsd_last_10_points.xvg #Extract 10 last points
			avg2=$( awk 'BEGIN {total = 0}{total+=$2} END {print total/NR}' md_100ns_run_${protein}_rmsd_last_10_points.xvg)			
			echo "$protein $fasta $avg2" >> ${SIM_DIR}/proteins_100ns_rmsd_avg_fasta.txt
		        echo "$protein $avg2" >> ${SIM_DIR}/proteins_100ns_rmsd_avg.txt
 		fi

done



