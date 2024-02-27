#!/bin/bash


sim_time=1500


SCRIPT_DIR=$(cd .. && pwd)
SIM_DIR=$(cd ../.. && pwd)
PROT_DIR=$SIM_DIR/Unst_hydrolase

magn_field=$(awk 'NR==1 {print $6}' "${SCRIPT_DIR}/MD_scripts/hydrolase_exp_data.txt" 2>/dev/null)
make_index=${SCRIPT_DIR}/MD_scripts/makeNHindex.awk
py_script=${SCRIPT_DIR}/PY_scripts/Old_Relaxations_for_Samuli.py
mdmat_plot=${SCRIPT_DIR}/PY_scripts/xpm_plot.py
relax_plot=${SCRIPT_DIR}/PY_scripts/plot_replicas_to_experiment.py

mkdir -p $SIM_DIR/results/$(basename $PROT_DIR)

list=${PROT_DIR}/model*/*/


for i in $list; do

	cd ${i}

	path=${PWD}
	echo $path
	base_dir=$(basename "$(dirname "$path")")
	replicas=$(basename "$base_dir")
	ff=$(basename $path)

	cp $py_script ${path}
	sed -i "s|PATH_TO_CORR|${path}/correlation_functions|" ${path}/Old_Relaxations_for_Samuli.py
	sed -i "s|magn_field=magn_field|magn_field=$magn_field|" ${path}/Old_Relaxations_for_Samuli.py
	TEMP_NAME=(md_${sim_time}ns.tpr)
	name=${TEMP_NAME%.tpr}


	mkdir correlation_functions

	#sed -i.bak 's/ H /HN /g' ${name}.gro
	#awk -f ${make_index} ${name}.gro > HN.ndx

	module purge
	module load gromacs-env

	echo 1 1 | gmx_mpi trjconv -f ${name}.xtc -s ${name}.tpr -pbc mol -center -dump 0 -o temp_${name}.gro
	echo 1 1 | gmx_mpi trjconv -f ${name}.xtc -s ${name}.tpr -pbc mol -center -o ${name}_noPBC.xtc
	gmx_mpi filter -f ${name}_noPBC.xtc -s temp_${name}.gro -nf 20 -all -ol ${name}_smooth.xtc
	echo 1 | gmx_mpi gyrate -s ${name}.tpr -f ${name}_noPBC.xtc -o ${name}_gyrate.xvg

	echo -e "Alpha\nAlpha" | gmx_mpi mdmat -f ${name}.xtc -s ${name}.tpr -mean ${name}_mdmat.xpm
	gmx_mpi xpm2ps -f ${name}_mdmat.xpm -o ${name}_mdmat.eps

	GRO_FILE=(temp_md_1500ns.gro)
	sed -i.bak 's/ H /HN /g' $GRO_FILE
	sed -i.bak 's/H1/HN/g' $GRO_FILE
	awk -f ${make_index} $GRO_FILE > HN.ndx
	numberOFfuncs=$(grep "\[" HN.ndx | tail -n 1 | awk '{print $2}')
	for ((i = 0; i <= $numberOFfuncs; i++))
	do
		if [ ! -s correlation_functions/NHrotaCF_$i.xvg ]; then
			echo $i | gmx_mpi rotacf -f ${name}_noPBC.xtc -s ${name}.tpr -n HN.ndx -o correlation_functions/NHrotaCF_$i.xvg -P 2 -d -xvg none  #-nice 20 &
		fi
	done
	

	module purge
	export PATH="/scratch/project_2003809/cmcajsa/env/bin:$PATH"

	python3 $mdmat_plot
	python3 ${path}/Old_Relaxations_for_Samuli.py > relaxation_data.txt
	module purge
	mkdir -p $SIM_DIR/results/$(basename $PROT_DIR)/$replicas/$ff
	sim_results=$SIM_DIR/results/$(basename $PROT_DIR)/$replicas/$ff
	cp $path/*gyrate*xvg $sim_results
	cp $path/*mdmat* $sim_results
	cp $path/*png $sim_results
	cp $path/*eps  $sim_results
	cp -r $path/correlation_functions $sim_results

done

export PATH="/scratch/project_2003809/cmcajsa/env/bin:$PATH"
python3 $relax_plot

module purge



find $SIM_DIR -name slurm* -exec mv -t $SIM_DIR/output_files {} +
find $SIM_DIR -name array* -exec mv -t $SIM_DIR/output_files {} +
