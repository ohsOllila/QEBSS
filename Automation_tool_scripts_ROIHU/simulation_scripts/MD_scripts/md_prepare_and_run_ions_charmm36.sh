#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --partition=small
#SBATCH --ntasks-per-node=64
#SBATCH --mem-per-cpu=400
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --array=0-2
#SBATCH --output=array_job_output_%A_%a.txt
#SBATCH --account=Project_462000199


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GMXLIB=/scratch/project_462000154/cmcajsa/CHARMM36
export GMX_MAXBACKUP=-1

module use /appl/local/csc/modulefiles
module load gromacs/2023.1-hipsycl

FORCEFIELD=CHARMM36

# set variables
SUFFIX=c36
water_model=tip3p
forcefield=charmm36

# set base directories
PARAM_DIR=/scratch/project_462000154/cmcajsa/CHARMM36
sim_time=sim_time

ions=(FE2P CU2P MN2P)
ions2=(FE2P CU2P MN2P)

for i in "${ions[@]}"; do
	mkdir ${i}
done

#cp *pdb ${i}
cd ${ions[${SLURM_ARRAY_TASK_ID}]}

i=$(basename $PWD)	

temp_name=(*.pdb)
PROTEIN=${temp_name%.pdb}

srun -n 1 gmx_mpi -quiet pdb2gmx -f ${PROTEIN}.pdb -o ${PROTEIN}.gro -water ${water_model} -ff ${forcefield} -ignh 
srun -n 1 gmx_mpi -quiet editconf -f ${PROTEIN}.gro -o ${PROTEIN}_newbox.gro -c -d 1.0 -bt cubic
srun -n 1 gmx_mpi -quiet solvate -cp ${PROTEIN}_newbox.gro -cs spc216.gro -o ${PROTEIN}_solv.gro -p topol.top

srun -n 1 gmx_mpi -quiet grompp -f ${PARAM_DIR}/ions_${SUFFIX}.mdp -c ${PROTEIN}_solv.gro -p topol.top -o ions_${SUFFIX}_${i}.tpr

if [[ "${ions1[*]}" =~ "${i}" ]]; then
	echo SOL | srun -n 1 gmx_mpi -quiet genion -s ions_${SUFFIX}_${i}.tpr -o ${PROTEIN}_solv_ions_${i}.gro -p topol.top  -pname ${i} -nname CLA -neutral -conc 1  
else
	echo SOL | srun -n 1 gmx_mpi -quiet genion -s ions_${SUFFIX}_${i}.tpr -o ${PROTEIN}_solv_ions_${i}.gro -p topol.top  -pname ${i} -pq 2 -nname CLA -neutral -conc 1 
fi

srun -n 1 gmx_mpi -quiet grompp -f ${PARAM_DIR}/em_${SUFFIX}.mdp -c ${PROTEIN}_solv_ions_${i}.gro -p topol.top -o em_${SUFFIX}_${i}.tpr 
srun gmx_mpi -quiet mdrun -deffnm em_${SUFFIX}_${i} 

srun -n 1 gmx_mpi -quiet grompp -f ${PARAM_DIR}/nvt_${SUFFIX}.mdp -c em_${SUFFIX}_${i}.gro -r em_${SUFFIX}_${i}.gro -p topol.top -o nvt_${SUFFIX}_${i}.tpr 
srun gmx_mpi -quiet mdrun -deffnm nvt_${SUFFIX}_${i} 

srun -n 1 gmx_mpi -quiet grompp -f ${PARAM_DIR}/npt_${SUFFIX}.mdp -c nvt_${SUFFIX}_${i}.gro -r nvt_${SUFFIX}_${i}.gro -t nvt_${SUFFIX}_${i}.cpt -p topol.top -o npt_${SUFFIX}_${i}.tpr 
srun gmx_mpi -quiet mdrun -deffnm npt_${SUFFIX}_${i}


if [ -f md*${sim_time}ns.xtc ] && [ ! -f md*${sim_time}ns.gro ]; then 
	srun gmx_mpi -quiet mdrun -deffnm md_run_${PROTEIN}_${SUFFIX}_${i}_${sim_time}ns -dlb yes -cpi md_run_${PROTEIN}_${SUFFIX}_${i}_${sim_time}ns.cpt
elif [ ! -f md_*${sim_time}ns.xtc ] && [ ! -f md_*${sim_time}ns.gro ]; then 
	srun gmx_mpi -quiet grompp -f ${PARAM_DIR}/md_diff_sim_time/md_${SUFFIX}_${sim_time}ns.mdp -c npt_${SUFFIX}.gro -t npt_${SUFFIX}.cpt -p topol.top -o md_run_${PROTEIN}_${SUFFIX}_${i}_${sim_time}ns.tpr
 	srun gmx_mpi -quiet mdrun -deffnm md_run_${PROTEIN}_${SUFFIX}_${i}_${sim_time}ns -dlb yes
fi

temp_name=(md_*${sim_time}ns.edr)
name=${temp_name%.edr}

echo 1 1 | gmx_mpi -quiet trjconv -s ${name}.tpr -f ${name}.xtc -o ${name}_noPBC.xtc -pbc mol -center #Correct for periodic boundary condition (PBC)
echo 4 4 | gmx_mpi -quiet rms -s ${name}.tpr -f ${name}_noPBC.xtc -o ${name}_rmsd.xvg -tu ns #compute RMSD

echo -e "Protein\nProtein" | gmx_mpi trjconv -f ${name}.xtc -s ${name}.tpr -pbc mol -center -dump 0 -o temp_${name}.gro
gmx_mpi filter -f ${name}_noPBC.xtc -s temp_${name}.gro -nf 20 -all -ol md_smooth_${name}.xtc

seff $SLURM_JOBID
