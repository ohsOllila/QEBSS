#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --partition=medium
#SBATCH --ntasks=25
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --account=Project_2001058

sim_time=1000

cd ..
cd ..
SIM_DIR=${PWD}

#FORCEFIELDS=(AMBER03WS AMBER99SB-DISP AMBER99SBWS CHARMM36 DESAMBER)
FORCEFIELDS=(AMBER03WS)

list=()

for pdb_file in *pdb; do
        pdb_filename=$(basename $pdb_file)
        pdb_name=${pdb_filename%.pdb}
        pdb_folder=$SIM_DIR/$pdb_name
        # Create the folder if it doesn't exist
    
        # Move the PDB file into the folder
        cd $pdb_folder
        for i in "${FORCEFIELDS[@]}"; do
                cd $pdb_folder/$i
       		list+=(${PWD})         
        done
done

srun gmx_mpi mdrun -multidir ${list[@]} -deffnm md_run_${sim_time}ns -cpi md_run_${sim_time}ns.cpt -dlb yes
