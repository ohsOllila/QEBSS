#!/bin/bash

cd ..
cd ..
BASE_DIR=$PWD

SIM_DIR=$BASE_DIR/Unst_prot


echo "Input simulation time in ns:"
read time_input


SCRIPTS=$BASE_DIR/simulation_scripts/MD_scripts
md_script=${SCRIPTS}/md_standard.sh
cp ${md_script} ${SCRIPTS}/batch_md.sh
JOB_SCRIPT=${SCRIPTS}/batch_md.sh

sed -i "s/sim_time=sim_time/sim_time=${time_input}/" "${JOB_SCRIPT}"
sed -i "s/nr_nodes/$NODE_COUNT/" "${JOB_SCRIPT}"

for i in $SIM_DIR/*/*/; do
    cd "$i" || continue
    id="${i}"

    # Check job status
    job_status=$(squeue -u "$USER" -n "$id" -h -o "%T")

    # If job is not running and no completed output file exists
    #if [[ "$job_status" != "R" && ! -e md_*ns.gro ]]; then
    if [[ "$job_status" != "R" ]]; then
        # If the job failed before, decrease node count
        if [[ "$job_status" == "FAILED" || "$job_status" == "CANCELLED" ]]; then
            ((NODE_COUNT--))  # Decrease node count by 1
            sed -i "s/10/${NODE_COUNT}/" "${JOB_SCRIPT}"
            echo "Job $id failed. Decreasing node count to ${NODE_COUNT} and restarting."
        fi

        # Submit job
        sbatch --job-name="$id" "${JOB_SCRIPT}"
    fi
done
