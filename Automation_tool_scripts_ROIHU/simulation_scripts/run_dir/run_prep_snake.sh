#!/bin/bash -l
#SBATCH --time=12:00:00
#SBATCH --partition=medium
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=128
#SBATCH --output=job_output_%A_%a.log
#SBATCH --error=job_error_%A_%a.log

# -------------------------
# Read arguments
# -------------------------

PROT_NAME="$1"
shift

FORCEFIELD=("$@")

if [[ -z "$PROT_NAME" ]]; then
    echo "Error: No system/protein name supplied."
    exit 1
fi

if (( ${#FORCEFIELD[@]} == 0 )); then
    echo "Error: No force fields supplied."
    exit 1
fi

# -------------------------
# Set paths
# -------------------------

BASE_DIR="$PWD"
PROT_DIR="$BASE_DIR/$PROT_NAME"

SCRIPTS="${BASE_DIR}/simulation_scripts/MD_scripts"
SNAKEFILE="${SCRIPTS}/Snakefile"

if [[ ! -f "$SNAKEFILE" ]]; then
    echo "Error: Snakefile not found: $SNAKEFILE"
    exit 1
fi

if [[ ! -d "$PROT_DIR" ]]; then
    echo "Error: System directory not found: $PROT_DIR"
    exit 1
fi

# -------------------------
# Load modules
# -------------------------

module purge --force

module load snakemake
module load gromacs-env


# -------------------------
# Find PDB files
# -------------------------

shopt -s nullglob
PDB_FILES=("$PROT_DIR"/rep*.pdb)
shopt -u nullglob

NUM_PDBS=${#PDB_FILES[@]}
NUM_FORCEFIELDS=${#FORCEFIELD[@]}
EXPECTED_JOBS=$((NUM_PDBS * NUM_FORCEFIELDS))

if (( NUM_PDBS == 0 )); then
    echo "Error: No rep*.pdb files found in $PROT_DIR"
    exit 1
fi

# -------------------------
# Print configuration
# -------------------------

echo "======================================"
echo "System:          $PROT_NAME"
echo "PDB files:       $NUM_PDBS"
echo "Force fields:    ${FORCEFIELD[*]}"
echo "Number of fields: $NUM_FORCEFIELDS"
echo "Expected jobs:   $EXPECTED_JOBS"
echo "Array task ID:   ${SLURM_ARRAY_TASK_ID:-not set}"
echo "======================================"

# -------------------------
# Prepare directories
# -------------------------

PROT_FOLDERS=()

for pdb_file in "${PDB_FILES[@]}"; do

    replicas="$(basename "${pdb_file%.pdb}")"
    directory_path="${pdb_file%/*}"

    for ff in "${FORCEFIELD[@]}"; do

        ff_dir="$directory_path/$replicas/$ff"

        mkdir -p "$ff_dir"

        cp -u -p "$pdb_file" "$ff_dir/"

        PROT_FOLDERS+=("$ff_dir")

    done
done

# -------------------------
# Validate array task ID
# -------------------------

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "Error: SLURM_ARRAY_TASK_ID is not set."
    echo "Submit this script through the helper script using sbatch."
    exit 1
fi

TASK_ID="$SLURM_ARRAY_TASK_ID"

if ! [[ "$TASK_ID" =~ ^[0-9]+$ ]]; then
    echo "Error: Invalid array task ID: $TASK_ID"
    exit 1
fi

if (( TASK_ID >= ${#PROT_FOLDERS[@]} )); then
    echo "Error: Array task $TASK_ID exceeds the number of available jobs (${#PROT_FOLDERS[@]})."
    exit 1
fi

# -------------------------
# Select working directory
# -------------------------

WORK_DIR="${PROT_FOLDERS[$TASK_ID]}"

cd "$WORK_DIR" || {
    echo "Error: Failed to enter directory: $WORK_DIR"
    exit 1
}

echo
echo "======================================"
echo "Job ID:          ${SLURM_JOB_ID:-unknown}"
echo "Array task:      $TASK_ID"
echo "Working dir:     $(pwd)"
echo "======================================"
echo

# -------------------------
# Run Snakemake
# -------------------------

snakemake \
    -s "$SNAKEFILE" \
    --cores "${SLURM_CPUS_PER_TASK:-128}"
