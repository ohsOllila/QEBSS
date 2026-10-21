#!/bin/bash

# -------------------------
# Set paths
# -------------------------

cd ../.. || exit 1

BASE_DIR="$PWD"
SCRIPTS="$BASE_DIR/simulation_scripts/run_dir"
SNAKEFILE="$SCRIPTS/run_prep_snake.sh"

# -------------------------
# Define force fields
# -------------------------

FORCEFIELD=(
    AMBER03WS
    AMBER99SB-DISP
    AMBER99SBWS
    CHARMM36M
    DESAMBER
)

# -------------------------
# Select system
# -------------------------

echo "Available systems:"
echo

SYSTEMS=()
num=1

while IFS= read -r dir; do
    [ -d "$dir" ] || continue

    system=$(basename "$dir")
    SYSTEMS+=("$system")

    echo "$num) $system"
    ((num++))
done < <(
    find "$BASE_DIR" \
        -mindepth 1 \
        -maxdepth 1 \
        -type d \
        -print | sort -V
)

if (( ${#SYSTEMS[@]} == 0 )); then
    echo "Error: No systems found in $BASE_DIR"
    exit 1
fi

echo
read -rp "Select system: " choice

if ! [[ "$choice" =~ ^[0-9]+$ ]] ||
   (( choice < 1 || choice > ${#SYSTEMS[@]} )); then
    echo "Invalid system selection."
    exit 1
fi

SYSTEM="${SYSTEMS[$((choice-1))]}"
PROT_DIR="$BASE_DIR/$SYSTEM"

# -------------------------
# Select project
# -------------------------

echo
echo "Available projects:"
echo

PROJECTS=()
num=1

your_projects=$(lumi-allocations | grep -o "project_.*" | awk '{print $1}')

while IFS= read -r project; do
    [ -n "$project" ] || continue

    PROJECTS+=("$project")

    echo "$num) $project"
    ((num++))
done <<< "$your_projects"

if (( ${#PROJECTS[@]} == 0 )); then
    echo "Error: No projects found."
    exit 1
fi

echo
read -rp "Select project: " choice

if ! [[ "$choice" =~ ^[0-9]+$ ]] ||
   (( choice < 1 || choice > ${#PROJECTS[@]} )); then
    echo "Invalid project selection."
    exit 1
fi

PROJECT="${PROJECTS[$((choice-1))]}"

# -------------------------
# Check job script
# -------------------------

if [[ ! -f "$SNAKEFILE" ]]; then
    echo "Error: Job script not found:"
    echo "$SNAKEFILE"
    exit 1
fi

# -------------------------
# Count PDB files
# -------------------------

shopt -s nullglob
PDB_FILES=("$PROT_DIR"/rep*.pdb)
shopt -u nullglob

NUM_PDBS=${#PDB_FILES[@]}
NUM_FORCEFIELDS=${#FORCEFIELD[@]}

if (( NUM_PDBS == 0 )); then
    echo "Error: No rep*.pdb files found in:"
    echo "$PROT_DIR"
    exit 1
fi

if (( NUM_FORCEFIELDS == 0 )); then
    echo "Error: No force fields defined."
    exit 1
fi

# -------------------------
# Calculate number of jobs
# -------------------------

NJOBS=$((NUM_PDBS * NUM_FORCEFIELDS))

ARRAY_RANGE="0-$((NJOBS - 1))"

# -------------------------
# Display configuration
# -------------------------

echo
echo "========================================"
echo "Submission configuration"
echo "========================================"
echo "System:           $SYSTEM"
echo "System directory: $PROT_DIR"
echo "Project:          $PROJECT"
echo
echo "PDB files found:  $NUM_PDBS"
echo "Force fields:     ${FORCEFIELD[*]}"
echo "Number of fields: $NUM_FORCEFIELDS"
echo "Total jobs:       $NJOBS"
echo "Array range:      $ARRAY_RANGE"
echo "========================================"
echo

echo "PDB files:"
for pdb_file in "${PDB_FILES[@]}"; do
    echo "  $(basename "$pdb_file")"
done

echo
read -rp "Submit these jobs? [y/N]: " confirm

if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    echo "Submission cancelled."
    exit 0
fi

# -------------------------
# Submit Slurm job array
# -------------------------

echo
echo "Submitting job array..."

sbatch \
    --account="$PROJECT" \
    --array="$ARRAY_RANGE" \
    "$SNAKEFILE" \
    "$SYSTEM" \
    "${FORCEFIELD[@]}"

STATUS=$?

if (( STATUS == 0 )); then
    echo
    echo "Job array submitted successfully."
else
    echo
    echo "Error: Job submission failed."
    exit "$STATUS"
fi
