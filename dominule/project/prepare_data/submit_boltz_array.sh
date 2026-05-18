#!/bin/bash
#PBS -N boltz_rerun
#PBS -q all
#PBS -l select=1:ncpus=4:ngpus=1:mem=32gb
#PBS -l walltime=1:00:00
#PBS -m ae
#
# DO NOT set -J here. Submit with the correct range from the shell:
#   N=$(wc -l < missing_runs.txt)
#   qsub -J 1-${N}%30 submit_boltz_array.sh
#
# %30 = 30 jobs running concurrently.
# Each job runs exactly ONE Boltz prediction (one ligand, one run),
# so walltime 1h is sufficient and there are no serial bottlenecks.
# Increase %30 if the cluster allows more concurrent GPU jobs.

MISSING_FILE="/home/senovad/boltz/missing_runs.txt"

# Read the (ligand_id, run_id) pair for this array index
LINE=$(sed -n "${PBS_ARRAY_INDEX}p" "$MISSING_FILE")
LIGAND_ID=$(echo "$LINE" | cut -d',' -f1)
RUN_ID=$(echo "$LINE" | cut -d',' -f2)

if [ -z "$LIGAND_ID" ] || [ -z "$RUN_ID" ]; then
    echo "ERROR: could not parse line ${PBS_ARRAY_INDEX} from ${MISSING_FILE}"
    exit 1
fi

# Environment
source /home/senovad/.bashrc
micromamba activate boltz_env
cd /home/senovad/boltz

INPUT_FILE="/home/senovad/boltz/inputs_cxcr4/cxcr4_ligand_${LIGAND_ID}.yaml"
OUTPUT_DIR="/home/senovad/boltz/all_results/results_ligand_${LIGAND_ID}_${RUN_ID}"

if [ -f "$INPUT_FILE" ]; then
    echo "Ligand ${LIGAND_ID}, run ${RUN_ID} — starting prediction"
    boltz predict "$INPUT_FILE" --out_dir "$OUTPUT_DIR" --override
    echo "Ligand ${LIGAND_ID}, run ${RUN_ID} — done"
else
    echo "ERROR: YAML not found: $INPUT_FILE"
    exit 1
fi
