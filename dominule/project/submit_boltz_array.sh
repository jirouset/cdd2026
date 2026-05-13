#!/bin/bash
#PBS -N boltz_cofolds
#PBS -q all
#PBS -l select=1:ncpus=4:ngpus=1:mem=32gb
#PBS -l walltime=5:00:00
#PBS -m ae
#PBS -J 1-2000%10

## Note: The %50 above means "Run 50 at a time."
## You can increase this to %50 if you want it even faster.

# 1. Environment Setup
source /home/senovad/.bashrc
micromamba activate boltz_env

# 2. Directory Setup
cd /home/senovad/boltz

# 3. Handle the Index
# This matches your file naming: cxcr4_ligand_1.yaml, etc.
INPUT_FILE="/home/senovad/boltz/inputs_cxcr4/cxcr4_ligand_${PBS_ARRAY_INDEX}.yaml"
OUTPUT_DIR="/home/senovad/boltz/all_results/results_ligand_${PBS_ARRAY_INDEX}"

# 4. Execution
if [ -f "$INPUT_FILE" ]; then
  for (( i = 0; i < 5; i++ )); do
    echo "Starting prediction for ligand ${PBS_ARRAY_INDEX} iteration $i"
    boltz predict "$INPUT_FILE" --out_dir "$OUTPUT_DIR"_"$i" --override
  done

else
    echo "File $INPUT_FILE not found, skipping."
    exit 1
fi