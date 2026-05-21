#!/bin/bash
#PBS -N boltz_cxcr4
#PBS -l select=1:ncpus=4:ngpus=1:mem=32gb
#PBS -l walltime=04:00:00
#PBS -m ae

source /home/senovad/.bashrc
micromamba activate boltz_env
cd /home/senovad/boltz

echo "Starting Boltz batch run from /home/senovad/boltz/inputs"

for input_file in /home/senovad/boltz/inputs/cxcr4_ligand_*.yaml; do

    base_name=$(basename "$input_file" .yaml)

    echo "Current file: $input_file"
    echo "Running prediction for: $base_name"

    # We output directly to a results folder in boltz/
    boltz predict "$input_file" --out_dir "/home/senovad/boltz/all_results/results_$base_name"

done

echo "Batch job finished."