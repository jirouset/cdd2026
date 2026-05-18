"""Collect Boltz cofolding results (CIF + JSON) into a single flat directory.

Run on the cluster after predictions finish:
    python collect_boltz_results.py

Skips files that already exist in TARGET_DIR so it is safe to re-run
incrementally (e.g. after a rerun of missing jobs).
"""

import os
import csv
import glob
import shutil
from tqdm import tqdm

# ── paths ────────────────────────────────────────────────────────────────────
BASE_DIR   = os.path.expanduser("~/CXCR4_cofoldings")
TARGET_DIR = os.path.join(BASE_DIR, "results")
SRC_DIR    = os.path.expanduser("~/boltz/all_results")
CSV_PATH   = os.path.join(BASE_DIR, "summary_results.csv")
# ─────────────────────────────────────────────────────────────────────────────

os.makedirs(TARGET_DIR, exist_ok=True)

run_dirs = glob.glob(os.path.join(SRC_DIR, "results_ligand_*_*"))
print(f"Found {len(run_dirs)} result directories.\n")

found_count    = 0
skipped_count  = 0
missing_count  = 0
csv_data       = []

for run_dir in tqdm(run_dirs, desc="Collecting results", unit="folder"):
    folder_name = os.path.basename(run_dir)
    parts = folder_name.split("_")
    if len(parts) < 4:
        continue

    ligand_id = parts[2]
    run_id    = parts[3]

    pred_dir = os.path.join(
        run_dir,
        f"boltz_results_cxcr4_ligand_{ligand_id}",
        "predictions",
        f"cxcr4_ligand_{ligand_id}",
    )
    cif_src  = os.path.join(pred_dir, f"cxcr4_ligand_{ligand_id}_model_0.cif")
    json_src = os.path.join(pred_dir, f"confidence_cxcr4_ligand_{ligand_id}_model_0.json")

    cif_dst  = os.path.join(TARGET_DIR, f"cxcr4_ligand_{ligand_id}_run_{run_id}_model_0.cif")
    json_dst = os.path.join(TARGET_DIR, f"confidence_cxcr4_ligand_{ligand_id}_run_{run_id}_model_0.json")

    if os.path.isfile(cif_src) and os.path.isfile(json_src):
        # Skip if already collected (incremental re-runs)
        if os.path.isfile(cif_dst) and os.path.isfile(json_dst):
            status = "Skipped/Exists"
            skipped_count += 1
        else:
            shutil.copy(cif_src, cif_dst)
            shutil.copy(json_src, json_dst)
            status = "Success"
            found_count += 1
    else:
        status = "Failed/Missing"
        missing_count += 1

    csv_data.append([ligand_id, run_id, status])

with open(CSV_PATH, mode="w", newline="", encoding="utf-8") as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(["Ligand_ID", "Run_ID", "Status"])
    writer.writerows(sorted(csv_data, key=lambda x: (int(x[0]), int(x[1]))))

print("\n" + "-" * 45)
print(f"Newly copied   : {found_count}")
print(f"Already existed: {skipped_count}")
print(f"Missing/failed : {missing_count}")
print(f"CSV saved to   : {CSV_PATH}")
print("-" * 45)
