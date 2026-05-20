"""Build the CXCR4 protein-ligand contact dataset from Boltz cofolding CIF files.

Run from dominule/project/:
    python build_contacts.py

Output: datasets/CXCR4_contacts.csv
"""

import time
import warnings
import pandas as pd

# Biotite uses label_atom_id as a fallback when auth_atom_id is absent in Boltz CIF files.
# This is expected and harmless — suppress to keep output readable.
warnings.filterwarnings("ignore", message=".*auth_atom_id.*", category=UserWarning)

from contact_detection import create_ligand_contact_dataset

BOLTZ_CSV   = "datasets/CXCR4_chembl_ligands_boltz.csv"
LABELS_CSV  = "datasets/CXCR4_for_ml.csv"
BOLTZ_DIR   = "boltz_cofoldings/"
OUTPUT_CSV  = "datasets/CXCR4_contacts.csv"

df_boltz = pd.read_csv(BOLTZ_CSV, sep=";")
df_labels = pd.read_csv(LABELS_CSV)

# Only process ligands that have a label — unlabelled ones are unused in analysis
labelled_smiles = set(df_labels["Original Smiles"])
df_boltz_filtered = df_boltz[df_boltz["Smiles"].isin(labelled_smiles)].copy()

n_all  = df_boltz["Ligand_ID"].nunique()
n_kept = df_boltz_filtered["Ligand_ID"].nunique()
print(f"Ligands in boltz CSV : {n_all}")
print(f"Ligands with labels  : {n_kept}  (skipping {n_all - n_kept})")

t0 = time.time()
df_contacts = create_ligand_contact_dataset(
    df_boltz_filtered,
    boltz_dir=BOLTZ_DIR,
    df_labels=df_labels,
)
elapsed = time.time() - t0

df_contacts.to_csv(OUTPUT_CSV, index=False)
print(f"\nSaved to {OUTPUT_CSV}  —  shape {df_contacts.shape}")
print(f"Elapsed: {elapsed/60:.1f} min")
