"""Contact detection for protein-small molecule Boltz cofolding structures (CXCR4)."""

import os
import warnings
from typing import List, Optional

import numpy as np
import pandas as pd
import biotite.structure.io as strucio
from tqdm import tqdm


def inspect_chains(cif_path: str) -> dict:
    """Return chain IDs, residue types, and atom counts found in a CIF file.

    Use this to verify that protein_chain and ligand_chain defaults are correct
    before running a full analysis.
    """
    structure = strucio.load_structure(cif_path)
    chain_info = {}
    for chain_id in np.unique(structure.chain_id):
        mask = structure.chain_id == chain_id
        res_names = np.unique(structure.res_name[mask]).tolist()
        chain_info[chain_id] = {
            "residue_types": res_names,
            "n_atoms": int(mask.sum()),
        }
    return chain_info


def detect_protein_ligand_contacts(
    cif_path: str,
    protein_chain: str = "A",
    ligand_chain: str = "B",
    distance_threshold: float = 4.5,
) -> List[int]:
    """Detect protein residue IDs in contact with a small-molecule ligand.

    Computes heavy-atom pairwise distances between protein_chain and ligand_chain.

    Returns sorted list of unique protein residue IDs where at least one heavy atom
    is within distance_threshold Angstroms of any ligand heavy atom.
    """
    structure = strucio.load_structure(cif_path)

    # Exclude hydrogens — their positions are model-dependent and often absent
    is_heavy = ~np.char.startswith(structure.atom_name.astype(str), "H")
    structure = structure[is_heavy]

    protein = structure[structure.chain_id == protein_chain]
    ligand = structure[structure.chain_id == ligand_chain]

    if len(protein) == 0 or len(ligand) == 0:
        raise ValueError(
            f"Chain '{protein_chain}' ({len(protein)} atoms) or "
            f"chain '{ligand_chain}' ({len(ligand)} atoms) is empty in {cif_path}"
        )

    # Distance matrix shape: (n_protein_atoms, n_ligand_atoms)
    diff = protein.coord[:, np.newaxis, :] - ligand.coord[np.newaxis, :, :]
    distmat = np.linalg.norm(diff, axis=2)

    # Protein atom is "in contact" if it is within threshold of any ligand atom
    in_contact = np.any(distmat <= distance_threshold, axis=1)
    contacted_residues = sorted(set(int(r) for r in protein.res_id[in_contact]))
    return contacted_residues


def analyze_ligand_run(
    cif_path: str,
    protein_chain: str = "A",
    ligand_chain: str = "B",
    distance_threshold: float = 4.5,
) -> Optional[dict]:
    """Analyze one Boltz run CIF file and return contact statistics.

    Returns None if the file is missing or parsing fails (graceful degradation for
    failed Boltz runs). Callers are responsible for tracking and reporting missing files.
    """
    if not os.path.isfile(cif_path):
        return None

    try:
        contacted = detect_protein_ligand_contacts(
            cif_path=cif_path,
            protein_chain=protein_chain,
            ligand_chain=ligand_chain,
            distance_threshold=distance_threshold,
        )
        return {
            "contacted_residues": contacted,
            "contact_count": len(contacted),
            "cif_path": cif_path,
        }
    except Exception as e:
        warnings.warn(f"Error processing {cif_path}: {e}")
        return None


def create_ligand_contact_dataset(
    df_boltz: pd.DataFrame,
    boltz_dir: str,
    df_labels: Optional[pd.DataFrame] = None,
    protein_chain: str = "A",
    ligand_chain: str = "B",
    distance_threshold: float = 4.5,
) -> pd.DataFrame:
    """Aggregate protein-ligand contacts across all Boltz runs for each ligand.

    Parameters
    ----------
    df_boltz : DataFrame loaded from CXCR4_chembl_ligands_boltz.csv (sep=';').
               Must contain columns: Ligand_ID, Run_ID, boltz_cif_file, Smiles.
    boltz_dir : Directory containing the .cif files referenced in boltz_cif_file.
    df_labels : DataFrame with 'Original Smiles' and 'Inhibitor' columns (e.g.
                loaded from CXCR4_for_ml.csv). If None, 'Inhibitor' must already
                be present in df_boltz.
    protein_chain : Chain ID for the CXCR4 protein in the CIF (default 'A').
    ligand_chain  : Chain ID for the small-molecule ligand in the CIF (default 'B').
    distance_threshold : Angstrom cutoff for heavy-atom contact detection (default 4.5).

    Returns
    -------
    DataFrame with columns:
        Ligand_ID | Inhibitor | contact_count_mean | runs_processed | <res_id> ...

    Each residue column contains the number of runs (0 to 5) in which that CXCR4
    residue was in contact with the ligand. Residues never contacted by any ligand
    are excluded from the output.
    """
    # --- Attach Inhibitor labels via SMILES join if df_labels provided ---
    if df_labels is not None:
        label_map = dict(zip(df_labels["Original Smiles"], df_labels["Inhibitor"]))
        df_boltz = df_boltz.copy()
        df_boltz["Inhibitor"] = df_boltz["Smiles"].map(label_map)
        n_matched = df_boltz["Inhibitor"].notna().sum()
        n_total = len(df_boltz)
        print(f"Label join: {n_matched}/{n_total} rows matched ({n_total - n_matched} unmatched).")

    if "Inhibitor" not in df_boltz.columns:
        raise ValueError(
            "No 'Inhibitor' column found. Provide df_labels with columns "
            "'Original Smiles' and 'Inhibitor', or ensure df_boltz already contains 'Inhibitor'."
        )

    grouped = df_boltz.groupby("Ligand_ID")
    all_residues: set = set()
    ligand_results: dict = {}

    # --- First pass: run contact detection, collect all residue IDs seen ---
    missing_files: List[str] = []
    for ligand_id, group in tqdm(grouped, desc="Detecting contacts"):
        inhibitor = group["Inhibitor"].iloc[0]
        run_results = []
        for _, row in group.iterrows():
            cif_path = os.path.join(boltz_dir, row["boltz_cif_file"])
            result = analyze_ligand_run(
                cif_path=cif_path,
                protein_chain=protein_chain,
                ligand_chain=ligand_chain,
                distance_threshold=distance_threshold,
            )
            if result is not None:
                run_results.append(result)
                all_residues.update(result["contacted_residues"])
            else:
                missing_files.append(cif_path)

        ligand_results[ligand_id] = {"inhibitor": inhibitor, "runs": run_results}

    if missing_files:
        n_ligands_missing = sum(
            1 for data in ligand_results.values() if len(data["runs"]) == 0
        )
        print(
            f"\nMissing CIF files: {len(missing_files)} runs across "
            f"{n_ligands_missing} ligands with no successful run at all."
        )

    all_residues_sorted = sorted(all_residues)

    # --- Second pass: build feature rows with consistent residue columns ---
    final_rows = []
    for ligand_id, data in ligand_results.items():
        runs = data["runs"]

        residue_counts = {res: 0 for res in all_residues_sorted}
        for run in runs:
            for res in run["contacted_residues"]:
                residue_counts[res] += 1

        mean_contacts = float(np.mean([r["contact_count"] for r in runs])) if runs else 0.0

        row = {
            "Ligand_ID": ligand_id,
            "Inhibitor": data["inhibitor"],
            "contact_count_mean": mean_contacts,
            "runs_processed": len(runs),
        }
        row.update(residue_counts)
        final_rows.append(row)

    output_df = pd.DataFrame(final_rows)
    meta_cols = ["Ligand_ID", "Inhibitor", "contact_count_mean", "runs_processed"]
    return output_df[meta_cols + all_residues_sorted]
