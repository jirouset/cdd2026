from ast import List
import pandas as pd
from tqdm import tqdm
from antibody_contact_detection import count_target_residue_contacts


def create_epitope_occurrence_dataset(
        df_refolds: pd.DataFrame,
        possible_epitope_ids: List[int],
        antigen_chain: str = 'B',
        antibody_chain: str = 'A'
):
    # 1. Initialize
    clean_epitope_ids = [int(res_id) for res_id in possible_epitope_ids]
    target_set = set(clean_epitope_ids)

    valid_types = ["functional_binder", "binder_combinatorial"]
    df_filtered = df_refolds[df_refolds['set'].isin(valid_types)].copy()
    grouped = df_filtered.groupby('refold_job_name')

    final_data = []

    for complex_id, group in tqdm(grouped, desc="Processing complexes"):
        status = group['set'].iloc[0]
        agonist_val = 1 if status == "functional_binder" else 0

        # Metadata
        anticonf_val = group['anticonf'].iloc[0]
        plddt_mean_val = group['plddt_mean'].iloc[0]
        flexibility_val = group['flexibility'].iloc[0]

        # 2. Initialize counts with INTEGER keys
        residue_counts = {res_id: 0 for res_id in clean_epitope_ids}

        pdb_files = ["../" + p for p in group['pdb_file'].tolist()]

        for pdb_path in pdb_files:
            try:
                results = count_target_residue_contacts(
                    file_path=pdb_path,
                    antigen_chain=antigen_chain,
                    antibody_chain=antibody_chain,
                    target_residues=target_set,
                    distance_threshold=4.5
                )

                # 3. Clean logic: result['target_residues_in_contact'] are already filtered
                # We just need to increment the counts
                for res_id in results['target_residues_in_contact']:
                    clean_id = int(res_id)
                    if clean_id in residue_counts:
                        residue_counts[clean_id] += 1

            except Exception as e:
                print(f"Error processing {pdb_path}: {e}")

        # Assemble row
        row = {
            "complex_id": complex_id,
            "agonist": agonist_val,
            "anticonf": anticonf_val,
            "plddt_mean": plddt_mean_val,
            "flexibility": flexibility_val
        }
        row.update(residue_counts)
        final_data.append(row)

    output_df = pd.DataFrame(final_data)

    # 4. Standardize Column Names
    meta_cols = ["complex_id", "agonist", "anticonf", "plddt_mean", "flexibility"]

    # Map all numeric columns to ints to avoid string vs int header issues
    new_cols = {col: int(col) for col in output_df.columns if col not in meta_cols}
    output_df.rename(columns=new_cols, inplace=True)

    # Sort the residue columns numerically
    sorted_res_cols = sorted([int(i) for i in clean_epitope_ids])

    return output_df[meta_cols + sorted_res_cols]