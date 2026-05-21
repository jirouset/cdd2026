import pandas as pd

# Paths
original_csv = "../data/CXCR4_chembl_ligands.csv"
output_csv = "../data/CXCR4_chembl_ligands_boltz.csv"

# Load original data
df = pd.read_csv(original_csv, sep=";")

# Create an empty list to hold the expanded rows
expanded_rows = []

print("Mapping Boltz file formats...")

# Iterate through each row in the original dataframe
for index, row in df.iterrows():
    # Mapping logic: ligand_id is 1-based index from prepare_yaml_files.py
    ligand_id = index + 1

    # Generate entries for all 5 runs (0 to 4) produced by submit_boltz_array.sh
    for run_id in range(5):
        cif_file = f"cxcr4_ligand_{ligand_id}_run_{run_id}_model_0.cif"
        json_file = f"confidence_cxcr4_ligand_{ligand_id}_run_{run_id}_model_0.json"

        # Copy the original row data and add the new columns
        new_row = row.copy()
        new_row['Ligand_ID'] = ligand_id
        new_row['Run_ID'] = run_id
        new_row['boltz_cif_file'] = cif_file
        new_row['boltz_confidence_file'] = json_file

        expanded_rows.append(new_row)

# Create the new expanded dataframe
output_df = pd.DataFrame(expanded_rows)

# Reorder columns slightly to keep IDs at the front
cols = ['Ligand_ID', 'Run_ID', 'boltz_cif_file', 'boltz_confidence_file'] + [c for c in df.columns]
output_df = output_df[cols]

# Save updated dataframe
output_df.to_csv(output_csv, index=False, sep=";")
print(f"Updated CSV saved to: {output_csv}")