import yaml
import os
import pandas as pd

# Load the ligands
df = pd.read_csv("../data/CXCR4_chembl_ligands.csv", sep=";")
df = df[["Smiles", "Assay Type", "Molecule Name", "Molecular Weight", "AlogP", "#RO5 Violations", "Standard Type", "Standard Relation", "Standard Value", "Standard Units"]]
smiles = df["Smiles"].tolist()
# smiles = smiles[:3]


# Configuration
target_sequence = "MEGISIYTSDNYTEEMGSGDYDSMKEPCFREENANFNKIFLPTIYSIIFLTGIVGNGLVILVMGYQKKLRSMTDKYRLHLSVADLLFVITLPFWAVDAVANWYFGNFLCKAVHVIYTVNLYSSVLILAFISLDRYLAIVHATNSQRPRKLLAEKVVYVGVWIPALLLTIPDFIFANVSEADDRYICDRFYPNDLWVVVFQFQHIMVGLILPGIVILSCYCIIISKLSHSKGHQKRKALKTTVILILAFFACWLPYYIGISIDSFILLEIIKQGCEFENTVHKWISITEALAFFHCCLNPILYAFLGAKFKTSAQHALTSVSRGSSLKILSKGKRGGHSSVSTESESSSFHSS"
msa_path = "/home/senovad/boltz/boltz_results_cxcr4_initial_01/msa/cxcr4_initial_01_0.csv"
ligands = smiles

# Create yaml files
# os.makedirs("inputs_cxcr4", exist_ok=True)

for i, smiles in enumerate(ligands):
    data = {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": target_sequence, "msa": msa_path}},
            {"ligand": {"id": "B", "smiles": smiles}}
        ]
    }
    filename = f"../inputs_cxcr4/cxcr4_ligand_{i+1}.yaml"
    with open(filename, 'w') as f:
        yaml.dump(data, f)