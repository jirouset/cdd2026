import pandas as pd
import io

# Loading your data
df = pd.read_csv("your_file.csv")

# 1. Melt the dataframe to turn columns into rows
# 2. Filter out rows where the ID value is missing (NaN or empty)
new_df = df.melt(
    id_vars=["Canonical_SMILES"],
    value_vars=["drugbank_id", "ChEMBL ID", "zincid"],
    var_name="Source",
    value_name="Source_ID"
).dropna(subset=["Source_ID"])

# Optional: If you only want the "Canonical_SMILES" and "Source" columns
result = new_df[["Canonical_SMILES", "Source"]]

print(result.head())