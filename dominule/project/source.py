import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import RDLogger
from rdkit.Chem import Draw


lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)


def get_murcko(mol):
    """ Gets murcko scaffold for molecule rdkit object. """
    if mol:
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            return Chem.MolToSmiles(scaffold)
        except:
            print("Error in get_murcko: mol is None")
            return None
    print("Error in get_murcko: mol is None")
    return None

def get_generic(mol):
    """ Gets generic scaffold for molecule rdkit object. """
    if mol:
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            generic = MurckoScaffold.MakeScaffoldGeneric(scaffold)
            return Chem.MolToSmiles(generic)
        except:
            print("Error in get_generic: mol is None")
            return None
    return None


def convert_ugml_to_nm(df, ugml_col, mw_col):
    """
    Converts concentration from ug/mL to nM.
    Formula: (ug/mL / MW) * 1,000,000
    """
    # Calculation: (ug/mL / MW) gives Molar, then * 1e6 for nanomolar
    df['nM'] = (df[ugml_col] / df[mw_col]) * 1_000_000

    return df

# Example Usage:
# df = convert_ugml_to_nm(filtered, 'Standard Value', 'MW_column_name')


def draw_scaffold_comparison(df, row_id):
    """
    Draws the original molecule, its Murcko scaffold, and its generic scaffold.
    """
    row = df.iloc[row_id]

    # Get SMILES strings
    smiles_list = [
        row['Smiles'],
        row['Scaffold'],
        row['Generic Scaffold']
    ]

    # Convert SMILES to Mol objects (filter out None values)
    mols = [Chem.MolFromSmiles(s) for s in smiles_list if s is not None]

    # Create labels
    labels = ["Original", "Murcko Scaffold", "Generic Scaffold"]

    # Draw the grid
    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=3,
        subImgSize=(300, 300),
        legends=labels[:len(mols)]
    )

    return img



def create_balanced_cxcr4_dataset(df, id_col='Smiles', val_col='Standard Value',
                                  rel_col='Standard Relation', mw_col='Molecular Weight',
                                  unit_col='Standard Units', margin_percent=0.1):

    """
    Creates a balanced dataset for CXCR4 inhibitors based on the provided DataFrame.

    :param id_col: Column name for unique compound identifiers (for duplicates searching). Default is 'Smiles'.
    :param margin_percent: Percentage of the equals pool to drop from the middle to create a clearance gap. Default is 10% (0.1).

    :return: Adjusted DataFrame with balanced classes and a new 'Inhibitor' column (1 for active, 0 for inactive).
    """

    df = df.copy()

    # Sanitize the Relation column
    df[rel_col] = df[rel_col].astype(str).str.replace("'", "").str.replace('"', "").str.strip()

    # 1. Direct Conversion for ug.mL-1
    ugml_mask = df[unit_col].astype(str).str.contains(r'ug\.mL-1', case=False, na=False)
    if ugml_mask.any():
        vals = pd.to_numeric(df.loc[ugml_mask, val_col], errors='coerce')
        mws = pd.to_numeric(df.loc[ugml_mask, mw_col], errors='coerce')
        df.loc[ugml_mask, val_col] = (vals / mws) * 1_000_000
        df.loc[ugml_mask, unit_col] = 'nM'

    df['nM'] = pd.to_numeric(df[val_col], errors='coerce')
    df = df.dropna(subset=['nM', id_col]).reset_index(drop=True)

    # 2. Duplicate Handling & Consistency Check
    stats = {'averaged': 0, 'deleted': 0}

    def process_duplicates(group):
        if len(group) == 1:
            return group

        # Calculate log10 difference
        min_val = group['nM'].min()
        max_val = group['nM'].max()

        # Using 1e-9 to avoid log(0)
        log_diff = np.log10(max_val + 1e-9) - np.log10(min_val + 1e-9)

        if log_diff > 1.0:
            # Difference > 1 log unit: return empty to delete compound
            stats['deleted'] += 1
            return pd.DataFrame()
        else:
            # Difference <= 1 log unit: take the average
            stats['averaged'] += 1
            first_row = group.iloc[0].copy()
            first_row['nM'] = group['nM'].mean()
            # If relations conflict within duplicates, we treat the averaged result as '='
            if group[rel_col].nunique() > 1:
                first_row[rel_col] = '='
            return pd.DataFrame([first_row])

    print("-" * 40)
    print("CONSISTENCY & DUPLICATE ANALYSIS")
    print("-" * 40)
    print(f"Initial number of rows: {len(df)}")

    df = df.groupby(id_col, group_keys=False).apply(process_duplicates).reset_index(drop=True)

    print(f"Compounds averaged (diff <= 1 log unit): {stats['averaged']}")
    print(f"Compounds deleted  (diff > 1 log unit):  {stats['deleted']}")
    print(f"Rows remaining after cleanup:      {len(df)}")

    # 3. Categorize
    actives_strict = df[df[rel_col] == '<'].copy()
    inactives_strict = df[df[rel_col] == '>'].copy()
    equals_pool = df[df[rel_col] == '='].sort_values(by='nM').reset_index(drop=True)

    actives_strict['Inhibitor'] = 1
    inactives_strict['Inhibitor'] = 0

    # 4. Balancing with Expanded Margin
    total_potential = len(actives_strict) + len(inactives_strict) + len(equals_pool)
    target_per_class = total_potential // 2

    needed_act = max(0, target_per_class - len(actives_strict))
    needed_inact = max(0, target_per_class - len(inactives_strict))

    drop_buffer = int(len(equals_pool) * margin_percent)
    reduction = drop_buffer // 2
    final_needed_act = max(0, needed_act - reduction)
    final_needed_inact = max(0, needed_inact - reduction)

    eq_act = equals_pool.iloc[:final_needed_act].copy()
    eq_act['Inhibitor'] = 1

    if final_needed_inact > 0:
        eq_inact = equals_pool.iloc[-final_needed_inact:].copy()
        eq_inact['Inhibitor'] = 0
    else:
        eq_inact = pd.DataFrame(columns=df.columns)

    # --- FINAL STATS ---
    print("-" * 40)
    print("FINAL DATASET SUMMARY")
    print("-" * 40)
    if not eq_act.empty and not eq_inact.empty:
        print(f"Active Range:      {eq_act['nM'].min():.2f} - {eq_act['nM'].max():.2f} nM")
        print(f"Inactive Range:    {eq_inact['nM'].min():.2f} - {eq_inact['nM'].max():.2f} nM")
        print(f"Clearance Gap:     {eq_inact['nM'].min() - eq_act['nM'].max():.2f} nM")
        print(f"Dropped Margin:    {drop_buffer} rows")

    combined = pd.concat([actives_strict, inactives_strict, eq_act, eq_inact], ignore_index=True)
    counts = combined['Inhibitor'].value_counts()

    if len(counts) < 2:
        return combined

    min_class_size = int(counts.min())
    final_df = pd.concat([
        combined[combined['Inhibitor'] == 1].sample(n=min_class_size, random_state=42),
        combined[combined['Inhibitor'] == 0].sample(n=min_class_size, random_state=42)
    ]).sample(frac=1, random_state=42).reset_index(drop=True)

    print(f"Balanced Classes:  {min_class_size} Actives / {min_class_size} Inactives")
    print("-" * 40)

    return final_df