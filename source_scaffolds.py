import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Scaffolds import MurckoScaffold
from source import get_generic, get_murcko
from IPython.display import display


from sklearn.manifold import TSNE
import seaborn as sns
from sklearn.decomposition import PCA


def extract_morgan_fingerprints(df, smiles_col='Smiles', radius=2, n_bits=1024):
    """
    Generates Morgan (ECFP) fingerprints from SMILES.
    

    :return: dataframe (pd.DataFrame) and list of fingerprints column names (List[str])
    """
    
    def mol_to_fp(smiles):
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)
            return np.fromiter(fp.ToBitString(), dtype=np.int8)
        else:
            return np.zeros(n_bits, dtype=np.int8)


    fps_array = np.vstack(df[smiles_col].apply(mol_to_fp))

    fp_columns = [f'FP_{i}' for i in range(n_bits)]
    
    fp_df = pd.DataFrame(fps_array, columns=fp_columns, index=df.index)
    
    return fp_df, fp_columns

def get_scaffolds(df, smiles_col='Smiles'):
    df = df.copy()
    
    df['Mol'] = df[smiles_col].apply(Chem.MolFromSmiles)

    df['Scaffold'] = df['Mol'].apply(get_murcko)
    df['Generic Scaffold'] = df['Mol'].apply(get_generic)

    return df.drop(columns=['Mol']), ['Scaffold', 'Generic Scaffold']


def draw_scaffolds(df, scaff_col):
    top_scaffolds = df[scaff_col].value_counts().head(6)

    scaffold_mols = [Chem.MolFromSmiles(smiles) for smiles in top_scaffolds.index if smiles]
    
    labels = [f"Count: {count}" for count in top_scaffolds.values]

    img = Draw.MolsToGridImage(
        scaffold_mols, 
        molsPerRow=3, 
        subImgSize=(300, 200), 
        legends=labels
    )
    
    display(img)

    return


def compute_tsne_from_fingerprints(df, n_components=2, perplexity=30):
    """
    Counts t-SNE coords for fingerprints using Jaccardo (Tanimoto) distance.
    """
    X = df.values.astype(np.float32)
    
    tsne = TSNE(
        n_components=n_components, 
        metric='jaccard', 
        perplexity=perplexity
    )
    
    tsne_results = tsne.fit_transform(X)

    df_tsne = pd.DataFrame(
        tsne_results, 
        columns=['tSNE_1', 'tSNE_2'], 
        index=df.index
    )
    df_tsne['Inhibitor']=df['Inhibitor']
    return df_tsne


def compute_pca_from_fingerprints(df, n_components=2):
    """
    Counts PCA coordinates.
    """
    X = df.values.astype(np.float32)
    
    pca = PCA(n_components=n_components)
    pca_results = pca.fit_transform(X)

    df_pca = pd.DataFrame(
        pca_results, 
        columns=['PC_1', 'PC_2'], 
        index=df.index
    )
    df_pca['Inhibitor']=df['Inhibitor']
    
    return df_pca


def draw_clusters(df, method='pca', target_col='Inhibitor', n_components=2, perplexity=30, colors=['#440154', '#fde725']):
    """
    Counts PCA or TSNE from data a draw results.
    """
    if method == 'pca':
        df_result = compute_pca_from_fingerprints(df)
        title='PCA'
        x='PC_1'
        y='PC_2'
        
    else:
        df_result = compute_tsne_from_fingerprints(df, n_components=n_components, perplexity=perplexity)
        title='TSNE'
        x='tSNE_1'
        y='tSNE_2'
    

    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        x=x, 
        y=y, 
        hue=target_col,
        palette=colors,
        data=df_result,
        alpha=0.7
    )
    plt.title(title)
    plt.xlabel(x)
    plt.ylabel(y)
    plt.show()