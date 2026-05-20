"""Analysis functions for CXCR4 Boltz cofolding contact data."""

import json
import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from rdkit import Chem
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

sns.set_style("whitegrid")

PALETTE = {"Active": "#e74c3c", "Inactive": "#3498db"}

_AMIDE_SMARTS = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")


def _is_peptide(smiles: str, amide_threshold: int = 3) -> bool:
    """Return True if the molecule has >= amide_threshold amide bonds (peptide heuristic)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    return len(mol.GetSubstructMatches(_AMIDE_SMARTS)) >= amide_threshold


def prepare_valid_df(df_contacts: pd.DataFrame):
    """Filter to labelled ligands with ≥1 successful run and drop all-zero residue columns.

    Returns (df_valid, res_cols).
    """
    res_cols_all = [c for c in df_contacts.columns if isinstance(c, int)]

    nonzero_cols = [c for c in res_cols_all if df_contacts[c].sum() > 0]
    n_dropped = len(res_cols_all) - len(nonzero_cols)
    if n_dropped:
        print(f"Dropped {n_dropped} all-zero residue columns. Remaining: {len(nonzero_cols)}")

    keep = ["Ligand_ID", "Inhibitor", "contact_count_mean", "runs_processed"] + nonzero_cols
    df_valid = df_contacts[
        df_contacts["Inhibitor"].notna() & (df_contacts["runs_processed"] > 0)
    ][keep].copy()

    df_valid["Group"] = df_valid["Inhibitor"].map({1.0: "Active", 0.0: "Inactive"})

    n_act = (df_valid["Inhibitor"] == 1.0).sum()
    n_ina = (df_valid["Inhibitor"] == 0.0).sum()
    print(f"Ligands with runs: {len(df_valid)}  (actives={n_act}, inactives={n_ina})")

    return df_valid, nonzero_cols


def plot_contact_count_boxplot(df_valid: pd.DataFrame, compare=False):
    """
    Boxplot of mean contact count per ligand, actives vs inactives.

    :param compare: if True, compare actives and inactives with Mann-Whitney U test (enable only if the general shape is similar and other assumptions are met).
    """

    act = df_valid[df_valid["Group"] == "Active"]["contact_count_mean"]
    ina = df_valid[df_valid["Group"] == "Inactive"]["contact_count_mean"]
    u_stat, p_val = mannwhitneyu(act, ina, alternative="two-sided")

    fig, ax = plt.subplots(figsize=(5, 5))
    sns.boxplot(data=df_valid, x="Group", y="contact_count_mean",
                hue="Group", palette=PALETTE, width=0.5, fliersize=0, legend=False, ax=ax)
    sns.stripplot(data=df_valid, x="Group", y="contact_count_mean",
                  color="black", alpha=0.35, size=3.5, jitter=True, ax=ax)
    ax.set_ylabel("Average count of contacts of ligand")
    ax.set_xlabel("")
    if compare:
        ax.set_title(f"Protein–Ligand Contact Count\nMann–Whitney U = {u_stat:.0f},  p = {p_val:.3g}")
    else:
        ax.set_title(f"Protein–Ligand Contact Count")
    plt.tight_layout()
    plt.show()
    print(f"Median — Active: {act.median():.1f},  Inactive: {ina.median():.1f}")


def plot_contact_count_boxplot_peptides(
    df_valid: pd.DataFrame,
    df_smiles: pd.DataFrame,
    smiles_col: str = "Smiles",
    id_col: str = "Ligand_ID",
    amide_threshold: int = 3,
):
    """Same boxplot as plot_contact_count_boxplot, but peptide-based ligands are highlighted.

    Peptides are identified by having >= amide_threshold amide bonds in their SMILES.
    In the strip overlay, peptides are shown as larger orange markers; non-peptides as small grey dots.

    :param df_smiles: DataFrame containing at least id_col and smiles_col (e.g. boltz CSV deduplicated).
    :param amide_threshold: minimum number of amide bonds to classify as peptide (default 3).
    """
    # Join SMILES onto df_valid
    smiles_map = (
        df_smiles[[id_col, smiles_col]]
        .drop_duplicates(subset=id_col)
        .set_index(id_col)[smiles_col]
    )
    df_plot = df_valid.copy()
    df_plot["Smiles"] = df_plot["Ligand_ID"].map(smiles_map)
    df_plot["Peptide"] = df_plot["Smiles"].apply(
        lambda s: _is_peptide(s, amide_threshold) if pd.notna(s) else False
    )

    n_peptides = df_plot["Peptide"].sum()

    act = df_plot[df_plot["Group"] == "Active"]["contact_count_mean"]
    ina = df_plot[df_plot["Group"] == "Inactive"]["contact_count_mean"]

    fig, ax = plt.subplots(figsize=(5, 5))
    sns.boxplot(
        data=df_plot, x="Group", y="contact_count_mean",
        hue="Group", palette=PALETTE, width=0.5, fliersize=0, legend=False, ax=ax,
    )

    # Non-peptides: small grey dots
    non_pep = df_plot[~df_plot["Peptide"]]
    sns.stripplot(
        data=non_pep, x="Group", y="contact_count_mean",
        color="black", alpha=0.25, size=3, jitter=True, ax=ax,
    )

    # Peptides: larger orange markers on top
    pep = df_plot[df_plot["Peptide"]]
    sns.stripplot(
        data=pep, x="Group", y="contact_count_mean",
        color="#e67e22", alpha=0.85, size=6, jitter=True, marker="D", ax=ax,
    )

    ax.set_ylabel("Average count of contacts of ligand")
    ax.set_xlabel("")
    ax.set_title(
        f"Protein–Ligand Contact Count (peptides highlighted)\n"
        f"Peptides: {n_peptides} / {len(df_plot)}"
    )

    legend_handles = [
        mpatches.Patch(color="#e74c3c", label="Active"),
        mpatches.Patch(color="#3498db", label="Inactive"),
        plt.Line2D([0], [0], marker="D", color="w", markerfacecolor="#e67e22",
                   markersize=7, label=f"Peptide-based (≥{amide_threshold} amide bonds)"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="black",
                   markersize=5, alpha=0.5, label="Small molecule"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=8)
    plt.tight_layout()
    plt.show()
    print(f"Peptide-based ligands: {n_peptides}  ({100 * n_peptides / len(df_plot):.1f}%)")
    print(f"  Active peptides  : {(df_plot['Peptide'] & (df_plot['Inhibitor'] == 1.0)).sum()}")
    print(f"  Inactive peptides: {(df_plot['Peptide'] & (df_plot['Inhibitor'] == 0.0)).sum()}")


def plot_contact_frequency_heatmap(df_valid: pd.DataFrame, res_cols: list, min_freq: float = 0.40):
    """Heatmap of mean contact frequency per residue for each group.

    Only residues with mean contact frequency ≥ min_freq in at least one group are shown.
    """
    df_freq = df_valid[res_cols].div(df_valid["runs_processed"], axis=0)
    df_freq["Group"] = df_valid["Group"].values
    heatmap_data = df_freq.groupby("Group")[res_cols].mean()

    top_res = sorted(heatmap_data.max(axis=0).pipe(lambda s: s[s >= min_freq]).index.tolist())

    fig, ax = plt.subplots(figsize=(max(12, len(top_res) * 0.35), 3))
    sns.heatmap(heatmap_data[top_res], cmap="YlOrRd", vmin=0, vmax=1,
                xticklabels=True, yticklabels=True,
                linewidths=0.3, linecolor="white",
                cbar_kws={"label": "Mean fraction of runs in contact", "shrink": 0.8},
                ax=ax)
    ax.set_title("CXCR4 Residue Contact Frequency — Actives vs Inactives")
    ax.set_xlabel("CXCR4 Residue ID")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelsize=8, rotation=90)
    plt.tight_layout()
    plt.show()
    print(f"Residues shown: {len(top_res)} of {len(res_cols)}  (threshold ≥ {min_freq})")


def run_mannwhitney_test(df_valid: pd.DataFrame, res_cols: list) -> pd.DataFrame:
    """Per-residue Mann-Whitney U test between actives and inactives.

    Y metric: contact frequency in percentage points
        freq_pct = sum(contact counts) / sum(runs_processed) * 100

    This normalises for ligands with fewer than 5 runs.

    Returns DataFrame sorted by p_value:
        residue | U_stat | p_value | p_fdr | freq_active_pct | freq_inactive_pct | diff_pct
    """
    actives_df = df_valid[df_valid["Inhibitor"] == 1.0]
    inactives_df = df_valid[df_valid["Inhibitor"] == 0.0]
    active_total_runs = actives_df["runs_processed"].sum()
    inactive_total_runs = inactives_df["runs_processed"].sum()

    rows = []
    for res in res_cols:
        a_vals = actives_df[res].values
        i_vals = inactives_df[res].values
        if a_vals.sum() == 0 and i_vals.sum() == 0:
            continue
        stat, p = mannwhitneyu(a_vals, i_vals, alternative="two-sided")
        rows.append({
            "residue": res,
            "U_stat": stat,
            "p_value": p,
            "freq_active_pct": a_vals.sum() / active_total_runs * 100,
            "freq_inactive_pct": i_vals.sum() / inactive_total_runs * 100,
            "diff_pct": a_vals.sum() / active_total_runs * 100
                        - i_vals.sum() / inactive_total_runs * 100,
        })

    df_mw = pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)
    _, p_fdr, _, _ = multipletests(df_mw["p_value"], method="fdr_bh")
    df_mw["p_fdr"] = p_fdr

    print(f"Residues tested:              {len(df_mw)}")
    print(f"  p < 0.05 (uncorrected):     {(df_mw['p_value'] < 0.05).sum()}")
    print(f"  q < 0.05 (FDR-corrected):   {(df_mw['p_fdr'] < 0.05).sum()}")
    return df_mw


def _draw_residue_bars(ax, df_plot: pd.DataFrame, y_lim: tuple, title: str):
    """Draw a single residue bar chart sorted by residue ID on ax."""
    df_sorted = df_plot.sort_values("residue").reset_index(drop=True)
    colors = ["#e74c3c" if d > 0 else "#3498db" for d in df_sorted["diff_pct"]]

    bars = ax.bar(range(len(df_sorted)), df_sorted["diff_pct"], color=colors, width=0.8)

    for bar, (_, row) in zip(bars, df_sorted.iterrows()):
        star = "**" if row["p_fdr"] < 0.01 else ("*" if row["p_fdr"] < 0.05 else "")
        if star:
            offset = 0.3 if bar.get_height() >= 0 else -0.3
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + offset, star,
                    ha="center", va="bottom" if bar.get_height() >= 0 else "top",
                    fontsize=9, fontweight="bold")

    # Show every Nth label to avoid crowding on large charts
    tick_step = max(1, len(df_sorted) // 30)
    labels = [str(r) if i % tick_step == 0 else ""
              for i, r in enumerate(df_sorted["residue"])]
    ax.set_xticks(range(len(df_sorted)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylim(y_lim)
    ax.set_ylabel("Contact frequency difference (%)\n(Active − Inactive)")
    ax.set_xlabel("CXCR4 Residue ID")
    ax.set_title(title)


def plot_residue_contact_diff(df_mw: pd.DataFrame, p_threshold: float = 0.05):
    """Two bar charts of residue contact frequency difference (Active% − Inactive%).

    Chart 1 — all residues sorted by residue ID on x-axis.
    Chart 2 — only residues with p < p_threshold, same y-axis limits.

    Y-axis unit: percentage points (sum of contacts / sum of runs × 100).
    Red = actives contact more,  Blue = inactives contact more.
    Stars mark FDR significance:  * q < 0.05,  ** q < 0.01.
    """
    sig = df_mw[df_mw["p_value"] < p_threshold].copy()
    y_max = df_mw["diff_pct"].abs().max() * 1.15
    y_lim = (-y_max, y_max)

    legend = [
        mpatches.Patch(color="#e74c3c", label="Active > Inactive"),
        mpatches.Patch(color="#3498db", label="Inactive > Active"),
    ]

    fig1, ax1 = plt.subplots(figsize=(max(14, len(df_mw) * 0.12), 5))
    _draw_residue_bars(ax1, df_mw, y_lim,
                       "All CXCR4 Residues — Contact Frequency Difference\n"
                       "* q < 0.05   ** q < 0.01  (FDR-corrected)")
    ax1.legend(handles=legend, loc="upper right")
    plt.tight_layout()
    plt.show()

    if len(sig) == 0:
        print(f"No residues with p < {p_threshold}.")
        return

    fig2, ax2 = plt.subplots(figsize=(max(8, len(sig) * 0.35), 5))
    _draw_residue_bars(ax2, sig, y_lim,
                       f"Significant Residues (p < {p_threshold}) — Contact Frequency Difference\n"
                       "* q < 0.05   ** q < 0.01  (FDR-corrected)")
    ax2.legend(handles=legend, loc="upper right")
    plt.tight_layout()
    plt.show()
    print(f"Significant residues: {len(sig)}")


# ── Confidence scores ─────────────────────────────────────────────────────────

CONFIDENCE_METRICS = {
    "confidence_score": "Confidence score",
    "iptm":             "ipTM  (interface pTM)",
    "ligand_iptm":      "Ligand ipTM",
    "complex_plddt":    "Complex pLDDT",
}


def load_confidence_scores(
    df_valid: pd.DataFrame,
    confidence_dir: str,
    n_runs: int = 5,
) -> pd.DataFrame:
    """Load Boltz confidence JSON files and return per-ligand mean scores.

    Iterates over all run files for each ligand in df_valid (which already
    carries the Inhibitor label).  Missing files are skipped silently.

    Returns DataFrame with columns:
        Ligand_ID | Inhibitor | Group | confidence_score | iptm | ligand_iptm | complex_plddt
    """
    rows = []
    for _, lig_row in tqdm(df_valid.iterrows(), total=len(df_valid),
                           desc="Loading confidence scores"):
        ligand_id = int(lig_row["Ligand_ID"])
        run_scores = {k: [] for k in CONFIDENCE_METRICS}

        for run_id in range(n_runs):
            fname = f"confidence_cxcr4_ligand_{ligand_id}_run_{run_id}_model_0.json"
            fpath = os.path.join(confidence_dir, fname)
            if not os.path.isfile(fpath):
                continue
            with open(fpath) as f:
                data = json.load(f)
            for k in CONFIDENCE_METRICS:
                if k in data:
                    run_scores[k].append(data[k])

        if not any(run_scores[k] for k in CONFIDENCE_METRICS):
            continue

        row = {
            "Ligand_ID": ligand_id,
            "Inhibitor": lig_row["Inhibitor"],
            "Group": lig_row["Group"],
        }
        for k in CONFIDENCE_METRICS:
            row[k] = float(np.mean(run_scores[k])) if run_scores[k] else np.nan
        rows.append(row)

    df_conf = pd.DataFrame(rows)
    n_act = (df_conf["Inhibitor"] == 1.0).sum()
    n_ina = (df_conf["Inhibitor"] == 0.0).sum()
    print(f"Loaded confidence for {len(df_conf)} ligands  (actives={n_act}, inactives={n_ina})")
    return df_conf


def plot_confidence_boxplot(
    df_valid: pd.DataFrame,
    confidence_dir: str,
    n_runs: int = 5,
):
    """Boxplot of mean Boltz confidence metrics for actives vs inactives.

    Loads confidence JSON files from confidence_dir, averages across runs,
    and plots one panel per metric with a Mann-Whitney U p-value.

    Parameters
    ----------
    df_valid       : filtered contacts DataFrame from prepare_valid_df()
    confidence_dir : path to the folder containing confidence JSON files
    n_runs         : maximum number of runs per ligand to look for (default 5)
    """
    df_conf = load_confidence_scores(df_valid, confidence_dir, n_runs)

    metrics = list(CONFIDENCE_METRICS.keys())
    labels  = list(CONFIDENCE_METRICS.values())
    n_metrics = len(metrics)

    # Shared y limits across all metrics
    all_vals = df_conf[metrics].values
    y_min = float(np.nanmin(all_vals))
    y_max = float(np.nanmax(all_vals))
    margin = (y_max - y_min) * 0.05
    y_lim = (y_min - margin, y_max + margin)

    fig, axes = plt.subplots(1, n_metrics, figsize=(4 * n_metrics, 5), sharey=True)

    for ax, metric, label in zip(axes, metrics, labels):
        df_plot = df_conf[["Group", metric]].dropna()
        act = df_plot[df_plot["Group"] == "Active"][metric]
        ina = df_plot[df_plot["Group"] == "Inactive"][metric]

        if len(act) > 0 and len(ina) > 0:
            _, p = mannwhitneyu(act, ina, alternative="two-sided")
            p_str = f"p = {p:.3g}"
        else:
            p_str = "n/a"

        sns.boxplot(data=df_plot, x="Group", y=metric,
                    hue="Group", palette=PALETTE, width=0.5, fliersize=0, legend=False, ax=ax)
        sns.stripplot(data=df_plot, x="Group", y=metric,
                      color="black", alpha=0.25, size=3, jitter=True, ax=ax)
        ax.set_title(f"{label}\n{p_str}", fontsize=10)
        ax.set_ylim(y_lim)
        ax.set_xlabel("")
        ax.set_ylabel("")

    fig.suptitle("Boltz Confidence Scores — Actives vs Inactives", fontsize=12, y=1.02)
    plt.tight_layout()
    plt.show()


def plot_confidence_boxplot_peptides(
    df_valid: pd.DataFrame,
    df_smiles: pd.DataFrame,
    confidence_dir: str,
    smiles_col: str = "Smiles",
    id_col: str = "Ligand_ID",
    amide_threshold: int = 3,
    n_runs: int = 5,
):
    """Same 4-panel confidence boxplot as plot_confidence_boxplot, with peptides highlighted.

    Peptide-based ligands (>= amide_threshold amide bonds) are shown as orange diamonds
    in the strip overlay; non-peptides as small grey dots.

    :param df_smiles:       DataFrame with id_col and smiles_col (e.g. boltz CSV).
    :param confidence_dir:  Path to the folder containing confidence JSON files.
    :param amide_threshold: Minimum amide bonds to classify as peptide (default 3).
    :param n_runs:          Maximum runs per ligand to look for (default 5).
    """
    df_conf = load_confidence_scores(df_valid, confidence_dir, n_runs)

    # Attach peptide flag
    smiles_map = (
        df_smiles[[id_col, smiles_col]]
        .drop_duplicates(subset=id_col)
        .set_index(id_col)[smiles_col]
    )
    df_conf["Smiles"] = df_conf["Ligand_ID"].map(smiles_map)
    df_conf["Peptide"] = df_conf["Smiles"].apply(
        lambda s: _is_peptide(s, amide_threshold) if pd.notna(s) else False
    )

    metrics = list(CONFIDENCE_METRICS.keys())
    labels = list(CONFIDENCE_METRICS.values())
    n_metrics = len(metrics)

    all_vals = df_conf[metrics].values
    y_min = float(np.nanmin(all_vals))
    y_max = float(np.nanmax(all_vals))
    margin = (y_max - y_min) * 0.05
    y_lim = (y_min - margin, y_max + margin)

    fig, axes = plt.subplots(1, n_metrics, figsize=(4 * n_metrics, 5), sharey=True)

    for ax, metric, label in zip(axes, metrics, labels):
        df_plot = df_conf[["Group", "Peptide", metric]].dropna(subset=[metric])
        act = df_plot[df_plot["Group"] == "Active"][metric]
        ina = df_plot[df_plot["Group"] == "Inactive"][metric]

        p_str = f"p = {mannwhitneyu(act, ina, alternative='two-sided')[1]:.3g}" if (len(act) > 0 and len(ina) > 0) else "n/a"

        sns.boxplot(data=df_plot, x="Group", y=metric,
                    hue="Group", palette=PALETTE, width=0.5, fliersize=0, legend=False, ax=ax)

        # Non-peptides: small grey dots
        sns.stripplot(data=df_plot[~df_plot["Peptide"]], x="Group", y=metric,
                      color="black", alpha=0.2, size=3, jitter=True, ax=ax)

        # Peptides: orange diamonds
        sns.stripplot(data=df_plot[df_plot["Peptide"]], x="Group", y=metric,
                      color="#e67e22", alpha=0.85, size=6, jitter=True, marker="D", ax=ax)

        ax.set_title(f"{label}\n{p_str}", fontsize=10)
        ax.set_ylim(y_lim)
        ax.set_xlabel("")
        ax.set_ylabel("")

    n_pep = df_conf["Peptide"].sum()
    legend_handles = [
        mpatches.Patch(color="#e74c3c", label="Active"),
        mpatches.Patch(color="#3498db", label="Inactive"),
        plt.Line2D([0], [0], marker="D", color="w", markerfacecolor="#e67e22",
                   markersize=7, label=f"Peptide (≥{amide_threshold} amide bonds, n={n_pep})"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="black",
                   markersize=5, alpha=0.5, label="Small molecule"),
    ]
    axes[-1].legend(handles=legend_handles, loc="lower right", fontsize=8)

    fig.suptitle("Boltz Confidence Scores — Actives vs Inactives (peptides highlighted)",
                 fontsize=12, y=1.02)
    plt.tight_layout()
    plt.show()
