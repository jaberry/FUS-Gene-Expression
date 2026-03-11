"""
Cell-type marker genes and enrichment analysis.

Canonical marker genes for major brain cell types, used to assess
co-expression / spatial correlation with mechanosensitive channels.

Marker genes are based on well-established transcriptomic cell-type signatures
from studies including Tasic et al. 2018, Hodge et al. 2019, and the
Allen Cell Types Database.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

# ── cell type marker genes ─────────────────────────────────────────────────

CELL_TYPE_MARKERS = {
    "PV (Parvalbumin)": [
        "PVALB",     # definitive marker
        "GAD1",      # GABAergic
        "GAD2",
        "SST",       # included as negative control (marks different interneuron)
    ],
    "SST (Somatostatin)": [
        "SST",
        "GAD1",
        "GAD2",
        "NPY",
    ],
    "VIP": [
        "VIP",
        "GAD1",
        "GAD2",
        "CALB2",     # calretinin
    ],
    "Pyramidal / Excitatory": [
        "SLC17A7",   # VGLUT1 — cortical excitatory
        "SLC17A6",   # VGLUT2 — subcortical excitatory
        "CAMK2A",    # CaMKII-alpha
        "SATB2",     # upper layer cortical
        "TBR1",      # deep layer cortical
    ],
    "Astrocyte": [
        "GFAP",
        "AQP4",
        "ALDH1L1",
        "SLC1A2",    # GLT-1 / EAAT2
        "SLC1A3",    # GLAST / EAAT1
        "S100B",
    ],
    "Oligodendrocyte": [
        "MBP",
        "MOG",
        "OLIG2",
        "PLP1",
    ],
    "Microglia": [
        "CX3CR1",
        "P2RY12",
        "TMEM119",
        "AIF1",      # IBA1
    ],
}

# Flat list of all marker genes
ALL_MARKER_GENES = sorted(set(
    g for genes in CELL_TYPE_MARKERS.values() for g in genes
))


def compute_spatial_correlation(
    channel_expr: pd.DataFrame,
    marker_expr: pd.DataFrame,
    method: str = "spearman",
) -> pd.DataFrame:
    """Compute correlation between channel genes and cell-type markers across brain regions.

    Parameters
    ----------
    channel_expr : DataFrame
        Mechanosensitive channel genes (rows) x brain regions (columns).
    marker_expr : DataFrame
        Cell-type marker genes (rows) x brain regions (columns).
    method : str
        'spearman' or 'pearson'.

    Returns
    -------
    DataFrame of correlation coefficients (channel genes x marker genes).
    """
    # Align columns (shared brain regions)
    shared = channel_expr.columns.intersection(marker_expr.columns)
    ch = channel_expr[shared]
    mk = marker_expr[shared]

    n_channels = ch.shape[0]
    n_markers = mk.shape[0]

    corr_matrix = np.zeros((n_channels, n_markers))
    pval_matrix = np.zeros((n_channels, n_markers))

    for i, (cg, crow) in enumerate(ch.iterrows()):
        for j, (mg, mrow) in enumerate(mk.iterrows()):
            valid = ~(np.isnan(crow.values) | np.isnan(mrow.values))
            if valid.sum() < 5:
                corr_matrix[i, j] = np.nan
                pval_matrix[i, j] = np.nan
            else:
                if method == "spearman":
                    r, p = spearmanr(crow.values[valid], mrow.values[valid])
                else:
                    from scipy.stats import pearsonr
                    r, p = pearsonr(crow.values[valid], mrow.values[valid])
                corr_matrix[i, j] = r
                pval_matrix[i, j] = p

    corr_df = pd.DataFrame(corr_matrix, index=ch.index, columns=mk.index)
    pval_df = pd.DataFrame(pval_matrix, index=ch.index, columns=mk.index)

    return corr_df, pval_df


def compute_celltype_enrichment_scores(
    channel_expr: pd.DataFrame,
    marker_expr: pd.DataFrame,
    method: str = "spearman",
) -> pd.DataFrame:
    """Compute mean spatial correlation of each channel gene with each cell type.

    Returns DataFrame (channel genes x cell types) with mean correlation.
    """
    corr_df, pval_df = compute_spatial_correlation(channel_expr, marker_expr, method)

    enrichment = {}
    for cell_type, markers in CELL_TYPE_MARKERS.items():
        present = [m for m in markers if m in corr_df.columns]
        if present:
            enrichment[cell_type] = corr_df[present].mean(axis=1)

    return pd.DataFrame(enrichment)


def plot_celltype_enrichment(
    enrichment_df: pd.DataFrame,
    title: str = "Cell-Type Enrichment of Mechanosensitive Channels",
    figsize: tuple = (12, 10),
    cmap: str = "RdBu_r",
) -> tuple[plt.Figure, plt.Axes]:
    """Heatmap of channel gene enrichment across cell types."""
    from . import gene_lists

    # Rename rows
    display_names = [gene_lists.get_display_name(g) for g in enrichment_df.index]

    plot_df = enrichment_df.copy()
    plot_df.index = display_names

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        plot_df,
        cmap=cmap,
        center=0,
        annot=True,
        fmt=".2f",
        linewidths=0.5,
        ax=ax,
        cbar_kws={"label": "Mean Spearman ρ", "shrink": 0.6},
    )
    ax.set_title(title, fontsize=13, pad=15)
    ax.set_xlabel("Cell Type", fontsize=11)
    ax.set_ylabel("Mechanosensitive Channel Gene", fontsize=11)
    ax.tick_params(axis="x", rotation=30, labelsize=10)
    ax.tick_params(axis="y", labelsize=9)
    plt.tight_layout()

    return fig, ax


def plot_correlation_detail(
    corr_df: pd.DataFrame,
    pval_df: pd.DataFrame,
    alpha: float = 0.05,
    title: str = "Channel–Marker Spatial Correlation",
    figsize: tuple = (16, 12),
) -> tuple[plt.Figure, plt.Axes]:
    """Detailed heatmap of channel × marker gene correlations with significance."""
    from . import gene_lists

    display_names = [gene_lists.get_display_name(g) for g in corr_df.index]
    plot_df = corr_df.copy()
    plot_df.index = display_names

    # Create annotation: show rho with asterisk if significant
    annot = plot_df.copy().astype(str)
    for i in range(pval_df.shape[0]):
        for j in range(pval_df.shape[1]):
            val = corr_df.iloc[i, j]
            p = pval_df.iloc[i, j]
            if np.isnan(val):
                annot.iloc[i, j] = ""
            elif p < 0.001:
                annot.iloc[i, j] = f"{val:.2f}***"
            elif p < 0.01:
                annot.iloc[i, j] = f"{val:.2f}**"
            elif p < alpha:
                annot.iloc[i, j] = f"{val:.2f}*"
            else:
                annot.iloc[i, j] = f"{val:.2f}"

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        plot_df,
        cmap="RdBu_r",
        center=0,
        annot=annot,
        fmt="",
        linewidths=0.3,
        ax=ax,
        cbar_kws={"label": "Spearman ρ", "shrink": 0.5},
        annot_kws={"fontsize": 7},
    )
    ax.set_title(title, fontsize=13, pad=15)
    ax.set_xlabel("Cell-Type Marker Gene", fontsize=11)
    ax.set_ylabel("Mechanosensitive Channel", fontsize=11)
    ax.tick_params(axis="x", rotation=45, labelsize=8)
    ax.tick_params(axis="y", labelsize=9)
    plt.tight_layout()

    return fig, ax
