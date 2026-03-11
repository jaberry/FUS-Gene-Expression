"""
Hierarchical clustering of brain regions by mechanosensitive ion channel expression.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler

from . import gene_lists


def zscore_expression(df: pd.DataFrame) -> pd.DataFrame:
    """Z-score normalize expression values per gene (row-wise)."""
    scaler = StandardScaler()
    scaled = scaler.fit_transform(df.values.T).T
    return pd.DataFrame(scaled, index=df.index, columns=df.columns)


def build_family_palette(genes: list[str]) -> dict[str, str]:
    """Map gene family labels to colors for annotation."""
    family_colors = {
        "Piezo": "#E74C3C",
        "K2P (mechano)": "#3498DB",
        "K2P (other)": "#85C1E9",
        "TRPV": "#E67E22",
        "TRPA": "#F39C12",
        "TRPC": "#27AE60",
        "TRPM": "#1ABC9C",
        "TRPP": "#8E44AD",
        "Unknown": "#95A5A6",
    }
    return family_colors


def build_row_colors(genes: list[str]) -> pd.Series:
    """Create a Series of colors for gene family annotation on heatmaps."""
    palette = build_family_palette(genes)
    families = [gene_lists.get_family(g) for g in genes]
    return pd.Series([palette.get(f, "#95A5A6") for f in families], index=genes, name="Family")


def clustermap(
    expression_df: pd.DataFrame,
    title: str = "Mechanosensitive Ion Channels — Brain Region Clustering",
    figsize: tuple = (18, 10),
    method: str = "ward",
    metric: str = "euclidean",
    z_score_rows: bool = True,
    row_family_colors: bool = True,
    cmap: str = "RdBu_r",
    **kwargs,
) -> sns.matrix.ClusterGrid:
    """Generate a clustered heatmap of gene expression across brain regions.

    Parameters
    ----------
    expression_df : DataFrame
        Genes (rows) x brain regions (columns).
    title : str
        Plot title.
    method, metric : str
        Linkage method and distance metric for clustering.
    z_score_rows : bool
        If True, z-score normalize per gene before clustering.
    row_family_colors : bool
        If True, add colored sidebar showing ion channel family.

    Returns
    -------
    seaborn ClusterGrid object.
    """
    df = expression_df.copy()

    # Drop columns with all NaN
    df = df.dropna(axis=1, how="all")
    # Drop rows with all NaN
    df = df.dropna(axis=0, how="all")
    # Fill remaining NaN with row mean
    df = df.T.fillna(df.mean(axis=1)).T

    if z_score_rows:
        df = zscore_expression(df)

    # Rename rows to display names
    display_names = [gene_lists.get_display_name(g) for g in df.index]
    df.index = display_names

    row_colors = None
    if row_family_colors:
        families = [gene_lists.get_family(g) for g in expression_df.index if g in expression_df.index]
        palette = build_family_palette(list(expression_df.index))
        row_colors_list = [palette.get(f, "#95A5A6") for f in families]
        row_colors = pd.Series(row_colors_list, index=display_names, name="Family")

    g = sns.clustermap(
        df,
        method=method,
        metric=metric,
        row_colors=row_colors,
        cmap=cmap,
        figsize=figsize,
        linewidths=0.2,
        xticklabels=True,
        yticklabels=True,
        dendrogram_ratio=(0.15, 0.15),
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        **kwargs,
    )

    g.ax_heatmap.set_xlabel("Brain Region", fontsize=12)
    g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.tick_params(axis="x", labelsize=7, rotation=90)
    g.ax_heatmap.tick_params(axis="y", labelsize=9)
    g.fig.suptitle(title, fontsize=14, y=1.02)

    # Add family legend
    if row_family_colors:
        unique_families = sorted(set(families))
        legend_patches = [
            Patch(facecolor=palette[f], label=f) for f in unique_families
        ]
        g.ax_heatmap.legend(
            handles=legend_patches,
            title="Channel Family",
            bbox_to_anchor=(1.15, 1),
            loc="upper left",
            fontsize=8,
            title_fontsize=9,
        )

    return g


def region_dendrogram(
    expression_df: pd.DataFrame,
    method: str = "ward",
    metric: str = "euclidean",
    z_score_rows: bool = True,
    n_clusters: int | None = None,
    figsize: tuple = (16, 8),
    title: str = "Brain Region Dendrogram — Mechanosensitive Channel Expression",
) -> tuple[plt.Figure, np.ndarray | None]:
    """Plot a standalone dendrogram of brain regions.

    Returns (fig, cluster_labels) where cluster_labels is None if n_clusters is None.
    """
    df = expression_df.copy().dropna(axis=1, how="all").dropna(axis=0, how="all")
    df = df.T.fillna(df.mean(axis=1)).T

    if z_score_rows:
        df = zscore_expression(df)

    # Cluster on columns (regions)
    dist = pdist(df.values.T, metric=metric)
    Z = linkage(dist, method=method)

    fig, ax = plt.subplots(figsize=figsize)

    color_threshold = None
    if n_clusters is not None:
        # Find threshold for desired cluster count
        from scipy.cluster.hierarchy import maxdists
        thresholds = sorted(Z[:, 2])
        color_threshold = thresholds[-(n_clusters - 1)] if n_clusters <= len(thresholds) else 0

    dendrogram(
        Z,
        labels=df.columns.tolist(),
        leaf_rotation=90,
        leaf_font_size=8,
        ax=ax,
        color_threshold=color_threshold,
    )
    ax.set_title(title, fontsize=13)
    ax.set_ylabel(f"Distance ({metric})", fontsize=11)
    ax.set_xlabel("Brain Region", fontsize=11)
    plt.tight_layout()

    cluster_labels = None
    if n_clusters is not None:
        cluster_labels = fcluster(Z, n_clusters, criterion="maxclust")

    return fig, cluster_labels
