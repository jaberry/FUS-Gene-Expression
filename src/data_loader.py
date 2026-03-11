"""
Load and preprocess Allen Human Brain Atlas microarray and RNA-seq datasets.

Expects raw zip files in the data directory with the naming convention:
    normalized_microarray_donor{id}.zip
    rnaseq_donor{id}.zip
"""
from __future__ import annotations

import zipfile
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

DATA_DIR = Path("/Users/Jack/Documents/UCSF/Research/Gene Expression")

MICROARRAY_DONORS = [9861, 10021, 12876, 14380, 15496, 15697]
RNASEQ_DONORS = [9861, 10021]


# ── helpers ────────────────────────────────────────────────────────────────

def _read_csv_from_zip(zip_path: Path, filename: str, **kwargs) -> pd.DataFrame:
    with zipfile.ZipFile(zip_path) as zf:
        with zf.open(filename) as f:
            return pd.read_csv(f, **kwargs)


# ── microarray ─────────────────────────────────────────────────────────────

def load_microarray_probes(donor_id: int = 9861) -> pd.DataFrame:
    """Load probe metadata for a microarray donor."""
    path = DATA_DIR / f"normalized_microarray_donor{donor_id}.zip"
    probes = _read_csv_from_zip(path, "Probes.csv")
    return probes


def load_microarray_samples(donor_id: int = 9861) -> pd.DataFrame:
    """Load sample annotations (brain region metadata) for a microarray donor."""
    path = DATA_DIR / f"normalized_microarray_donor{donor_id}.zip"
    return _read_csv_from_zip(path, "SampleAnnot.csv")


def load_microarray_expression(donor_id: int = 9861) -> pd.DataFrame:
    """Load full microarray expression matrix (probes x samples).

    Returns DataFrame with probe_id as index. Column order matches SampleAnnot row order.
    """
    path = DATA_DIR / f"normalized_microarray_donor{donor_id}.zip"
    expr = _read_csv_from_zip(path, "MicroarrayExpression.csv", header=None)
    expr = expr.set_index(0)
    expr.index.name = "probe_id"
    return expr


def load_microarray_ontology(donor_id: int = 9861) -> pd.DataFrame:
    """Load brain region ontology."""
    path = DATA_DIR / f"normalized_microarray_donor{donor_id}.zip"
    return _read_csv_from_zip(path, "Ontology.csv")


def load_microarray_pacall(donor_id: int = 9861) -> pd.DataFrame:
    """Load present/absent call matrix for microarray."""
    path = DATA_DIR / f"normalized_microarray_donor{donor_id}.zip"
    pac = _read_csv_from_zip(path, "PACall.csv", header=None)
    pac = pac.set_index(0)
    pac.index.name = "probe_id"
    return pac


def get_microarray_gene_expression(
    gene_symbols: list[str],
    donor_id: int = 9861,
    collapse_probes: str = "mean",
) -> pd.DataFrame:
    """Extract expression for specific genes from microarray data.

    Parameters
    ----------
    gene_symbols : list of str
        Gene symbols to extract (as they appear in Probes.csv).
    donor_id : int
        Donor ID.
    collapse_probes : str
        How to collapse multiple probes per gene: 'mean', 'median', or 'max'.

    Returns
    -------
    DataFrame with genes as rows and samples as columns.
    """
    probes = load_microarray_probes(donor_id)
    samples = load_microarray_samples(donor_id)
    expr = load_microarray_expression(donor_id)

    # filter probes for our genes
    mask = probes["gene_symbol"].isin(gene_symbols)
    probe_subset = probes.loc[mask, ["probe_id", "gene_symbol"]]

    if probe_subset.empty:
        return pd.DataFrame()

    # extract expression rows for these probes
    expr_subset = expr.loc[expr.index.isin(probe_subset["probe_id"])]

    # merge gene_symbol onto expression
    expr_subset = expr_subset.merge(
        probe_subset.set_index("probe_id")[["gene_symbol"]],
        left_index=True,
        right_index=True,
    )

    # collapse probes per gene
    agg_func = {"mean": "mean", "median": "median", "max": "max"}[collapse_probes]
    gene_expr = expr_subset.groupby("gene_symbol").agg(agg_func)

    # set column names to structure acronyms
    gene_expr.columns = samples["structure_acronym"].values[: gene_expr.shape[1]]

    return gene_expr


# ── RNA-seq ────────────────────────────────────────────────────────────────

def load_rnaseq_genes(donor_id: int = 9861) -> pd.DataFrame:
    """Load gene metadata for RNA-seq dataset."""
    path = DATA_DIR / f"rnaseq_donor{donor_id}.zip"
    return _read_csv_from_zip(path, "Genes.csv")


def load_rnaseq_samples(donor_id: int = 9861) -> pd.DataFrame:
    """Load sample annotations for RNA-seq dataset."""
    path = DATA_DIR / f"rnaseq_donor{donor_id}.zip"
    return _read_csv_from_zip(path, "SampleAnnot.csv")


def load_rnaseq_tpm(donor_id: int = 9861) -> pd.DataFrame:
    """Load RNA-seq TPM expression matrix (genes x samples).

    First column is gene_symbol used as row index.
    """
    path = DATA_DIR / f"rnaseq_donor{donor_id}.zip"
    tpm = _read_csv_from_zip(path, "RNAseqTPM.csv", header=None)
    tpm = tpm.set_index(0)
    tpm.index.name = "gene_symbol"
    return tpm


def load_rnaseq_ontology(donor_id: int = 9861) -> pd.DataFrame:
    """Load brain region ontology from RNA-seq zip."""
    path = DATA_DIR / f"rnaseq_donor{donor_id}.zip"
    return _read_csv_from_zip(path, "Ontology.csv")


def get_rnaseq_gene_expression(
    gene_symbols: list[str],
    donor_id: int = 9861,
) -> pd.DataFrame:
    """Extract TPM expression for specific genes from RNA-seq data.

    Returns DataFrame with genes as rows and brain region acronyms as columns.
    """
    tpm = load_rnaseq_tpm(donor_id)
    samples = load_rnaseq_samples(donor_id)

    mask = tpm.index.isin(gene_symbols)
    gene_expr = tpm.loc[mask].copy()

    # set column names to structure acronyms
    gene_expr.columns = samples["ontology_structure_acronym"].values[: gene_expr.shape[1]]

    return gene_expr


# ── multi-donor aggregation ────────────────────────────────────────────────

def get_microarray_expression_all_donors(
    gene_symbols: list[str],
    collapse_probes: str = "mean",
) -> pd.DataFrame:
    """Load and concatenate microarray expression across all 6 donors.

    Returns a DataFrame with genes as rows, and columns are
    '{structure_acronym}_{donor_id}' identifiers.
    """
    frames = []
    for donor_id in MICROARRAY_DONORS:
        df = get_microarray_gene_expression(gene_symbols, donor_id, collapse_probes)
        if df.empty:
            continue
        df.columns = [f"{c}_d{donor_id}" for c in df.columns]
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, axis=1)


def get_microarray_region_means(
    gene_symbols: list[str],
    collapse_probes: str = "mean",
    level: str = "structure_acronym",
) -> pd.DataFrame:
    """Average microarray expression across donors per brain region.

    Parameters
    ----------
    level : str
        Column from SampleAnnot to group by. Common choices:
        'structure_acronym', 'structure_name', or use the ontology for
        coarser groupings.

    Returns
    -------
    DataFrame (genes x regions) with mean expression across donors and
    samples within each region.
    """
    all_data = []
    for donor_id in MICROARRAY_DONORS:
        probes = load_microarray_probes(donor_id)
        samples = load_microarray_samples(donor_id)
        expr = load_microarray_expression(donor_id)

        mask = probes["gene_symbol"].isin(gene_symbols)
        probe_subset = probes.loc[mask, ["probe_id", "gene_symbol"]]
        if probe_subset.empty:
            continue

        expr_subset = expr.loc[expr.index.isin(probe_subset["probe_id"])]
        expr_subset = expr_subset.merge(
            probe_subset.set_index("probe_id")[["gene_symbol"]],
            left_index=True,
            right_index=True,
        )

        agg_func = {"mean": "mean", "median": "median", "max": "max"}[collapse_probes]
        gene_expr = expr_subset.groupby("gene_symbol").agg(agg_func)
        gene_expr.columns = samples[level].values[: gene_expr.shape[1]]

        # Melt to long format for averaging across regions
        long = gene_expr.T.reset_index()
        long = long.rename(columns={"index": "region"})
        long["donor_id"] = donor_id
        all_data.append(long)

    if not all_data:
        return pd.DataFrame()

    combined = pd.concat(all_data, ignore_index=True)
    region_means = combined.groupby("region")[gene_symbols].mean()

    return region_means.T  # genes x regions
