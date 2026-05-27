"""Input readers for PerturbFlow.

Three matrix shapes are supported:

- 10x Genomics ``filtered_feature_bc_matrix.h5`` (recommended; preserves
  feature types so CRISPR Guide Capture features are split out for you).
- 10x Genomics MTX directory (legacy CellRanger output).
- An existing ``.h5ad`` (lets you re-enter the pipeline after a custom
  upstream step — your QC, demultiplexing, batch correction, etc.).

Guide calls and guide metadata are CSV/TSV with a documented schema. We
validate column presence aggressively because silently dropping a malformed
``umi_count`` column wrecks everything downstream.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

logger = logging.getLogger(__name__)


REQUIRED_GUIDE_CALL_COLUMNS = ("cell_barcode", "guide_id", "umi_count")
REQUIRED_GUIDE_METADATA_COLUMNS = ("guide_id", "target_gene", "is_control")


def read_10x_h5(path: str | Path, *, gex_only: bool = True) -> ad.AnnData:
    """Read a CellRanger ``filtered_feature_bc_matrix.h5`` file.

    Parameters
    ----------
    path
        Path to the ``.h5`` file.
    gex_only
        If True (default), keep only Gene Expression features. CRISPR Guide
        Capture features are returned separately by :func:`read_guide_calls`
        — keeping them in ``X`` confuses the QC and DE modules. Set False
        only if you've already verified the matrix is gene-only.

    Returns
    -------
    AnnData
        With ``obs_names`` = cell barcodes and ``var_names`` = gene symbols
        (CellRanger ``features`` "name" column). The Ensembl IDs are kept in
        ``var['gene_ids']`` for unambiguous joins.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"10x H5 file not found: {p}")
    logger.info("Reading 10x H5 matrix: %s", p)
    adata = sc.read_10x_h5(str(p), gex_only=gex_only)
    adata.var_names_make_unique()
    return adata


def read_10x_mtx(directory: str | Path, *, var_names: str = "gene_symbols") -> ad.AnnData:
    """Read a CellRanger MTX directory (``matrix.mtx``, ``barcodes.tsv``, ``features.tsv``)."""
    p = Path(directory)
    if not p.is_dir():
        raise FileNotFoundError(f"10x MTX directory not found: {p}")
    logger.info("Reading 10x MTX directory: %s", p)
    adata = sc.read_10x_mtx(str(p), var_names=var_names, make_unique=True)
    return adata


def read_h5ad(path: str | Path) -> ad.AnnData:
    """Read an AnnData ``.h5ad`` file."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"H5AD file not found: {p}")
    logger.info("Reading h5ad: %s", p)
    return ad.read_h5ad(str(p))


def read_guide_calls(path: str | Path) -> pd.DataFrame:
    """Read a guide-call table.

    Required columns
    ----------------
    cell_barcode : str
        Must match ``adata.obs_names`` after any barcode rewriting (e.g., a
        ``-1`` suffix added by CellRanger).
    guide_id : str
        The guide identifier; will be looked up in the guide metadata table.
    umi_count : int
        Raw UMI count for this (cell, guide) pair. We do *not* normalize
        this — guide assignment uses absolute thresholds because guide
        capture libraries are sequenced shallower than GEX.

    Cells with no guide calls above threshold are added back as
    ``assignment_status == "unassigned"`` by :func:`assign_guides`; you do
    not need to include zero-count rows here.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Guide call file not found: {p}")
    sep = "\t" if p.suffix.lower() in {".tsv", ".txt"} else ","
    df = pd.read_csv(p, sep=sep)
    _validate_columns(df, REQUIRED_GUIDE_CALL_COLUMNS, source=str(p))
    df["umi_count"] = pd.to_numeric(df["umi_count"], errors="raise").astype(int)
    df["cell_barcode"] = df["cell_barcode"].astype(str)
    df["guide_id"] = df["guide_id"].astype(str)
    logger.info(
        "Loaded %d guide calls covering %d cells and %d unique guides from %s",
        len(df),
        df["cell_barcode"].nunique(),
        df["guide_id"].nunique(),
        p,
    )
    return df


def read_guide_metadata(path: str | Path) -> pd.DataFrame:
    """Read a guide metadata table mapping ``guide_id`` to a target perturbation.

    Required columns
    ----------------
    guide_id : str
    target_gene : str
        The gene symbol the guide targets. For non-targeting guides, leave
        blank or set ``is_control = True``; ``target_gene`` is ignored when
        ``is_control`` is true.
    is_control : bool
        Marks non-targeting / scrambled guides. These define the control set
        the rest of the pipeline contrasts against.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Guide metadata file not found: {p}")
    sep = "\t" if p.suffix.lower() in {".tsv", ".txt"} else ","
    df = pd.read_csv(p, sep=sep)
    _validate_columns(df, REQUIRED_GUIDE_METADATA_COLUMNS, source=str(p))
    df["guide_id"] = df["guide_id"].astype(str)
    df["target_gene"] = df["target_gene"].fillna("").astype(str)
    df["is_control"] = df["is_control"].astype(bool)
    if not df["is_control"].any():
        raise ValueError(
            f"{p}: no control guides found (no rows with is_control=True). "
            "A control set is required for perturbation analysis and DE."
        )
    if df["guide_id"].duplicated().any():
        dups = df.loc[df["guide_id"].duplicated(), "guide_id"].unique().tolist()
        raise ValueError(f"{p}: duplicate guide_id values: {dups}")
    return df


def _validate_columns(df: pd.DataFrame, required: tuple[str, ...], *, source: str) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{source}: missing required columns {missing}; got {list(df.columns)}")
