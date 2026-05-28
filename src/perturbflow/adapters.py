"""Adapters: turn upstream tool output into PerturbFlow's canonical schema.

Most CRISPR-screen pipelines either go through CellRanger (10x screens, by
far the most common) or land at one of a few public-release schemas
(Replogle 2022, Dixit 2016). PerturbFlow's documented input is a
deliberately minimal long-format CSV, but no real lab actually produces
that file by hand — so we ship adapters that take what their upstream
tooling already produced.

If you're starting a new analysis from CellRanger output, look at
:func:`read_cellranger_protospacer_calls` first.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def read_cellranger_protospacer_calls(
    path: str | Path,
    *,
    drop_multifeature: bool = False,
) -> pd.DataFrame:
    """Read the CellRanger 7.x ``protospacer_calls_per_cell.csv``.

    CellRanger's CRISPR Guide Capture pipeline emits one row per cell with
    pipe-separated guide assignments::

        cell_barcode,num_features,feature_call,num_umis
        AAACCTGAGAAACGAG-1,1,Geneset_A_g1,42
        AAACCTGAGTTAACGA-1,2,Geneset_A_g1|Geneset_A_g2,23|17

    We pivot this into our long format::

        cell_barcode, guide_id, umi_count

    Multi-feature cells are emitted as one row per (cell, guide). The
    downstream :func:`perturbflow.guide_assignment.assign_guides` flags
    them as ``multi-guide`` based on ``max_guides`` — this function does
    NOT pre-filter them out (set ``drop_multifeature=True`` to do that
    here instead, but flagging downstream is the recommended path because
    it preserves the QC count).

    Parameters
    ----------
    path
        Path to ``crispr_analysis/protospacer_calls_per_cell.csv``.
    drop_multifeature
        If True, only emit single-feature cells. Default False.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"CellRanger protospacer file not found: {p}")
    df = pd.read_csv(p)
    expected = {"cell_barcode", "num_features", "feature_call", "num_umis"}
    missing = expected - set(df.columns)
    if missing:
        raise ValueError(
            f"{p}: not a CellRanger protospacer_calls_per_cell.csv "
            f"(missing columns {sorted(missing)}; got {list(df.columns)})"
        )

    if drop_multifeature:
        df = df[df["num_features"] == 1].copy()

    df["feature_call"] = df["feature_call"].astype(str)
    df["num_umis"] = df["num_umis"].astype(str)
    rows: list[dict[str, object]] = []
    for _, row in df.iterrows():
        guides = row["feature_call"].split("|")
        umis = row["num_umis"].split("|")
        if len(guides) != len(umis):
            raise ValueError(
                f"Inconsistent feature_call / num_umis split for cell "
                f"{row['cell_barcode']!r}: {len(guides)} guides vs {len(umis)} UMI counts"
            )
        for g, u in zip(guides, umis, strict=True):
            try:
                umi_int = int(u)
            except ValueError as e:
                raise ValueError(
                    f"Non-integer UMI count {u!r} for cell {row['cell_barcode']!r}"
                ) from e
            rows.append({"cell_barcode": row["cell_barcode"], "guide_id": g, "umi_count": umi_int})

    out = pd.DataFrame(rows, columns=["cell_barcode", "guide_id", "umi_count"])
    logger.info(
        "CellRanger protospacer: %d cells, %d unique guides, %d (cell, guide) rows from %s",
        out["cell_barcode"].nunique(),
        out["guide_id"].nunique(),
        len(out),
        p,
    )
    return out


def guide_metadata_from_cellranger_features(
    features_path: str | Path,
    *,
    nontargeting_pattern: str = "non-targeting",
) -> pd.DataFrame:
    """Build a guide-metadata table from CellRanger's ``features.tsv``.

    Parameters
    ----------
    features_path
        Path to the CellRanger ``features.tsv`` (or ``features.tsv.gz``).
        Expected columns are tab-separated: ``feature_id, feature_name,
        feature_type, [target_gene_id, target_gene_name]``. Only rows
        with ``feature_type == 'CRISPR Guide Capture'`` are kept.
    nontargeting_pattern
        Substring (case-insensitive) on ``target_gene_name`` that marks a
        non-targeting / scrambled control guide. Adjust if your library
        uses ``"NT"``, ``"scrambled"``, or a custom prefix.
    """
    p = Path(features_path)
    if not p.exists():
        raise FileNotFoundError(f"CellRanger features file not found: {p}")
    df = pd.read_csv(p, sep="\t", header=None)
    if df.shape[1] < 5:
        raise ValueError(
            f"{p}: CRISPR libraries need features.tsv with at least 5 columns "
            "(feature_id, feature_name, feature_type, target_gene_id, "
            f"target_gene_name); got {df.shape[1]}"
        )
    df.columns = [
        "feature_id",
        "feature_name",
        "feature_type",
        "target_gene_id",
        "target_gene_name",
    ] + [f"extra_{i}" for i in range(df.shape[1] - 5)]
    df = df[df["feature_type"] == "CRISPR Guide Capture"].copy()
    if df.empty:
        raise ValueError(f"{p}: no rows with feature_type == 'CRISPR Guide Capture'")
    is_control = (
        df["target_gene_name"].astype(str).str.contains(nontargeting_pattern, case=False, na=False)
    )
    out = pd.DataFrame(
        {
            "guide_id": df["feature_name"].astype(str),
            "target_gene": np.where(is_control, "", df["target_gene_name"].astype(str)),
            "is_control": is_control.values,
        }
    )
    if out["guide_id"].duplicated().any():
        dups = out.loc[out["guide_id"].duplicated(), "guide_id"].unique().tolist()
        raise ValueError(f"{p}: duplicate guide feature_name values: {dups}")
    return out.reset_index(drop=True)


def from_replogle_2022_anndata(
    adata: ad.AnnData,
    *,
    perturbation_obs_col: str = "gene_target",
    nontargeting_label: str = "non-targeting",
) -> ad.AnnData:
    """Adapt a Replogle 2022 figshare AnnData to PerturbFlow's expected layout.

    Replogle's public release uses ``adata.obs[perturbation_obs_col]`` to
    record the perturbed gene (or ``"non-targeting"`` for controls) and
    keeps cells already filtered. We translate this to the conventions
    PerturbFlow's downstream modules expect:

    - ``obs['perturbation']`` — gene name or ``"NT"`` for controls
    - ``obs['guide']`` — synthesized as ``"<gene>_g1"`` since the public
      release doesn't preserve individual guide identity per cell
    - ``obs['assignment_status']`` — ``"assigned"`` for all cells
    - ``obs['is_control']`` — boolean

    After this, you can run perturbation_analysis / de directly. You do
    NOT need to run ``assign_guides`` because the upstream release already
    did it for you. The guide_metadata table is not needed either.
    """
    if perturbation_obs_col not in adata.obs.columns:
        raise KeyError(
            f"adata.obs[{perturbation_obs_col!r}] not found; "
            f"available: {list(adata.obs.columns)[:8]}..."
        )
    out = adata.copy()
    gene = out.obs[perturbation_obs_col].astype(str)
    is_control = gene.str.lower() == nontargeting_label.lower()
    out.obs["perturbation"] = pd.Categorical(np.where(is_control, "NT", gene))
    out.obs["guide"] = pd.Categorical(np.where(is_control, "NT_g1", gene + "_g1"))
    out.obs["assignment_status"] = pd.Categorical(
        ["assigned"] * out.n_obs,
        categories=["assigned", "ambiguous", "multi-guide", "unassigned"],
    )
    out.obs["is_control"] = is_control.values
    out.uns.setdefault("perturbflow", {})
    out.uns["perturbflow"]["adapter"] = {
        "source": "replogle_2022",
        "perturbation_obs_col": perturbation_obs_col,
        "n_cells": int(out.n_obs),
        "n_perturbations": int(out.obs["perturbation"].nunique()),
    }
    logger.info(
        "Replogle 2022 adapter: %d cells, %d perturbations (control = NT, %d cells)",
        out.n_obs,
        out.obs["perturbation"].nunique(),
        int(is_control.sum()),
    )
    return out
