"""QC tables for Perturb-seq runs.

Three levels:

- **per_cell** — the usual single-cell QC stats (n_genes, total_counts,
  pct_mito) plus the Perturb-seq-specific (guide, guide_umi, dominance
  ratio, mixscape_class).
- **per_guide** — coverage of every guide in the library: how many cells
  it was called in, what the mean UMI count is, what fraction of carrying
  cells were called as "escaped" (NP) by Mixscape. A guide called in zero
  cells almost always means the guide library FASTQ wasn't actually
  sequenced; we surface this loudly.
- **per_perturbation** — collapsed to the perturbation level: cells, NP
  fraction, top-DE-gene knockdown magnitude (so you can sanity-check
  on-target activity per perturbation independent of pseudobulk DE).
"""

from __future__ import annotations

import logging

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from perturbflow.config import QCConfig

logger = logging.getLogger(__name__)


def per_cell_qc(
    adata: ad.AnnData,
    *,
    config: QCConfig | None = None,
) -> pd.DataFrame:
    """Build a per-cell QC table.

    Computes total_counts, n_genes, pct_mito on the fly and joins them
    with whatever Perturb-seq-specific columns are present (guide,
    guide_umi, second_guide_umi, perturbation, assignment_status,
    mixscape_class).
    """
    cfg = config or QCConfig()
    X = adata.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)

    total_counts = np.asarray(X.sum(axis=1)).ravel()
    n_genes = np.asarray((X > 0).sum(axis=1)).ravel()

    mito_mask = adata.var_names.str.startswith(cfg.mito_prefix)
    if mito_mask.any():
        mito_counts = np.asarray(X[:, mito_mask.values].sum(axis=1)).ravel()
        with np.errstate(divide="ignore", invalid="ignore"):
            pct_mito = np.where(total_counts > 0, 100 * mito_counts / total_counts, 0.0)
    else:
        pct_mito = np.zeros(adata.n_obs)
        logger.warning(
            "No genes matched mito_prefix=%r; pct_mito set to 0 across all cells",
            cfg.mito_prefix,
        )

    df = pd.DataFrame(
        {
            "cell_barcode": adata.obs_names,
            "total_counts": total_counts.astype(int),
            "n_genes": n_genes.astype(int),
            "pct_mito": pct_mito.astype(float),
        }
    )
    for col in (
        "guide",
        "guide_umi",
        "second_guide_umi",
        "perturbation",
        "is_control",
        "assignment_status",
        "mixscape_class",
        "mixscape_class_global",
        "mixscape_perturbed",
    ):
        if col in adata.obs.columns:
            df[col] = adata.obs[col].values

    df["passes_qc"] = (df["n_genes"] >= cfg.min_genes_per_cell) & (
        df["pct_mito"] <= cfg.max_pct_mito
    )
    return df


def per_guide_qc(
    adata: ad.AnnData,
    guide_metadata: pd.DataFrame,
    *,
    guide_key: str = "guide",
) -> pd.DataFrame:
    """Per-guide coverage table.

    Every guide from the metadata appears in the table even if it was
    called in zero cells, so you can spot dropouts. ``escape_fraction`` is
    NaN when Mixscape hasn't been run.
    """
    counts = adata.obs.loc[adata.obs[guide_key].notna(), guide_key].astype(str).value_counts()
    have_mixscape = (
        "mixscape_class_global" in adata.obs.columns or "mixscape_perturbed" in adata.obs.columns
    )

    rows = []
    for _, row in guide_metadata.iterrows():
        gid = row["guide_id"]
        n_cells = int(counts.get(gid, 0))
        mean_umi = (
            float(adata.obs.loc[adata.obs[guide_key] == gid, "guide_umi"].mean())
            if "guide_umi" in adata.obs.columns and n_cells > 0
            else 0.0
        )
        escape = float("nan")
        if have_mixscape and n_cells > 0:
            sub = adata.obs[adata.obs[guide_key] == gid]
            if "mixscape_class_global" in sub.columns:
                escape = float((sub["mixscape_class_global"].astype(str) == "NP").mean())
            elif "mixscape_perturbed" in sub.columns and not bool(row["is_control"]):
                escape = float(1.0 - sub["mixscape_perturbed"].astype(bool).mean())
        rows.append(
            {
                "guide_id": gid,
                "target_gene": row["target_gene"],
                "is_control": bool(row["is_control"]),
                "n_cells": n_cells,
                "mean_guide_umi": mean_umi,
                "escape_fraction": escape,
            }
        )
    return (
        pd.DataFrame(rows)
        .sort_values(["is_control", "target_gene", "guide_id"])
        .reset_index(drop=True)
    )


def per_perturbation_qc(
    adata: ad.AnnData,
    *,
    perturbation_key: str = "perturbation",
    control_label: str = "NT",
) -> pd.DataFrame:
    """Per-perturbation summary: cells, escape rate, on-target knockdown.

    The on-target log2FC is the canonical "did the screen work?" check.
    When Mixscape has run we compute it on the KO arm only — escaped (NP)
    cells dilute the signal toward the null and would make a real
    knockdown look like noise. The ``on_target_log2fc_all_cells`` column
    is provided alongside for transparency: a big gap between the two
    columns is a sign that escape is the dominant story.
    """
    if perturbation_key not in adata.obs.columns:
        raise KeyError(f"adata.obs[{perturbation_key!r}] not found")

    perts = adata.obs[perturbation_key].astype(str)
    ctrl_mask = (perts == control_label).values
    if not ctrl_mask.any():
        raise ValueError(f"No control cells with perturbation == {control_label!r}")

    # Build a log-normalized view on the fly without mutating caller's adata.
    X = adata.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    lib = np.asarray(X.sum(axis=1)).ravel()
    lib[lib == 0] = 1
    norm = X.multiply(1e4 / lib[:, None]).log1p()  # log(1+CPM_10k)
    norm = norm.toarray() if sparse.issparse(norm) else np.asarray(norm)

    var_index = pd.Index(adata.var_names)
    ctrl_mean = norm[ctrl_mask].mean(axis=0)
    have_mixscape = (
        "mixscape_class_global" in adata.obs.columns or "mixscape_perturbed" in adata.obs.columns
    )

    rows = []
    for pert in sorted(set(perts)):
        if pert == control_label:
            continue
        mask = (perts == pert).values
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        escape = float("nan")
        ko_mask = mask.copy()
        if "mixscape_class_global" in adata.obs.columns:
            sub = adata.obs.loc[mask, "mixscape_class_global"].astype(str)
            escape = float((sub == "NP").mean())
            ko_mask = mask & (adata.obs["mixscape_class_global"].astype(str) == "KO").values
        elif "mixscape_perturbed" in adata.obs.columns:
            escape = float(1.0 - adata.obs.loc[mask, "mixscape_perturbed"].astype(bool).mean())
            ko_mask = mask & adata.obs["mixscape_perturbed"].astype(bool).values

        n_ko_cells = int(ko_mask.sum())
        on_target_lfc = float("nan")
        on_target_lfc_all = float("nan")
        if pert in var_index:
            g = var_index.get_loc(pert)
            ref_mean = float(ctrl_mean[g])
            on_target_lfc_all = float(norm[mask, g].mean()) - ref_mean
            if have_mixscape and n_ko_cells > 0:
                # KO-only on-target: this is the right number to report as
                # "did the perturbation work in cells that perturbed?"
                on_target_lfc = float(norm[ko_mask, g].mean()) - ref_mean
            else:
                on_target_lfc = on_target_lfc_all

        rows.append(
            {
                "perturbation": pert,
                "n_cells": n_cells,
                "n_ko_cells": n_ko_cells if have_mixscape else n_cells,
                "escape_fraction": escape,
                "on_target_log2fc": on_target_lfc,
                "on_target_log2fc_all_cells": on_target_lfc_all,
            }
        )
    return pd.DataFrame(rows).sort_values("perturbation").reset_index(drop=True)
