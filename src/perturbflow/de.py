"""Pseudobulk differential expression for Perturb-seq.

We run DE on pseudobulk count matrices, not on single cells. Single-cell DE
inflates p-values because every cell is treated as an independent replicate
even when the biological unit (the perturbed clone or the donor sample) is
much smaller. Pseudobulking — summing counts over groups of cells that
share a biological replicate — gives DESeq2/edgeR-style statistics the kind
of replication they were designed for. See `Squair et al. 2021
<https://www.nature.com/articles/s41467-021-25960-2>`_ for the empirical
case against single-cell DE.

When a sample column is configured (e.g. donor or biological replicate),
each (perturbation, sample) pair is one pseudobulk row. When no sample
column exists — common for single-donor screens — we synthesize pseudo-
replicates by deterministic hashing of cell barcodes. Pseudo-replicates do
not buy you statistical power that wasn't there to begin with; they just
give DESeq2 a sane dispersion estimate. Real biological replicates are
strictly better and you should add a sample column when you have one.
"""

from __future__ import annotations

import hashlib
import logging
from importlib.util import find_spec

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

from perturbflow.config import DEConfig, InputConfig
from perturbflow.validation import assert_raw_counts, warn_if_pseudoreplicates

logger = logging.getLogger(__name__)


def make_pseudo_replicates(
    adata: ad.AnnData,
    *,
    n_replicates: int,
    perturbation_key: str = "perturbation",
    out_col: str = "pseudo_replicate",
    seed: int = 0,
) -> ad.AnnData:
    """Add a deterministic pseudo-replicate column to ``adata.obs``.

    Cells within the same perturbation are split into ``n_replicates`` bins
    by hashing the cell barcode. Hashing makes the split deterministic
    across re-runs and immune to row order changes, which matters when the
    pipeline is restarted and pseudobulk counts must reproduce.
    """
    if perturbation_key not in adata.obs.columns:
        raise KeyError(f"adata.obs[{perturbation_key!r}] not found")
    if n_replicates < 1:
        raise ValueError(f"n_replicates must be >= 1, got {n_replicates}")
    salt = str(seed).encode()
    reps = np.empty(adata.n_obs, dtype=np.int32)
    for i, bc in enumerate(adata.obs_names):
        h = hashlib.sha1(salt + str(bc).encode()).digest()
        reps[i] = int.from_bytes(h[:4], "little") % n_replicates
    adata.obs[out_col] = pd.Series(reps, index=adata.obs_names).astype("category")
    return adata


def make_pseudobulk(
    adata: ad.AnnData,
    *,
    perturbation_key: str = "perturbation",
    sample_key: str | None = None,
    layer: str | None = None,
    min_cells: int = 1,
) -> ad.AnnData:
    """Sum raw counts within each (perturbation, sample) group.

    Returns a new AnnData of shape ``(n_groups, n_genes)`` whose ``obs``
    table records ``perturbation``, ``sample``, and ``n_cells``. ``X`` is
    integer counts (DESeq2 needs raw counts, not log-normalized).
    """
    counts = adata.X if layer is None else adata.layers[layer]
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)

    if sample_key is None:
        raise ValueError(
            "make_pseudobulk requires a sample_key (use make_pseudo_replicates to "
            "synthesize one if you have no biological replicates)"
        )

    perts = adata.obs[perturbation_key].astype(str).values
    samples = adata.obs[sample_key].astype(str).values
    pair = pd.DataFrame({"pert": perts, "sample": samples})
    pair["row"] = np.arange(len(pair))

    groups = pair.groupby(["pert", "sample"], sort=True)
    group_idx = list(groups.groups.keys())

    pb_rows: list[np.ndarray] = []
    n_cells_per_row: list[int] = []
    for key in group_idx:
        rows = groups.get_group(key)["row"].values
        if len(rows) < min_cells:
            continue
        summed = np.asarray(counts[rows].sum(axis=0)).ravel()
        pb_rows.append(summed)
        n_cells_per_row.append(len(rows))

    kept_keys = [k for k, n in zip(group_idx, n_cells_per_row, strict=True)]
    if not pb_rows:
        raise ValueError(
            f"No pseudobulk groups passed min_cells={min_cells}. Check guide assignment yield."
        )

    pb_matrix = np.vstack(pb_rows).astype(np.int64)
    pb_obs = pd.DataFrame(
        {
            "perturbation": [k[0] for k in kept_keys],
            "sample": [k[1] for k in kept_keys],
            "n_cells": n_cells_per_row,
        },
        index=[f"{p}__{s}" for p, s in kept_keys],
    )
    pb = ad.AnnData(X=pb_matrix, obs=pb_obs, var=adata.var.copy())
    return pb


def run_pseudobulk_de(
    adata: ad.AnnData,
    *,
    config: DEConfig | None = None,
    input_config: InputConfig | None = None,
    perturbation_key: str = "perturbation",
    control_label: str = "NT",
    n_threads: int = 1,
) -> dict[str, pd.DataFrame]:
    """Run pseudobulk DE for each non-control perturbation vs the control set.

    Returns a dict ``{perturbation: results_df}``. Each ``results_df`` has
    one row per gene with columns ``log2FoldChange``, ``pvalue``, ``padj``,
    ``baseMean``, ``stat``, ``significant``.

    Uses :mod:`pydeseq2` when available; otherwise falls back to a per-gene
    Welch's t-test on log-CPM with BH multiple-testing correction. The
    fallback is fine for sanity checks but you should install pydeseq2
    for any analysis that's going into a manuscript.
    """
    cfg = config or DEConfig()
    icfg = input_config or InputConfig()

    if perturbation_key not in adata.obs.columns:
        raise KeyError(f"adata.obs[{perturbation_key!r}] not found")

    # Gate 1: DESeq2 needs raw integer counts. Fail loud if X looks
    # log-normalized — silently producing a garbage volcano is the worst
    # possible failure mode.
    assert_raw_counts(adata, strict=True)

    work = adata
    if cfg.use_mixscape_filter and "mixscape_perturbed" in adata.obs.columns:
        # Keep KO cells (true perturbed) and all control cells; drop escaped.
        ctrl_mask = adata.obs[perturbation_key].astype(str) == control_label
        keep = adata.obs["mixscape_perturbed"].astype(bool) | ctrl_mask
        if not keep.any():
            raise ValueError("Mixscape filter eliminated all cells; check thresholds.")
        work = adata[keep.values].copy()
        logger.info(
            "DE: filtered to Mixscape KO + control cells (%d -> %d cells)",
            adata.n_obs,
            work.n_obs,
        )

    sample_key = icfg.sample_col
    if sample_key is None:
        n_perts_with_ko = work.obs[perturbation_key].nunique()
        warn_if_pseudoreplicates(
            sample_col=None,
            n_pseudo_replicates=icfg.n_pseudo_replicates,
            n_perturbations=int(n_perts_with_ko),
        )
        work = make_pseudo_replicates(
            work,
            n_replicates=icfg.n_pseudo_replicates,
            perturbation_key=perturbation_key,
        )
        sample_key = "pseudo_replicate"

    pb = make_pseudobulk(
        work,
        perturbation_key=perturbation_key,
        sample_key=sample_key,
        min_cells=cfg.min_cells_per_replicate,
    )

    pert_levels = sorted(set(pb.obs["perturbation"]) - {control_label})
    n_ctrl_reps = int((pb.obs["perturbation"] == control_label).sum())
    if n_ctrl_reps < cfg.min_replicates_per_group:
        raise ValueError(
            f"Need at least {cfg.min_replicates_per_group} control pseudobulk "
            f"replicates; got {n_ctrl_reps}"
        )

    have_pydeseq2 = find_spec("pydeseq2") is not None
    if not have_pydeseq2:
        logger.warning(
            "pydeseq2 not installed; using Welch's t-test fallback. "
            "Install with `pip install perturbflow[de]` for the published DESeq2 implementation."
        )

    results: dict[str, pd.DataFrame] = {}
    skipped_all_np: list[str] = []
    for pert in pert_levels:
        n_treat = int((pb.obs["perturbation"] == pert).sum())
        if n_treat == 0:
            # All cells of this perturbation were Mixscape-NP and dropped.
            skipped_all_np.append(pert)
            continue
        if n_treat < cfg.min_replicates_per_group:
            logger.warning(
                "Skipping DE for perturbation %r: only %d pseudobulk replicates (<%d)",
                pert,
                n_treat,
                cfg.min_replicates_per_group,
            )
            continue
        sub_mask = pb.obs["perturbation"].isin([pert, control_label])
        sub = pb[sub_mask].copy()
        if have_pydeseq2:
            df = _run_pydeseq2(sub, pert=pert, control=control_label, n_threads=n_threads)
        else:
            df = _run_ttest(sub, pert=pert, control=control_label)
        df["significant"] = (df["padj"] < cfg.padj_threshold) & (
            df["log2FoldChange"].abs() >= cfg.lfc_threshold
        )
        df["perturbation"] = pert
        results[pert] = df
    if skipped_all_np:
        logger.warning(
            "Skipped DE for %d perturbations whose treatment arm was empty after "
            "Mixscape filtering (all cells classified NP): %s",
            len(skipped_all_np),
            sorted(skipped_all_np)[:5],
        )

    logger.info("Pseudobulk DE complete for %d perturbations", len(results))
    return results


def _run_pydeseq2(pb: ad.AnnData, *, pert: str, control: str, n_threads: int) -> pd.DataFrame:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    counts = pd.DataFrame(pb.X, index=pb.obs_names, columns=pb.var_names)
    # DESeq2 needs integer counts; pseudobulk sums already are.
    counts = counts.astype(int)
    metadata = pb.obs[["perturbation"]].copy()
    # DESeq2 prefers explicit, alphabetically-sorted factor levels — name the
    # control "AAA_control" so it gets used as the reference level.
    metadata["condition"] = np.where(metadata["perturbation"] == control, "control", "treatment")
    metadata = metadata.drop(columns=["perturbation"])
    inference = DefaultInference(n_cpus=n_threads)
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design_factors="condition",
        refit_cooks=True,
        inference=inference,
        quiet=True,
    )
    dds.deseq2()
    stats_runner = DeseqStats(
        dds,
        contrast=["condition", "treatment", "control"],
        inference=inference,
        quiet=True,
    )
    stats_runner.summary()
    df = stats_runner.results_df.copy()
    df.index.name = "gene"
    return df.reset_index()


def _run_ttest(pb: ad.AnnData, *, pert: str, control: str) -> pd.DataFrame:
    """CPM-based Welch's t-test fallback (used when pydeseq2 is unavailable)."""
    counts = np.asarray(pb.X, dtype=float)
    lib_size = counts.sum(axis=1, keepdims=True)
    lib_size[lib_size == 0] = 1
    cpm = counts / lib_size * 1e6
    logcpm = np.log2(cpm + 1)

    mask_t = (pb.obs["perturbation"] == pert).values
    mask_c = (pb.obs["perturbation"] == control).values
    treat = logcpm[mask_t]
    ctrl = logcpm[mask_c]

    base_mean = counts.mean(axis=0)
    lfc = treat.mean(axis=0) - ctrl.mean(axis=0)
    # Welch's t per gene, vectorized.
    tt = stats.ttest_ind(treat, ctrl, equal_var=False, axis=0, nan_policy="omit")
    pvals = np.asarray(tt.pvalue)
    pvals = np.where(np.isnan(pvals), 1.0, pvals)
    padj = _bh_fdr(pvals)
    return pd.DataFrame(
        {
            "gene": pb.var_names,
            "baseMean": base_mean,
            "log2FoldChange": lfc,
            "stat": np.asarray(tt.statistic),
            "pvalue": pvals,
            "padj": padj,
        }
    )


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction (no external dep)."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / (np.arange(n) + 1)
    # Enforce monotonicity (the standard BH step-up).
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty_like(adj)
    out[order] = np.clip(adj, 0, 1)
    return out
