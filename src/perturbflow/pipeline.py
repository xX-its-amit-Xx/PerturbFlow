"""End-to-end pipeline orchestration.

Composes :mod:`perturbflow.io`, :mod:`.guide_assignment`,
:mod:`.perturbation_analysis`, :mod:`.qc`, :mod:`.de`, :mod:`.downstream`,
and :mod:`.report` into a single ``run()`` function. The CLI and the
Snakemake workflow both call into here so that there's one canonical
ordering of steps and one place to maintain provenance.
"""

from __future__ import annotations

import logging
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from perturbflow import (
    assign_guides,
    compute_cell_state_effects,
    per_cell_qc,
    per_guide_qc,
    per_perturbation_qc,
    read_10x_h5,
    read_10x_mtx,
    read_guide_calls,
    read_guide_metadata,
    read_h5ad,
    run_mixscape,
    run_pseudobulk_de,
    score_pathways,
    write_html_report,
)
from perturbflow._version import __version__
from perturbflow.config import PerturbFlowConfig
from perturbflow.perturbation_analysis import compute_perturbation_signature
from perturbflow.provenance import collect as collect_provenance
from perturbflow.provenance import write as write_provenance
from perturbflow.qc import detect_doublets

logger = logging.getLogger(__name__)


@dataclass
class PipelineArtifacts:
    """Handles returned by :func:`run` so callers can post-process programmatically."""

    adata: ad.AnnData
    per_cell: pd.DataFrame
    per_guide: pd.DataFrame
    per_perturbation: pd.DataFrame
    de_results: dict[str, pd.DataFrame]
    pathway_scores: pd.DataFrame
    cell_state: pd.DataFrame
    report_path: Path | None
    provenance: dict[str, Any]


def run(
    config: PerturbFlowConfig,
    *,
    config_path: str | Path | None = None,
    outdir: str | Path | None = None,
) -> PipelineArtifacts:
    """Execute the full pipeline.

    Parameters
    ----------
    config
        Validated :class:`PerturbFlowConfig`.
    config_path
        Optional path the config was loaded from, for provenance.
    outdir
        Override for ``config.run.outdir``. Useful when running on AWS Batch
        with a scratch path different from the YAML default.
    """
    _set_seeds(config.run.seed)

    out_root = Path(outdir or config.run.outdir)
    out_root.mkdir(parents=True, exist_ok=True)
    (out_root / "qc").mkdir(exist_ok=True)
    (out_root / "de").mkdir(exist_ok=True)
    (out_root / "figures").mkdir(exist_ok=True)

    provenance = collect_provenance(config, config_path=config_path)
    write_provenance(provenance, out_root / "provenance.json")

    # --- 1. Load inputs --------------------------------------------------
    kind, path = config.input.matrix_source()
    if kind == "h5":
        adata = read_10x_h5(path)
    elif kind == "mtx":
        adata = read_10x_mtx(path)
    else:
        adata = read_h5ad(path)

    guide_calls = read_guide_calls(config.input.guide_calls)
    guide_metadata = read_guide_metadata(config.input.guide_metadata)

    # --- 2a. Optional GEX-side doublet detection ------------------------
    if config.qc.detect_doublets:
        adata = detect_doublets(
            adata,
            expected_doublet_rate=config.qc.expected_doublet_rate,
            random_state=config.run.seed,
        )
        if config.qc.drop_doublets:
            keep = ~adata.obs["predicted_doublet"].astype(bool)
            n_drop = int((~keep).sum())
            if n_drop:
                logger.info("Dropping %d predicted doublets before guide assignment", n_drop)
                adata = adata[keep.values].copy()

    # --- 2b. Guide assignment -------------------------------------------
    adata = assign_guides(adata, guide_calls, guide_metadata, config=config.guide_assignment)

    # --- 3. QC tables (pre-mixscape; per_cell rerun later for final report)
    per_cell = per_cell_qc(adata, config=config.qc)
    per_cell.to_csv(out_root / "qc" / "per_cell.csv", index=False)
    per_guide = per_guide_qc(adata, guide_metadata)
    per_guide.to_csv(out_root / "qc" / "per_guide.csv", index=False)

    # --- 4. Mixscape -----------------------------------------------------
    if config.perturbation_analysis.enable_mixscape:
        compute_perturbation_signature(adata, config=config.perturbation_analysis)
        adata = run_mixscape(adata, config=config.perturbation_analysis)

    # Now compute per-perturbation QC (depends on Mixscape).
    per_pert = per_perturbation_qc(adata, control_label=config.perturbation_analysis.control_label)
    per_pert.to_csv(out_root / "qc" / "per_perturbation.csv", index=False)

    # --- 5. Embedding (needed for UMAP overlays in the report) -----------
    _ensure_embedding(adata, seed=config.run.seed)
    cell_state = compute_cell_state_effects(
        adata, control_label=config.perturbation_analysis.control_label
    )
    cell_state.to_csv(out_root / "qc" / "cell_state_effects.csv", index=False)

    # --- 6. DE -----------------------------------------------------------
    de_results: dict[str, pd.DataFrame] = {}
    if config.de.enable:
        de_results = run_pseudobulk_de(
            adata,
            config=config.de,
            input_config=config.input,
            control_label=config.perturbation_analysis.control_label,
        )
        for pert, df in de_results.items():
            safe = pert.replace("/", "_").replace(" ", "_")
            df.to_csv(out_root / "de" / f"{safe}.csv", index=False)
            df.to_parquet(out_root / "de" / f"{safe}.parquet", index=False)

    # --- 7. Pathway scoring ---------------------------------------------
    pathway_scores = score_pathways(de_results, config=config.downstream)
    if not pathway_scores.empty:
        pathway_scores.to_csv(out_root / "de" / "pathway_scores.csv", index=False)

    # --- 8. Final per-cell QC (now with Mixscape columns) + report ------
    per_cell_final = per_cell_qc(adata, config=config.qc)
    per_cell_final.to_csv(out_root / "qc" / "per_cell.csv", index=False)

    report_path: Path | None = None
    if config.report.enable and config.report.bundle_html:
        umap = adata.obsm.get("X_umap")
        labels = adata.obs["perturbation"].astype(str) if "perturbation" in adata.obs else None
        report_path = write_html_report(
            out_path=out_root / "report.html",
            run_name=config.run.name,
            seed=config.run.seed,
            version=__version__,
            git_rev=str(provenance["git"].get("revision", "unknown"))[:12],
            per_cell=per_cell_final,
            per_perturbation=per_pert,
            de_results=de_results,
            pathway_scores=pathway_scores,
            cell_state=cell_state,
            umap_coords=np.asarray(umap) if umap is not None else None,
            umap_labels=labels,
            provenance=provenance,
        )

    # Persist final AnnData for downstream notebook exploration.
    adata.write_h5ad(out_root / "perturbflow.h5ad")

    return PipelineArtifacts(
        adata=adata,
        per_cell=per_cell_final,
        per_guide=per_guide,
        per_perturbation=per_pert,
        de_results=de_results,
        pathway_scores=pathway_scores,
        cell_state=cell_state,
        report_path=report_path,
        provenance=provenance,
    )


def _set_seeds(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    # pertpy's neighbor graph uses torch in some paths; seed it if present.
    try:
        import torch

        torch.manual_seed(seed)
    except ImportError:
        pass


def _ensure_embedding(adata: ad.AnnData, *, seed: int) -> None:
    """Make sure ``adata.obsm['X_umap']`` exists.

    The downstream cell-state effects plot needs a UMAP. We don't overwrite
    one the user already computed (they may have applied batch correction).
    """
    if "X_umap" in adata.obsm:
        return
    if "X_pca" not in adata.obsm:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=min(2000, adata.n_vars - 1))
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=min(30, min(adata.shape) - 1), random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=15, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
