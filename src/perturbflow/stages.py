"""Stage runners used by the Snakemake DAG.

Each stage takes one or more h5ad / CSV inputs and writes a single
artifact. The Snakefile composes them into a DAG so a re-run can resume
from the last successful stage, and a Snakemake report renders the per-
stage timing.

The pipeline.run() entrypoint and the CLI ``perturbflow run`` still call
the high-level orchestration; this module is what the per-stage rules
shell out to.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

from perturbflow.config import PerturbFlowConfig
from perturbflow.de import run_pseudobulk_de
from perturbflow.downstream import compute_cell_state_effects, score_pathways
from perturbflow.guide_assignment import assign_guides
from perturbflow.io import (
    read_10x_h5,
    read_10x_mtx,
    read_guide_calls,
    read_guide_metadata,
    read_h5ad,
)
from perturbflow.perturbation_analysis import (
    compute_perturbation_signature,
    run_mixscape,
)
from perturbflow.qc import per_cell_qc, per_guide_qc, per_perturbation_qc

logger = logging.getLogger(__name__)


def stage_load(cfg: PerturbFlowConfig, out_h5ad: Path) -> None:
    """Stage 1: load the input matrix into an AnnData, write to disk."""
    kind, path = cfg.input.matrix_source()
    if kind == "h5":
        adata = read_10x_h5(path)
    elif kind == "mtx":
        adata = read_10x_mtx(path)
    else:
        adata = read_h5ad(path)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad)
    logger.info("stage_load: wrote %s (%d cells, %d genes)", out_h5ad, adata.n_obs, adata.n_vars)


def stage_assign(
    cfg: PerturbFlowConfig,
    in_h5ad: Path,
    out_h5ad: Path,
    per_guide_csv: Path,
) -> None:
    """Stage 2: guide assignment + per-guide QC table."""
    adata = read_h5ad(in_h5ad)
    guide_calls = read_guide_calls(cfg.input.guide_calls)
    guide_metadata = read_guide_metadata(cfg.input.guide_metadata)
    adata = assign_guides(adata, guide_calls, guide_metadata, config=cfg.guide_assignment)
    per_guide_csv.parent.mkdir(parents=True, exist_ok=True)
    per_guide_qc(adata, guide_metadata).to_csv(per_guide_csv, index=False)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad)


def stage_mixscape(cfg: PerturbFlowConfig, in_h5ad: Path, out_h5ad: Path) -> None:
    """Stage 3: Mixscape KO vs NP classification."""
    adata = read_h5ad(in_h5ad)
    if cfg.perturbation_analysis.enable_mixscape:
        compute_perturbation_signature(adata, config=cfg.perturbation_analysis)
        adata = run_mixscape(adata, config=cfg.perturbation_analysis)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad)


def stage_embedding(
    cfg: PerturbFlowConfig, in_h5ad: Path, out_h5ad: Path, cell_state_csv: Path
) -> None:
    """Stage 4: compute UMAP + cell-state effect table."""
    adata = read_h5ad(in_h5ad)
    _ensure_embedding(adata, seed=cfg.run.seed)
    cell_state = compute_cell_state_effects(
        adata, control_label=cfg.perturbation_analysis.control_label
    )
    cell_state_csv.parent.mkdir(parents=True, exist_ok=True)
    cell_state.to_csv(cell_state_csv, index=False)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad)


def stage_qc(cfg: PerturbFlowConfig, in_h5ad: Path, per_cell_csv: Path, per_pert_csv: Path) -> None:
    """Stage 5: per-cell + per-perturbation QC (depends on Mixscape having run)."""
    adata = read_h5ad(in_h5ad)
    per_cell_csv.parent.mkdir(parents=True, exist_ok=True)
    per_cell_qc(adata, config=cfg.qc).to_csv(per_cell_csv, index=False)
    per_perturbation_qc(adata, control_label=cfg.perturbation_analysis.control_label).to_csv(
        per_pert_csv, index=False
    )


def stage_de(cfg: PerturbFlowConfig, in_h5ad: Path, out_dir: Path) -> None:
    """Stage 6: pseudobulk DE per perturbation. Writes one CSV+parquet per pert."""
    if not cfg.de.enable:
        logger.info("stage_de: DE disabled in config; writing empty marker")
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / ".empty").write_text("DE disabled in config\n")
        return
    adata = read_h5ad(in_h5ad)
    results = run_pseudobulk_de(
        adata,
        config=cfg.de,
        input_config=cfg.input,
        control_label=cfg.perturbation_analysis.control_label,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    for pert, df in results.items():
        safe = pert.replace("/", "_").replace(" ", "_")
        df.to_csv(out_dir / f"{safe}.csv", index=False)
        df.to_parquet(out_dir / f"{safe}.parquet", index=False)
    logger.info("stage_de: wrote %d perturbation DE tables to %s", len(results), out_dir)


def stage_pathways(cfg: PerturbFlowConfig, de_dir: Path, out_csv: Path) -> None:
    """Stage 7: pathway scoring from the per-perturbation DE tables."""
    if not cfg.downstream.enable_pathway_scoring:
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=["perturbation", "pathway", "score", "pvalue"]).to_csv(
            out_csv, index=False
        )
        return
    de_results: dict[str, pd.DataFrame] = {}
    for csv in sorted(de_dir.glob("*.csv")):
        de_results[csv.stem] = pd.read_csv(csv)
    pathway_scores = score_pathways(de_results, config=cfg.downstream)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pathway_scores.to_csv(out_csv, index=False)


def _ensure_embedding(adata: ad.AnnData, *, seed: int) -> None:
    """Same embedding logic as :func:`perturbflow.pipeline._ensure_embedding`.

    Kept here as a module-private helper so the Snakemake rules can run
    independently of the orchestrator.
    """
    if "X_umap" in adata.obsm:
        return
    if "X_pca" not in adata.obsm:
        tmp = adata.copy()
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
        sc.pp.scale(tmp, max_value=10)
        sc.tl.pca(tmp, n_comps=min(30, min(tmp.shape) - 1), random_state=seed)
        adata.obsm["X_pca"] = tmp.obsm["X_pca"]
    sc.pp.neighbors(adata, n_neighbors=15, random_state=seed)
    sc.tl.umap(adata, random_state=seed)


STAGE_REGISTRY = {
    "load": stage_load,
    "assign": stage_assign,
    "mixscape": stage_mixscape,
    "embedding": stage_embedding,
    "qc": stage_qc,
    "de": stage_de,
    "pathways": stage_pathways,
}
"""Mapping of stage name -> runner. Used by the CLI's ``stage`` subcommand."""


__all__ = [
    "STAGE_REGISTRY",
    "stage_assign",
    "stage_de",
    "stage_embedding",
    "stage_load",
    "stage_mixscape",
    "stage_pathways",
    "stage_qc",
]
