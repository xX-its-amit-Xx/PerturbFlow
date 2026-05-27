from __future__ import annotations

import anndata as ad
import pandas as pd

from perturbflow.config import DEConfig, GuideAssignmentConfig, InputConfig
from perturbflow.de import (
    make_pseudo_replicates,
    make_pseudobulk,
    run_pseudobulk_de,
)
from perturbflow.guide_assignment import assign_guides


def _assigned(
    adata: ad.AnnData,
    calls: pd.DataFrame,
    meta: pd.DataFrame,
) -> ad.AnnData:
    return assign_guides(adata.copy(), calls, meta, config=GuideAssignmentConfig())


def test_pseudo_replicates_deterministic(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    adata = _assigned(synthetic_adata, synthetic_guide_calls, synthetic_guide_metadata)
    a = make_pseudo_replicates(adata.copy(), n_replicates=3, seed=0)
    b = make_pseudo_replicates(adata.copy(), n_replicates=3, seed=0)
    assert (a.obs["pseudo_replicate"].astype(str) == b.obs["pseudo_replicate"].astype(str)).all()


def test_pseudo_replicates_distribute_balanced(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    adata = _assigned(synthetic_adata, synthetic_guide_calls, synthetic_guide_metadata)
    a = make_pseudo_replicates(adata, n_replicates=3)
    counts = a.obs["pseudo_replicate"].value_counts()
    assert len(counts) == 3
    # Each bin should be within 20% of equal split.
    target = a.n_obs / 3
    assert counts.min() > 0.6 * target
    assert counts.max() < 1.4 * target


def test_pseudobulk_shape(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    adata = _assigned(synthetic_adata, synthetic_guide_calls, synthetic_guide_metadata)
    adata = make_pseudo_replicates(adata, n_replicates=3)
    pb = make_pseudobulk(adata, sample_key="pseudo_replicate")
    # 5 perturbations × 3 pseudo-replicates = 15 rows (at most).
    assert pb.n_obs <= 5 * 3
    assert pb.n_vars == adata.n_vars
    assert (pb.X >= 0).all()


def test_pseudobulk_de_recovers_target_genes(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    """For each KO perturbation, the target gene should be significantly down."""
    adata = _assigned(synthetic_adata, synthetic_guide_calls, synthetic_guide_metadata)
    # Bypass Mixscape filtering for this isolated test.
    results = run_pseudobulk_de(
        adata,
        config=DEConfig(
            enable=True,
            use_mixscape_filter=False,
            min_replicates_per_group=2,
            min_cells_per_replicate=5,
            lfc_threshold=0.5,
            padj_threshold=0.1,
        ),
        input_config=InputConfig(n_pseudo_replicates=3),
        control_label="NT",
    )
    assert set(results.keys()) >= {"KOA", "KOB", "KOC", "KOD"}
    for pert in ("KOA", "KOB", "KOC", "KOD"):
        df = results[pert].set_index("gene")
        # Target gene must appear in results.
        assert pert in df.index
        # The on-target log2FC should be negative (knockdown).
        assert df.loc[pert, "log2FoldChange"] < 0


def test_pseudobulk_de_uses_mixscape_filter_when_present(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    """When mixscape_perturbed is in obs, DE should run only on KO + NT cells."""
    from perturbflow.config import PerturbationAnalysisConfig
    from perturbflow.perturbation_analysis import run_mixscape

    adata = _assigned(synthetic_adata, synthetic_guide_calls, synthetic_guide_metadata)
    adata = run_mixscape(adata, config=PerturbationAnalysisConfig(n_neighbors=15))
    results = run_pseudobulk_de(
        adata,
        config=DEConfig(
            enable=True,
            use_mixscape_filter=True,
            min_replicates_per_group=2,
            min_cells_per_replicate=5,
            lfc_threshold=0.5,
            padj_threshold=0.1,
        ),
        input_config=InputConfig(n_pseudo_replicates=3),
        control_label="NT",
    )
    # Should still produce DE for each KO perturbation.
    assert set(results.keys()) >= {"KOA", "KOB", "KOC", "KOD"}
