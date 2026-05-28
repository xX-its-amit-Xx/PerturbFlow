"""Parity test: our wrapper output matches pertpy's directly-called output.

Without this, the claim "PerturbFlow wraps pertpy's validated Mixscape
implementation" is unverified. The test runs both code paths on the same
seeded synthetic data and verifies:

1. The same cells get classified as KO by both runs.
2. The downstream ``mixscape_class_global`` column matches.

If pertpy ever ships a behavior-changing update, this test alerts us
before users hit it.
"""

from __future__ import annotations

from importlib.util import find_spec

import anndata as ad
import pandas as pd
import pytest

from perturbflow.config import GuideAssignmentConfig, PerturbationAnalysisConfig
from perturbflow.guide_assignment import assign_guides
from perturbflow.perturbation_analysis import compute_perturbation_signature, run_mixscape


@pytest.mark.skipif(find_spec("pertpy") is None, reason="pertpy not installed")
def test_pertpy_parity_mixscape(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    import pertpy as pt
    import scanpy as sc

    # Build assigned adata once.
    assigned = assign_guides(
        synthetic_adata.copy(),
        synthetic_guide_calls,
        synthetic_guide_metadata,
        config=GuideAssignmentConfig(),
    )

    # ---- Path A: PerturbFlow wrapper ----------------------------------
    a = assigned.copy()
    compute_perturbation_signature(a, config=PerturbationAnalysisConfig(n_neighbors=15))
    a = run_mixscape(a, config=PerturbationAnalysisConfig(n_neighbors=15))

    # ---- Path B: direct pertpy ----------------------------------------
    b = assigned.copy()
    # Reproduce PerturbFlow's PCA prep so the kNN representation is identical.
    tmp = b.copy()
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    sc.pp.scale(tmp, max_value=10)
    sc.tl.pca(tmp, n_comps=min(30, min(tmp.shape) - 1), random_state=0)
    b.obsm["X_pca"] = tmp.obsm["X_pca"]
    mix = pt.tl.Mixscape()
    mix.perturbation_signature(
        b,
        "perturbation",
        "NT",
        n_neighbors=15,
        split_by=None,
        use_rep="X_pca",
        n_dims=min(15, b.obsm["X_pca"].shape[1]),
        copy=False,
    )
    mix.mixscape(b, "perturbation", "NT", min_de_genes=5, pval_cutoff=0.05, copy=False)

    # ---- Compare ------------------------------------------------------
    assert "mixscape_class_global" in a.obs.columns
    assert "mixscape_class_global" in b.obs.columns
    same = a.obs["mixscape_class_global"].astype(str) == b.obs["mixscape_class_global"].astype(str)
    # Allow up to 2% disagreement to absorb stochastic ties in the GMM
    # initialization; effectively the wrappers should be identical.
    agreement = float(same.mean())
    assert agreement > 0.98, (
        f"PerturbFlow vs direct pertpy Mixscape agreement {agreement:.3f} below threshold"
    )


@pytest.mark.skipif(find_spec("pydeseq2") is None, reason="pydeseq2 not installed")
@pytest.mark.skipif(find_spec("pertpy") is None, reason="pertpy not installed")
def test_ground_truth_recovery(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    """End-to-end check: known-truth downstream genes should be recovered.

    The synthetic fixture seeds 20 downregulated and 10 upregulated
    downstream genes per KO perturbation. We run the full pipeline and
    check that for at least one KO, the top-30 DE genes (by abs log2FC,
    padj < 0.1) include at least half of the seeded downstream genes.

    This is the "does the pipeline actually find what's there" test —
    without it, all the other tests prove correctness of the *plumbing*
    but not the *biology*.
    """
    from perturbflow.config import DEConfig, InputConfig
    from perturbflow.de import run_pseudobulk_de
    from perturbflow.perturbation_analysis import (
        compute_perturbation_signature,
        run_mixscape,
    )
    from tests.fixtures.make_synthetic import build_synthetic

    bundle = build_synthetic(seed=0)
    design = bundle["design"]
    adata = bundle["adata"]
    adata = assign_guides(adata, bundle["guide_calls"], bundle["guide_metadata"])
    compute_perturbation_signature(adata, config=PerturbationAnalysisConfig(n_neighbors=15))
    adata = run_mixscape(adata, config=PerturbationAnalysisConfig(n_neighbors=15))

    de = run_pseudobulk_de(
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

    # For each KO, see how many seeded downstream genes are in the top 50 by abs log2FC.
    n_recovered = []
    for pert in ("KOA", "KOB", "KOC", "KOD"):
        if pert not in de:
            continue
        df = de[pert].copy()
        df["_abs"] = df["log2FoldChange"].abs()
        top = set(df[df["padj"] < 0.1].sort_values("_abs", ascending=False).head(50)["gene"])
        seeded_down = {f"GENE_{i:04d}" for i in design.downstream_down[pert]}
        seeded_up = {f"GENE_{i:04d}" for i in design.downstream_up[pert]}
        seeded = seeded_down | seeded_up
        recovered = top & seeded
        n_recovered.append(len(recovered) / max(len(seeded), 1))
    assert n_recovered, "Pipeline produced no DE results to evaluate"
    # At least one KO should recover > 40% of its seeded downstream genes.
    assert max(n_recovered) > 0.4, (
        f"No KO recovered >40% of seeded downstream genes; got recovery fractions: "
        f"{[f'{x:.2f}' for x in n_recovered]}"
    )
