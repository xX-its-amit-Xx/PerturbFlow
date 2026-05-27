"""Sanity-check the synthetic fixture itself.

If the fixture is broken, every downstream test fails for the wrong reason.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from tests.fixtures.make_synthetic import (
    ESCAPE_RATE,
    GUIDES_PER_PERT,
    N_CELLS,
    N_GENES,
    build_synthetic,
)


def test_fixture_shape() -> None:
    bundle = build_synthetic(seed=0)
    adata: ad.AnnData = bundle["adata"]
    assert adata.shape == (N_CELLS, N_GENES)
    assert adata.X.dtype == np.int32
    assert (adata.X >= 0).all()


def test_fixture_guide_count() -> None:
    bundle = build_synthetic(seed=0)
    meta: pd.DataFrame = bundle["guide_metadata"]
    # 5 perturbations × 2 guides each = 10 guides.
    assert len(meta) == 10
    assert meta["is_control"].sum() == GUIDES_PER_PERT
    # 4 KO perturbations are non-control.
    assert (~meta["is_control"]).sum() == 8


def test_fixture_target_genes_present() -> None:
    """KOA..KOD must be real gene symbols so on_target_lfc has something to find."""
    bundle = build_synthetic(seed=0)
    adata = bundle["adata"]
    for target in ("KOA", "KOB", "KOC", "KOD"):
        assert target in adata.var_names


def test_fixture_escape_rate_reasonable() -> None:
    """The configured escape rate should be respected within tolerance."""
    bundle = build_synthetic(seed=0)
    adata = bundle["adata"]
    ko_mask = adata.obs["expected_perturbation"] != "NT"
    escape_actual = adata.obs.loc[ko_mask, "expected_escaped"].mean()
    assert abs(escape_actual - ESCAPE_RATE) < 0.05


def test_fixture_knockdown_in_perturbed_cells() -> None:
    """In truly perturbed cells, the target gene count should be ~zero."""
    bundle = build_synthetic(seed=0)
    adata = bundle["adata"]
    obs = adata.obs
    for pert, target_gene in (("KOA", "KOA"), ("KOB", "KOB")):
        idx_ko = np.where((obs["expected_perturbation"] == pert) & (~obs["expected_escaped"]))[0]
        idx_nt = np.where(obs["expected_perturbation"] == "NT")[0]
        target_idx = adata.var_names.get_loc(target_gene)
        ko_mean = float(adata.X[idx_ko, target_idx].mean())
        nt_mean = float(adata.X[idx_nt, target_idx].mean())
        # KO arm should be at least 5× lower than NT mean.
        assert ko_mean < nt_mean / 5, f"{pert}: ko={ko_mean} nt={nt_mean}"


def test_fixture_deterministic() -> None:
    a = build_synthetic(seed=0)["adata"]
    b = build_synthetic(seed=0)["adata"]
    assert np.array_equal(a.X, b.X)
