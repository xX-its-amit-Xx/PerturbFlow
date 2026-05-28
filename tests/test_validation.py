from __future__ import annotations

import warnings

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from perturbflow.validation import (
    RawCountsAssumptionError,
    assert_raw_counts,
    multi_guide_concordance,
    warn_if_pseudoreplicates,
)


def _make_count_adata(*, log_normalized: bool = False) -> ad.AnnData:
    rng = np.random.default_rng(0)
    X = rng.poisson(5, size=(100, 50)).astype(np.int32)
    if log_normalized:
        X = np.log1p(X / X.sum(axis=1, keepdims=True).clip(min=1) * 1e4).astype(np.float32)
    return ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"c{i}" for i in range(100)]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(50)]),
    )


def test_assert_raw_counts_accepts_integer() -> None:
    a = _make_count_adata(log_normalized=False)
    assert_raw_counts(a)  # no exception


def test_assert_raw_counts_rejects_log_normalized() -> None:
    a = _make_count_adata(log_normalized=True)
    with pytest.raises(RawCountsAssumptionError, match="does not look like raw counts"):
        assert_raw_counts(a)


def test_assert_raw_counts_rejects_negative() -> None:
    a = _make_count_adata()
    a.X = a.X.astype(np.int32) - 100
    with pytest.raises(RawCountsAssumptionError, match="negative"):
        assert_raw_counts(a)


def test_assert_raw_counts_strict_false_warns() -> None:
    a = _make_count_adata(log_normalized=True)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        assert_raw_counts(a, strict=False)
    assert any("does not look like raw counts" in str(wi.message) for wi in w)


def test_warn_if_pseudoreplicates_emits_when_sample_col_none() -> None:
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        warn_if_pseudoreplicates(sample_col=None, n_pseudo_replicates=3, n_perturbations=10)
    assert any("pseudo-replicate" in str(wi.message).lower() for wi in w)


def test_warn_if_pseudoreplicates_silent_when_sample_col_set() -> None:
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        warn_if_pseudoreplicates(sample_col="donor", n_pseudo_replicates=3, n_perturbations=10)
    assert not any("pseudo-replicate" in str(wi.message).lower() for wi in w)


def test_multi_guide_concordance_high_for_concordant_guides() -> None:
    """Two guides targeting the same gene with similar effect sizes should give r ~ 1."""
    rng = np.random.default_rng(0)
    base = rng.normal(0, 1, size=500)
    de_a = pd.DataFrame(
        {"gene": [f"g{i}" for i in range(500)], "log2FoldChange": base + rng.normal(0, 0.1, 500)}
    )
    de_b = pd.DataFrame(
        {"gene": [f"g{i}" for i in range(500)], "log2FoldChange": base + rng.normal(0, 0.1, 500)}
    )
    meta = pd.DataFrame(
        {
            "guide_id": ["G1_g1", "G1_g2"],
            "target_gene": ["G1", "G1"],
            "is_control": [False, False],
        }
    )
    out = multi_guide_concordance({"G1_g1": de_a, "G1_g2": de_b}, meta)
    assert len(out) == 1
    assert out["pearson_r"].iloc[0] > 0.9


def test_multi_guide_concordance_low_for_discordant_guides() -> None:
    """Two guides with anti-correlated effects should flag low concordance."""
    rng = np.random.default_rng(0)
    base = rng.normal(0, 1, size=500)
    de_a = pd.DataFrame({"gene": [f"g{i}" for i in range(500)], "log2FoldChange": base})
    de_b = pd.DataFrame({"gene": [f"g{i}" for i in range(500)], "log2FoldChange": -base})
    meta = pd.DataFrame(
        {
            "guide_id": ["G1_g1", "G1_g2"],
            "target_gene": ["G1", "G1"],
            "is_control": [False, False],
        }
    )
    out = multi_guide_concordance({"G1_g1": de_a, "G1_g2": de_b}, meta)
    assert out["pearson_r"].iloc[0] < -0.9


def test_multi_guide_concordance_skips_single_guide_targets() -> None:
    de_a = pd.DataFrame({"gene": ["g1"], "log2FoldChange": [1.0]})
    meta = pd.DataFrame({"guide_id": ["G1_g1"], "target_gene": ["G1"], "is_control": [False]})
    out = multi_guide_concordance({"G1_g1": de_a}, meta)
    assert out.empty


def test_multi_guide_concordance_skips_when_too_few_common_genes() -> None:
    """Genes too few in common -> emit a NaN row, don't silently drop."""
    de_a = pd.DataFrame({"gene": [f"g{i}" for i in range(5)], "log2FoldChange": [1.0] * 5})
    de_b = pd.DataFrame({"gene": [f"g{i}" for i in range(5)], "log2FoldChange": [1.0] * 5})
    meta = pd.DataFrame(
        {
            "guide_id": ["G1_g1", "G1_g2"],
            "target_gene": ["G1", "G1"],
            "is_control": [False, False],
        }
    )
    out = multi_guide_concordance({"G1_g1": de_a, "G1_g2": de_b}, meta, min_significant_genes=20)
    assert len(out) == 1
    assert pd.isna(out["pearson_r"].iloc[0])
    assert int(out["n_genes"].iloc[0]) == 5
