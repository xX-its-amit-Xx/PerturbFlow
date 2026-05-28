"""Tests for the doublet-detection hook."""

from __future__ import annotations

import warnings

import anndata as ad
import numpy as np
import pandas as pd

from perturbflow.qc import detect_doublets


def _make_count_adata(n: int = 200) -> ad.AnnData:
    rng = np.random.default_rng(0)
    X = rng.poisson(2, size=(n, 200)).astype(np.int32)
    return ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(200)]),
    )


def test_detect_doublets_writes_columns_with_explicit_threshold() -> None:
    """Explicit threshold skips the auto-detection so skimage isn't required."""
    adata = _make_count_adata()
    out = detect_doublets(
        adata,
        expected_doublet_rate=0.05,
        random_state=0,
        threshold=0.25,
    )
    assert "doublet_score" in out.obs
    assert "predicted_doublet" in out.obs
    scores = out.obs["doublet_score"].astype(float)
    # All scores should be in [0, 1]; sanity check on the predicted column dtype.
    assert (scores >= 0).all() and (scores <= 1).all()
    assert out.obs["predicted_doublet"].astype(bool).dtype == bool


def test_detect_doublets_falls_back_when_skimage_missing(monkeypatch) -> None:
    """When auto-threshold isn't available, fall back with a loud warning."""
    import perturbflow.qc as qc

    real_find = qc.find_spec

    def fake_find_spec(name: str):
        return None if name == "skimage" else real_find(name)

    monkeypatch.setattr(qc, "find_spec", fake_find_spec)
    adata = _make_count_adata()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        out = detect_doublets(adata)
    assert any("scikit-image is not installed" in str(wi.message) for wi in w)
    assert "doublet_score" in out.obs
    # Fallback runs the actual scrublet; columns are populated.
    assert out.obs["doublet_score"].notna().any()
