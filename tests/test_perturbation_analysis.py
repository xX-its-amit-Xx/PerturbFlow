from __future__ import annotations

import anndata as ad
import pandas as pd

from perturbflow.config import GuideAssignmentConfig, PerturbationAnalysisConfig
from perturbflow.guide_assignment import assign_guides
from perturbflow.perturbation_analysis import (
    compute_perturbation_signature,
    run_mixscape,
)


def _assigned(
    adata: ad.AnnData,
    calls: pd.DataFrame,
    meta: pd.DataFrame,
) -> ad.AnnData:
    return assign_guides(adata.copy(), calls, meta, config=GuideAssignmentConfig())


def test_perturbation_signature_layer_added(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    adata = _assigned(synthetic_adata, synthetic_guide_calls, synthetic_guide_metadata)
    compute_perturbation_signature(adata, config=PerturbationAnalysisConfig(n_neighbors=10))
    assert "perturbation_signature" in adata.layers
    assert adata.layers["perturbation_signature"].shape == adata.shape


def test_mixscape_separates_ko_from_np(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    adata = _assigned(synthetic_adata, synthetic_guide_calls, synthetic_guide_metadata)
    cfg = PerturbationAnalysisConfig(n_neighbors=15)
    adata = run_mixscape(adata, config=cfg)
    assert "mixscape_perturbed" in adata.obs.columns

    # Cells we engineered as escaped should mostly be classified NP (mixscape_perturbed == False).
    # Restrict to non-NT cells.
    ko_mask = adata.obs["expected_perturbation"] != "NT"
    sub = adata.obs[ko_mask]
    truly_escaped = sub["expected_escaped"].astype(bool)
    classified_perturbed = sub["mixscape_perturbed"].astype(bool)

    # Among truly-escaped cells, the majority should NOT be flagged as perturbed.
    if truly_escaped.any():
        np_rate_in_escaped = 1.0 - classified_perturbed[truly_escaped].mean()
        # Mixscape isn't perfect; we ask only for better-than-coin.
        assert np_rate_in_escaped > 0.55, f"Mixscape NP rate in escaped: {np_rate_in_escaped:.2f}"

    # Among truly-perturbed cells, the majority should be flagged as perturbed.
    truly_perturbed = ~truly_escaped
    if truly_perturbed.any():
        ko_rate_in_perturbed = classified_perturbed[truly_perturbed].mean()
        assert ko_rate_in_perturbed > 0.55, (
            f"Mixscape KO rate in perturbed: {ko_rate_in_perturbed:.2f}"
        )


def test_mixscape_requires_control(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    """If perturbation key has no cells matching the control_label, signature errors."""
    adata = _assigned(synthetic_adata, synthetic_guide_calls, synthetic_guide_metadata)
    cfg = PerturbationAnalysisConfig(control_label="DOES_NOT_EXIST")
    import pytest

    with pytest.raises(ValueError, match="absent"):
        compute_perturbation_signature(adata, config=cfg)
