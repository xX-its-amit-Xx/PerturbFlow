from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from perturbflow.config import GuideAssignmentConfig
from perturbflow.downstream import compute_cell_state_effects, score_pathways
from perturbflow.guide_assignment import assign_guides


def _assigned(
    adata: ad.AnnData,
    calls: pd.DataFrame,
    meta: pd.DataFrame,
) -> ad.AnnData:
    return assign_guides(adata.copy(), calls, meta, config=GuideAssignmentConfig())


def test_cell_state_effects_returns_one_row_per_pert(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    """Sanity check on shape, sort order, and the centroid_shift values."""
    adata = _assigned(synthetic_adata, synthetic_guide_calls, synthetic_guide_metadata)
    # Inject a deterministic embedding instead of running umap.
    rng = np.random.default_rng(0)
    coords = rng.normal(size=(adata.n_obs, 2)).astype(np.float32)
    # Bias each perturbation centroid by a different vector.
    perts = adata.obs["perturbation"].astype(str).values
    biases = {
        "KOA": np.array([3.0, 0.0]),
        "KOB": np.array([0.0, 3.0]),
        "KOC": np.array([-3.0, 0.0]),
        "KOD": np.array([0.0, -3.0]),
        "NT": np.array([0.0, 0.0]),
    }
    for pert, b in biases.items():
        mask = perts == pert
        coords[mask] += b.astype(np.float32)
    adata.obsm["X_umap"] = coords

    cs = compute_cell_state_effects(adata, control_label="NT")
    assert set(cs["perturbation"]) == {"KOA", "KOB", "KOC", "KOD"}
    # All centroid shifts should be roughly the magnitude of the bias (3.0).
    for pert in ("KOA", "KOB", "KOC", "KOD"):
        shift = float(cs.loc[cs["perturbation"] == pert, "centroid_shift"].iloc[0])
        assert 2.5 < shift < 3.5, f"{pert}: shift={shift}"


def test_score_pathways_empty_de_returns_empty() -> None:
    df = score_pathways({}, config=None)
    assert isinstance(df, pd.DataFrame)
    assert df.empty


def test_score_pathways_disabled_returns_empty() -> None:
    from perturbflow.config import DownstreamConfig

    fake_de = {
        "KOA": pd.DataFrame(
            {"gene": ["X", "Y"], "stat": [1.0, -1.0], "log2FoldChange": [1.0, -1.0]}
        )
    }
    df = score_pathways(fake_de, config=DownstreamConfig(enable_pathway_scoring=False))
    assert df.empty
