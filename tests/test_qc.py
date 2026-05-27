from __future__ import annotations

import anndata as ad
import pandas as pd

from perturbflow.config import GuideAssignmentConfig
from perturbflow.guide_assignment import assign_guides
from perturbflow.qc import per_cell_qc, per_guide_qc, per_perturbation_qc


def _assigned(adata: ad.AnnData, calls: pd.DataFrame, meta: pd.DataFrame) -> ad.AnnData:
    return assign_guides(adata.copy(), calls, meta, config=GuideAssignmentConfig())


def test_per_cell_qc_columns(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    adata = _assigned(synthetic_adata, synthetic_guide_calls, synthetic_guide_metadata)
    qc = per_cell_qc(adata)
    for col in ("cell_barcode", "total_counts", "n_genes", "pct_mito", "passes_qc"):
        assert col in qc.columns
    assert len(qc) == adata.n_obs


def test_per_guide_qc_includes_zero_count_guides(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    """Even a guide called in zero cells should appear in the per-guide table."""
    extra = synthetic_guide_metadata.copy()
    extra = pd.concat(
        [
            extra,
            pd.DataFrame([{"guide_id": "GHOST_g1", "target_gene": "GHOST", "is_control": False}]),
        ],
        ignore_index=True,
    )
    adata = _assigned(synthetic_adata, synthetic_guide_calls, extra)
    qc = per_guide_qc(adata, extra)
    ghost = qc[qc["guide_id"] == "GHOST_g1"]
    assert len(ghost) == 1
    assert int(ghost["n_cells"].iloc[0]) == 0


def test_per_perturbation_qc_on_target_lfc_is_negative(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    """KO perturbations should show strongly negative on-target log2FC."""
    adata = _assigned(synthetic_adata, synthetic_guide_calls, synthetic_guide_metadata)
    qc = per_perturbation_qc(adata, control_label="NT")
    assert set(qc["perturbation"]) == {"KOA", "KOB", "KOC", "KOD"}
    # Even with ~30% escape, mean knockdown should still register clearly < 0.
    for pert in ("KOA", "KOB", "KOC", "KOD"):
        lfc = float(qc.loc[qc["perturbation"] == pert, "on_target_log2fc"].iloc[0])
        assert lfc < -0.5, f"{pert} on_target_log2fc={lfc} (expected < -0.5)"
