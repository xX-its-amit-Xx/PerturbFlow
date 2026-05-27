from __future__ import annotations

import anndata as ad
import pandas as pd
import pytest

from perturbflow.config import GuideAssignmentConfig
from perturbflow.guide_assignment import assign_guides


def test_assign_guides_majority_assigned(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    adata = assign_guides(
        synthetic_adata.copy(),
        synthetic_guide_calls,
        synthetic_guide_metadata,
        config=GuideAssignmentConfig(drop_unassigned=False),
    )
    # With UMI mean of 25 and threshold of 5, almost all cells should be assigned.
    assigned_frac = (adata.obs["assignment_status"] == "assigned").mean()
    assert assigned_frac > 0.90


def test_assign_guides_recovers_truth(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    """The assigned guide should match the synthetic ground truth for assigned cells."""
    adata = assign_guides(
        synthetic_adata.copy(),
        synthetic_guide_calls,
        synthetic_guide_metadata,
    )
    assigned = adata.obs["assignment_status"] == "assigned"
    truth = adata.obs.loc[assigned, "expected_guide"].astype(str)
    called = adata.obs.loc[assigned, "guide"].astype(str)
    # Background calls are well below the UMI floor, so concordance should be very high.
    assert (truth == called).mean() > 0.95


def test_assign_guides_perturbation_label(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    adata = assign_guides(
        synthetic_adata.copy(),
        synthetic_guide_calls,
        synthetic_guide_metadata,
    )
    assert set(adata.obs["perturbation"].astype(str).unique()) == {
        "KOA",
        "KOB",
        "KOC",
        "KOD",
        "NT",
    }
    # is_control matches the NT label exactly.
    nt_mask = adata.obs["perturbation"].astype(str) == "NT"
    assert nt_mask.equals(adata.obs["is_control"])


def test_assign_guides_unknown_guide_errors(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    bad_calls = pd.concat(
        [
            synthetic_guide_calls,
            pd.DataFrame(
                [{"cell_barcode": "CELL_00000", "guide_id": "MYSTERY_g1", "umi_count": 999}]
            ),
        ],
        ignore_index=True,
    )
    with pytest.raises(ValueError, match="guide_ids not present"):
        assign_guides(synthetic_adata.copy(), bad_calls, synthetic_guide_metadata)


def test_assign_guides_multi_guide_detection(
    synthetic_adata: ad.AnnData,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    """When two guides are above threshold, cell should be multi-guide."""
    calls = pd.DataFrame(
        [
            {"cell_barcode": "CELL_00000", "guide_id": "KOA_g1", "umi_count": 20},
            {"cell_barcode": "CELL_00000", "guide_id": "KOB_g1", "umi_count": 20},
        ]
    )
    adata = assign_guides(
        synthetic_adata[:5].copy(),
        calls,
        synthetic_guide_metadata,
        config=GuideAssignmentConfig(drop_unassigned=False, max_guides=1),
    )
    assert adata.obs.loc["CELL_00000", "assignment_status"] == "multi-guide"


def test_assign_guides_ambiguous_when_no_dominance(
    synthetic_adata: ad.AnnData,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    """When top guide does not dominate runner-up by ratio, status is ambiguous."""
    calls = pd.DataFrame(
        [
            {"cell_barcode": "CELL_00001", "guide_id": "KOA_g1", "umi_count": 10},
            {"cell_barcode": "CELL_00001", "guide_id": "KOB_g1", "umi_count": 9},
        ]
    )
    adata = assign_guides(
        synthetic_adata[:5].copy(),
        calls,
        synthetic_guide_metadata,
        config=GuideAssignmentConfig(drop_unassigned=False, max_guides=2, dominance_ratio=2.0),
    )
    assert adata.obs.loc["CELL_00001", "assignment_status"] == "ambiguous"


def test_assign_guides_drops_unassigned(
    synthetic_adata: ad.AnnData,
    synthetic_guide_calls: pd.DataFrame,
    synthetic_guide_metadata: pd.DataFrame,
) -> None:
    kept = assign_guides(
        synthetic_adata.copy(),
        synthetic_guide_calls,
        synthetic_guide_metadata,
        config=GuideAssignmentConfig(drop_unassigned=True),
    )
    assert (kept.obs["assignment_status"] == "assigned").all()
    # Provenance should still reflect pre-drop totals.
    prov = kept.uns["perturbflow"]["guide_assignment"]
    assert prov["n_cells_assigned"] <= prov["n_cells_total"]
