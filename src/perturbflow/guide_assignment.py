"""Guide assignment for Perturb-seq.

The job here is to take a long-format guide-call table — one row per (cell,
guide) above some detection floor — and produce a per-cell call.

We use an opinionated three-rule procedure:

1. **UMI floor.** A guide must have at least ``min_guide_umi`` reads in a
   cell to be considered called. Guide capture is shallow; without a floor
   you get noise calls.
2. **Dominance.** For singly-assigned cells we require the top guide to be
   at least ``dominance_ratio``× more abundant than the second-best guide.
   This catches "ambiguous" cells where two guides are nearly tied — a
   common artifact of guide-library cross-contamination.
3. **Multi-guide cap.** Cells with more than ``max_guides`` guides above
   the UMI floor are flagged as multi-guide (typically doublets in
   pooled screens). The default ``max_guides=1`` reflects the usual MOI~0.3
   design.

Each cell ends up labelled with exactly one of:

- ``assigned``     — single guide passes all rules
- ``unassigned``   — no guide above ``min_guide_umi``
- ``ambiguous``    — top guide does not dominate second-best
- ``multi-guide``  — too many guides above floor (likely doublet)
"""

from __future__ import annotations

import logging
from typing import Literal

import anndata as ad
import numpy as np
import pandas as pd

from perturbflow.config import GuideAssignmentConfig

logger = logging.getLogger(__name__)


AssignmentStatus = Literal["assigned", "unassigned", "ambiguous", "multi-guide"]


def assign_guides(
    adata: ad.AnnData,
    guide_calls: pd.DataFrame,
    guide_metadata: pd.DataFrame,
    *,
    config: GuideAssignmentConfig | None = None,
) -> ad.AnnData:
    """Assign each cell to at most one guide and the corresponding perturbation.

    Modifies a copy of ``adata`` in-place semantics — returns the modified
    AnnData. Cells failing assignment are kept by default with their status
    recorded; set ``config.drop_unassigned=True`` to filter them out at the
    end.

    The following columns are added to ``obs``:

    - ``guide`` (str | NaN) — the single called guide_id, or NaN if not assigned
    - ``guide_umi`` (int) — UMI count of the called guide (0 if unassigned)
    - ``second_guide_umi`` (int) — UMI count of the runner-up
    - ``perturbation`` (str) — target_gene from the guide metadata, or the
      control_label for non-targeting guides, or "unassigned"
    - ``is_control`` (bool) — True iff the called guide is non-targeting
    - ``assignment_status`` (str) — one of the four AssignmentStatus values

    Parameters
    ----------
    adata
        Gene-expression AnnData. Cells not appearing in the guide-call table
        are still kept (as unassigned) so that QC can see total cell counts.
    guide_calls
        Long-format guide call table (see :func:`perturbflow.io.read_guide_calls`).
    guide_metadata
        Guide -> target gene mapping (see :func:`perturbflow.io.read_guide_metadata`).
    config
        Knobs for the assignment rules. Defaults if None.
    """
    cfg = config or GuideAssignmentConfig()

    # 1. Filter calls below the UMI floor.
    above = guide_calls.loc[guide_calls["umi_count"] >= cfg.min_guide_umi].copy()

    # 2. Validate all guide IDs against metadata.
    known_guides = set(guide_metadata["guide_id"])
    unknown = set(above["guide_id"]) - known_guides
    if unknown:
        # We treat unknown guides as a hard error rather than silently dropping —
        # it's almost always a guide-library version mismatch.
        raise ValueError(
            f"Guide call table contains {len(unknown)} guide_ids not present in "
            f"guide metadata (e.g. {sorted(unknown)[:3]}). Check that you loaded "
            "the right guide library version."
        )

    # 3. For each cell, find the top and runner-up guides.
    above = above.sort_values(
        ["cell_barcode", "umi_count"], ascending=[True, False], kind="mergesort"
    )
    per_cell = above.groupby("cell_barcode", sort=False).agg(
        n_guides=("guide_id", "size"),
        top_guide=("guide_id", "first"),
        top_umi=("umi_count", "first"),
        second_umi=("umi_count", lambda x: int(x.iloc[1]) if len(x) > 1 else 0),
    )

    # 4. Classify.
    status = pd.Series("assigned", index=per_cell.index, dtype=object)
    multi_mask = per_cell["n_guides"] > cfg.max_guides
    status[multi_mask] = "multi-guide"
    # Ambiguity only matters for cells that would otherwise be 'assigned'
    # (single-guide). For multi-guide cells the dominance ratio is irrelevant —
    # they're rejected by the multi-guide rule.
    second_floor = np.maximum(per_cell["second_umi"], 1)  # avoid /0
    ambig_mask = (
        (~multi_mask)
        & (per_cell["second_umi"] > 0)
        & ((per_cell["top_umi"] / second_floor) < cfg.dominance_ratio)
    )
    status[ambig_mask] = "ambiguous"
    per_cell["assignment_status"] = status

    # 5. Reattach to AnnData, including cells with NO guide above floor.
    obs = adata.obs.copy()
    obs_index_str = obs.index.astype(str)
    obs.index = obs_index_str

    joined = per_cell.reindex(obs_index_str)
    joined["assignment_status"] = joined["assignment_status"].fillna("unassigned")
    joined["top_guide"] = joined["top_guide"].astype("object")
    joined.loc[joined["assignment_status"] != "assigned", "top_guide"] = pd.NA

    obs["guide"] = joined["top_guide"]
    obs["guide_umi"] = joined["top_umi"].fillna(0).astype(int)
    obs["second_guide_umi"] = joined["second_umi"].fillna(0).astype(int)
    obs["assignment_status"] = joined["assignment_status"].astype("category")

    # 6. Map guide -> perturbation. Control guides collapse to a single label
    # (the perturbation_analysis module reads this label as the control set).
    guide_to_pert: dict[str, str] = {}
    guide_to_control: dict[str, bool] = {}
    for _, row in guide_metadata.iterrows():
        if bool(row["is_control"]):
            guide_to_pert[row["guide_id"]] = "NT"
            guide_to_control[row["guide_id"]] = True
        else:
            tgt = str(row["target_gene"]).strip()
            if not tgt:
                raise ValueError(
                    f"Guide {row['guide_id']!r} is not a control but has no target_gene"
                )
            guide_to_pert[row["guide_id"]] = tgt
            guide_to_control[row["guide_id"]] = False

    obs["perturbation"] = obs["guide"].map(guide_to_pert).fillna("unassigned").astype("category")
    obs["is_control"] = obs["guide"].map(guide_to_control).fillna(False).astype(bool)

    adata = adata.copy()
    adata.obs = obs

    # 7. Record provenance for the report.
    adata.uns.setdefault("perturbflow", {})
    adata.uns["perturbflow"]["guide_assignment"] = {
        "min_guide_umi": cfg.min_guide_umi,
        "dominance_ratio": cfg.dominance_ratio,
        "max_guides": cfg.max_guides,
        "n_cells_total": int(adata.n_obs),
        "n_cells_assigned": int((obs["assignment_status"] == "assigned").sum()),
        "n_cells_unassigned": int((obs["assignment_status"] == "unassigned").sum()),
        "n_cells_ambiguous": int((obs["assignment_status"] == "ambiguous").sum()),
        "n_cells_multi_guide": int((obs["assignment_status"] == "multi-guide").sum()),
    }
    logger.info(
        "Guide assignment: %d/%d assigned (%.1f%%), %d unassigned, %d ambiguous, %d multi-guide",
        adata.uns["perturbflow"]["guide_assignment"]["n_cells_assigned"],
        adata.n_obs,
        100.0
        * adata.uns["perturbflow"]["guide_assignment"]["n_cells_assigned"]
        / max(adata.n_obs, 1),
        adata.uns["perturbflow"]["guide_assignment"]["n_cells_unassigned"],
        adata.uns["perturbflow"]["guide_assignment"]["n_cells_ambiguous"],
        adata.uns["perturbflow"]["guide_assignment"]["n_cells_multi_guide"],
    )

    if cfg.drop_unassigned:
        keep = obs["assignment_status"] == "assigned"
        adata = adata[keep.values].copy()
        logger.info("Dropped non-assigned cells: %d remaining", adata.n_obs)

    return adata
