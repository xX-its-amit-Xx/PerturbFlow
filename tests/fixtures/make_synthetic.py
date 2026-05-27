"""Build a small synthetic Perturb-seq AnnData fixture.

Design
------
- 2000 cells, 500 genes.
- 10 guides organized as 4 KO perturbations × 2 guides + 1 NT × 2 guides,
  so the dataset has 5 perturbations total (KOA, KOB, KOC, KOD, NT) —
  matching the spec in the README ("10 guides, 5 perturbations").
- Each KO perturbation has a designated target gene whose expression is
  driven to ~0 in 70% of carrying cells ("KO" arm) and left at baseline in
  the remaining 30% ("escaped" arm). This gives Mixscape something
  meaningful to separate.
- Each KO perturbation also has 20 downregulated and 10 upregulated
  "downstream" genes (multiplicative effect on the Poisson mean) so the
  pseudobulk DE has on-target + downstream hits to detect.
- Counts are drawn from a NegBin so dispersion estimation in DESeq2 works.

The fixture is *deterministic* — seeded numpy RNG — so test failures
don't shift between runs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

N_CELLS = 2000
N_GENES = 500
N_PERTURBATIONS = 4  # plus NT control
GUIDES_PER_PERT = 2
ESCAPE_RATE = 0.30
GUIDE_UMI_MEAN_ASSIGNED = 25
GUIDE_UMI_MEAN_OFF_TARGET = 2  # very low background
BASELINE_MEAN_COUNTS = 1.5  # per gene per cell


@dataclass(frozen=True)
class SyntheticDesign:
    perturbation_names: tuple[str, ...]
    guide_ids: tuple[str, ...]
    guide_to_pert: dict[str, str]
    guide_is_control: dict[str, bool]
    target_gene_index: dict[str, int]
    downstream_down: dict[str, np.ndarray]
    downstream_up: dict[str, np.ndarray]


def build_synthetic(*, seed: int = 0) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    design = _design(rng)

    # First 4 gene names are the KO target symbols so per_perturbation_qc
    # can compute on-target log2FC. Remaining genes get generic names.
    gene_list = ["KOA", "KOB", "KOC", "KOD"] + [f"GENE_{i:04d}" for i in range(4, N_GENES)]
    var_names = pd.Index(gene_list, name="gene")

    # Assign each cell to a perturbation (roughly balanced).
    cell_perts = np.repeat(
        list(design.perturbation_names), N_CELLS // len(design.perturbation_names) + 1
    )[:N_CELLS]
    rng.shuffle(cell_perts)
    # For each cell with a KO perturbation, pick one of the 2 guides for that pert.
    cell_guides = np.empty(N_CELLS, dtype=object)
    for pert in design.perturbation_names:
        guides_for_pert = [g for g, p in design.guide_to_pert.items() if p == pert]
        mask = cell_perts == pert
        cell_guides[mask] = rng.choice(guides_for_pert, size=int(mask.sum()))

    # Escape labels: True means escaped (guide present but no perturbation effect).
    is_escaped = np.zeros(N_CELLS, dtype=bool)
    for pert in design.perturbation_names:
        if pert == "NT":
            continue
        mask = cell_perts == pert
        n = int(mask.sum())
        n_escape = round(ESCAPE_RATE * n)
        escape_idx = rng.choice(np.where(mask)[0], size=n_escape, replace=False)
        is_escaped[escape_idx] = True

    # Build the per-cell Poisson mean matrix.
    base = np.full((N_CELLS, N_GENES), BASELINE_MEAN_COUNTS, dtype=float)
    # Add small per-cell library size variability (lognormal).
    cell_factor = rng.lognormal(mean=0.0, sigma=0.3, size=N_CELLS)
    base *= cell_factor[:, None]
    # Add per-gene baseline variability so DE has a non-trivial null.
    gene_factor = rng.lognormal(mean=0.0, sigma=0.4, size=N_GENES)
    base *= gene_factor[None, :]

    # Apply perturbation effects to truly-perturbed cells only.
    perturbed_mask = (~is_escaped) & (cell_perts != "NT")
    for pert in design.perturbation_names:
        if pert == "NT":
            continue
        idx = np.where(perturbed_mask & (cell_perts == pert))[0]
        if len(idx) == 0:
            continue
        # Knock down the target gene to ~5% of baseline.
        tgt = design.target_gene_index[pert]
        base[idx, tgt] *= 0.05
        # Downstream effects.
        base[np.ix_(idx, design.downstream_down[pert])] *= 0.35
        base[np.ix_(idx, design.downstream_up[pert])] *= 3.0

    # Sample counts. NegBin via Gamma-Poisson: shape = 1/dispersion.
    dispersion = 0.2
    gamma_shape = 1.0 / dispersion
    rates = rng.gamma(shape=gamma_shape, scale=base / gamma_shape)
    counts = rng.poisson(rates).astype(np.int32)

    # Build AnnData.
    obs = pd.DataFrame(
        {
            "expected_perturbation": cell_perts,
            "expected_guide": cell_guides,
            "expected_escaped": is_escaped,
        },
        index=pd.Index([f"CELL_{i:05d}" for i in range(N_CELLS)], name="cell"),
    )
    adata = ad.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=var_names))

    # Build guide call CSV (long format). Each cell gets a high UMI for its true
    # guide and some background calls for unrelated guides (mimicking ambient).
    rows: list[dict[str, Any]] = []
    for i, bc in enumerate(adata.obs_names):
        true_guide = str(cell_guides[i])
        umi = int(rng.poisson(GUIDE_UMI_MEAN_ASSIGNED) + 5)
        rows.append({"cell_barcode": bc, "guide_id": true_guide, "umi_count": umi})
        # Background: a few random low-count calls (well below the default floor of 5).
        n_bg = int(rng.poisson(0.5))
        for bg in rng.choice(design.guide_ids, size=n_bg, replace=False):
            if bg == true_guide:
                continue
            bg_umi = int(rng.poisson(GUIDE_UMI_MEAN_OFF_TARGET))
            if bg_umi > 0:
                rows.append({"cell_barcode": bc, "guide_id": str(bg), "umi_count": bg_umi})
    guide_calls = pd.DataFrame(rows)

    # Guide metadata table.
    guide_metadata = pd.DataFrame(
        {
            "guide_id": list(design.guide_ids),
            "target_gene": [
                design.guide_to_pert[g] if not design.guide_is_control[g] else ""
                for g in design.guide_ids
            ],
            "is_control": [design.guide_is_control[g] for g in design.guide_ids],
        }
    )

    return {
        "adata": adata,
        "guide_calls": guide_calls,
        "guide_metadata": guide_metadata,
        "design": design,
    }


def _design(rng: np.random.Generator) -> SyntheticDesign:
    """Pick guide IDs, perturbation names, target gene indices, downstream gene sets."""
    pert_names = ("KOA", "KOB", "KOC", "KOD", "NT")
    guide_to_pert: dict[str, str] = {}
    guide_is_control: dict[str, bool] = {}
    guide_ids: list[str] = []
    for pert in pert_names:
        for k in range(GUIDES_PER_PERT):
            gid = f"{pert}_g{k + 1}"
            guide_ids.append(gid)
            guide_to_pert[gid] = pert
            guide_is_control[gid] = pert == "NT"

    # Target genes for KO perturbations: rename gene_0001..gene_0004 to match.
    target_gene_index = {
        "KOA": 0,
        "KOB": 1,
        "KOC": 2,
        "KOD": 3,
    }

    # Downstream gene sets: 20 down + 10 up per KO, sampled from gene indices 10..N_GENES.
    candidates = np.arange(10, N_GENES)
    downstream_down: dict[str, np.ndarray] = {}
    downstream_up: dict[str, np.ndarray] = {}
    for pert in ("KOA", "KOB", "KOC", "KOD"):
        picks = rng.choice(candidates, size=30, replace=False)
        downstream_down[pert] = picks[:20]
        downstream_up[pert] = picks[20:]

    return SyntheticDesign(
        perturbation_names=pert_names,
        guide_ids=tuple(guide_ids),
        guide_to_pert=guide_to_pert,
        guide_is_control=guide_is_control,
        target_gene_index=target_gene_index,
        downstream_down=downstream_down,
        downstream_up=downstream_up,
    )


if __name__ == "__main__":
    bundle = build_synthetic(seed=0)
    print(bundle["adata"])
    print("guide_calls:", bundle["guide_calls"].shape)
    print("guide_metadata:")
    print(bundle["guide_metadata"])
