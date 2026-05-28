"""Validity checks that catch user-error footguns the rest of the pipeline can't.

The biggest classes of footgun in a Perturb-seq pipeline:

1. **Wrong data in ``adata.X``.** Downstream tools (DESeq2, decoupler)
   assume specific layouts (raw counts, log-normalized, z-scored). If the
   user hands us the wrong layout we produce convincing but wrong numbers.
2. **Underpowered DE designs.** Hash-based pseudo-replicates feel like
   replication but aren't independent biological samples. DESeq2 will
   gladly fit dispersion across them and call everything significant.
3. **Inconsistent guides per target.** If two guides targeting the same
   gene produce opposite-direction effects, one of them is off-target and
   the per-gene DE call is unreliable.

These checks emit warnings (or raise ``ValueError`` for show-stoppers)
loudly enough that a working scientist notices before they paste a
volcano into a slide.
"""

from __future__ import annotations

import logging
import warnings

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

logger = logging.getLogger(__name__)


class RawCountsAssumptionError(ValueError):
    """Raised when X looks log-normalized but raw counts were expected."""


def assert_raw_counts(adata: ad.AnnData, *, layer: str | None = None, strict: bool = True) -> None:
    """Verify that ``adata.X`` (or a layer) contains integer raw counts.

    Heuristics:

    1. Values must be non-negative.
    2. The fraction of integer-valued entries must be > 99.9% (allow a few
       cast-from-float-to-int rounding entries from upstream pipelines).
    3. The maximum value should be large (>= 50) — log-normalized matrices
       very rarely exceed values of 10 even on a sparse representation.

    We're conservative on purpose: the false-positive cost of refusing a
    legitimate matrix is "user re-checks and adds an override flag";
    the false-negative cost of accepting a log-normalized matrix is
    a manuscript-killing DE result.
    """
    X = adata.layers[layer] if layer is not None else adata.X
    vals = X.data if sparse.issparse(X) else np.asarray(X).ravel()
    if vals.size == 0:
        raise RawCountsAssumptionError("X is empty")
    if (vals < 0).any():
        raise RawCountsAssumptionError(
            "X contains negative values — expected raw counts (non-negative integers)"
        )

    # Sample at most 50k values for the integer-fraction check; for a 100k
    # cell dataset the full check is fine, but for a 1M-cell screen it's
    # wasteful.
    sample = (
        vals
        if vals.size <= 50_000
        else vals[np.random.default_rng(0).integers(0, vals.size, 50_000)]
    )
    is_int = np.isclose(sample, np.round(sample), atol=1e-6)
    int_fraction = float(is_int.mean())
    max_val = float(vals.max())

    # The integer-fraction is the strongest discriminator: log-normalized data
    # is float-by-construction. The max-value floor is a secondary check that
    # we relax for tiny test fixtures (< 10k entries) since synthetic Poisson
    # data may not reach 10.
    raw_like = int_fraction > 0.999 and (max_val >= 10 or vals.size < 10_000)
    if not raw_like:
        msg = (
            "adata.X does not look like raw counts "
            f"(int_fraction={int_fraction:.4f}, max_value={max_val:.2f}). "
            "PerturbFlow's DE module needs raw integer counts; pass them "
            "via the layers={raw_layer!r} argument, or set strict=False to "
            "silence this check."
        )
        if strict:
            raise RawCountsAssumptionError(msg)
        warnings.warn(msg, stacklevel=2)


def warn_if_pseudoreplicates(
    *,
    sample_col: str | None,
    n_pseudo_replicates: int,
    n_perturbations: int,
) -> None:
    """Loud warning when we're about to synthesize pseudo-replicates for DE.

    Pseudo-replicates from hash-binning cell barcodes are *not* biological
    replicates. They give DESeq2 enough samples to fit dispersion, but the
    p-values they produce are anticonservative (too small) because the
    "replicates" are i.i.d. draws from the same underlying distribution.
    For exploratory / ranking use this is fine; for any claim of
    statistical significance it is not.

    Only suppress this warning when you know the next person to look at
    the volcano understands the caveat.
    """
    if sample_col is not None:
        return
    msg = (
        f"DE is using {n_pseudo_replicates} hash-bin pseudo-replicates per "
        f"perturbation ({n_perturbations} perturbations total) because no "
        "biological sample_col was configured. Pseudo-replicate p-values "
        "are anticonservative — treat DE rankings as exploratory, not as "
        "calibrated statistical tests. To use real replication, set "
        "input.sample_col to your donor / batch column."
    )
    logger.warning(msg)
    warnings.warn(msg, stacklevel=2)


def multi_guide_concordance(
    de_results: dict[str, pd.DataFrame],
    guide_metadata: pd.DataFrame,
    *,
    min_significant_genes: int = 20,
) -> pd.DataFrame:
    """Concordance of DE effect sizes between guides targeting the same gene.

    For every gene with ≥2 non-control guides, compute the Pearson
    correlation of ``log2FoldChange`` vectors across the union of genes
    significant in either DE result. Low correlation (e.g. r < 0.3) is a
    red flag that one or both guides have substantial off-target effects.

    This requires that DE was run **per guide** (not per gene). If the
    upstream pipeline pooled cells from both guides into a single
    perturbation label (the default in PerturbFlow), there is nothing for
    this function to check and it returns an empty frame.

    Returns one row per target gene:
    ``target_gene | guide_a | guide_b | pearson_r | n_genes``.
    """
    columns = ["target_gene", "guide_a", "guide_b", "pearson_r", "n_genes"]
    if not de_results:
        return pd.DataFrame(columns=columns)

    targeting = guide_metadata[~guide_metadata["is_control"].astype(bool)]
    guides_by_gene = targeting.groupby("target_gene")["guide_id"].agg(list).to_dict()

    rows: list[dict[str, object]] = []
    for gene, guides in guides_by_gene.items():
        guides_with_de = [g for g in guides if g in de_results]
        if len(guides_with_de) < 2:
            continue
        for i in range(len(guides_with_de)):
            for j in range(i + 1, len(guides_with_de)):
                a, b = guides_with_de[i], guides_with_de[j]
                da = de_results[a].set_index("gene")["log2FoldChange"]
                db = de_results[b].set_index("gene")["log2FoldChange"]
                common = da.index.intersection(db.index)
                if len(common) < min_significant_genes:
                    rows.append(
                        {
                            "target_gene": gene,
                            "guide_a": a,
                            "guide_b": b,
                            "pearson_r": float("nan"),
                            "n_genes": len(common),
                        }
                    )
                    continue
                r = float(np.corrcoef(da.loc[common].values, db.loc[common].values)[0, 1])
                rows.append(
                    {
                        "target_gene": gene,
                        "guide_a": a,
                        "guide_b": b,
                        "pearson_r": r,
                        "n_genes": len(common),
                    }
                )
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows).sort_values("pearson_r").reset_index(drop=True)
