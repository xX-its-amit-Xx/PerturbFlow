"""Mixscape-style perturbation signal detection.

This module wraps :mod:`pertpy`'s implementation of the Mixscape procedure
(Papalexi *et al.* 2021) when it's available, and falls back to an
equivalent two-component Gaussian mixture on a local perturbation signature
when it is not. The fallback exists so the pipeline still runs in an
``--extra mixscape``-free install, but for production work install pertpy.

Why Mixscape matters
--------------------
In any pooled CRISPR screen a substantial fraction of cells receive the
guide RNA but fail to lose the target protein — incomplete editing, in-frame
indels, escape via the alternate allele. These "escaped" cells dilute the
DE signal toward the null. Mixscape models the per-cell perturbation
signature as a mixture of "knock-out" and "non-perturbed" components and
labels each cell, letting downstream DE drop the escaped cells from the
treatment arm.

The output column ``mixscape_class`` takes values:

- ``"<target> KO"``     — cell is in the perturbed component
- ``"<target> NP"``     — cell received the guide but is non-perturbed (escaped)
- ``"<control_label>"`` — cell carries a non-targeting guide
"""

from __future__ import annotations

import logging
from importlib.util import find_spec

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from perturbflow.config import PerturbationAnalysisConfig

logger = logging.getLogger(__name__)


def compute_perturbation_signature(
    adata: ad.AnnData,
    *,
    config: PerturbationAnalysisConfig | None = None,
    perturbation_key: str = "perturbation",
    layer: str | None = None,
) -> ad.AnnData:
    """Compute the per-cell perturbation signature.

    Implements the Mixscape "local perturbation signature": for each
    perturbed cell, subtract the average expression of its ``n_neighbors``
    nearest control cells (in PCA space). The signature isolates the
    perturbation-attributable expression change from background variability
    across cells.

    After this runs, ``adata.layers["X_pert"]`` (pertpy's convention) and
    ``adata.layers["perturbation_signature"]`` (alias for back-compat) both
    contain a dense (cells × genes) matrix of log-normalized residuals.

    We log-normalize and PCA-embed if those steps haven't already been done,
    because the signature is meaningful only on a comparable space. If you
    have your own preferred preprocessing, run it first and pass
    ``layer="your_layer_name"`` to compute the signature on that.
    """
    cfg = config or PerturbationAnalysisConfig()
    if perturbation_key not in adata.obs.columns:
        raise KeyError(
            f"adata.obs[{perturbation_key!r}] not found — did you forget to run assign_guides?"
        )
    if cfg.control_label not in adata.obs[perturbation_key].astype(str).unique():
        raise ValueError(
            f"Control label {cfg.control_label!r} absent from "
            f"adata.obs[{perturbation_key!r}] — at least one control cell is required."
        )

    if find_spec("pertpy") is not None:
        import pertpy as pt

        # pertpy expects a PCA representation in adata.obsm['X_pca']; if absent,
        # compute one on a log-normalized copy so the kNN signature is meaningful.
        _ensure_pca(adata, layer=layer)
        mix = pt.tl.Mixscape()
        mix.perturbation_signature(
            adata,
            perturbation_key,
            cfg.control_label,
            n_neighbors=cfg.n_neighbors,
            split_by=None,
            use_rep="X_pca",
            n_dims=min(15, adata.obsm["X_pca"].shape[1]),
            copy=False,
        )
        # pertpy writes to layers['X_pert']; expose under the perturbflow name too.
        if "X_pert" in adata.layers:
            adata.layers["perturbation_signature"] = adata.layers["X_pert"]
        logger.info(
            "Computed perturbation signature via pertpy.tl.Mixscape (k=%d neighbors)",
            cfg.n_neighbors,
        )
        return adata

    logger.warning(
        "pertpy is not installed; falling back to an in-house perturbation signature. "
        "Install with `pip install perturbflow[mixscape]` for the validated implementation."
    )
    _fallback_perturbation_signature(adata, cfg=cfg, perturbation_key=perturbation_key, layer=layer)
    return adata


def run_mixscape(
    adata: ad.AnnData,
    *,
    config: PerturbationAnalysisConfig | None = None,
    perturbation_key: str = "perturbation",
) -> ad.AnnData:
    """Run Mixscape KO vs NP classification.

    Requires :func:`compute_perturbation_signature` to have populated
    ``adata.layers["perturbation_signature"]`` (we call it for you if it
    is missing).

    Adds ``adata.obs["mixscape_class"]`` and a boolean
    ``adata.obs["mixscape_perturbed"]`` (True for KO, False for NP/NT) that
    the DE module reads when ``de.use_mixscape_filter`` is True.
    """
    cfg = config or PerturbationAnalysisConfig()
    if "perturbation_signature" not in adata.layers:
        compute_perturbation_signature(adata, config=cfg, perturbation_key=perturbation_key)

    if find_spec("pertpy") is not None:
        import pertpy as pt

        mix = pt.tl.Mixscape()
        # pertpy's mixscape defaults to reading layers["X_pert"], which our
        # compute_perturbation_signature populates (see :func:`compute_perturbation_signature`).
        mix.mixscape(
            adata,
            perturbation_key,
            cfg.control_label,
            min_de_genes=5,
            pval_cutoff=cfg.mixscape_pval_cutoff,
            copy=False,
        )
        # pertpy writes 'mixscape_class' (e.g. "GENE_A KO") and 'mixscape_class_global'
        # ('KO' | 'NP' | '<control_label>'). Standardize to a single boolean column.
        if "mixscape_class_global" in adata.obs.columns:
            adata.obs["mixscape_perturbed"] = adata.obs["mixscape_class_global"].astype(str) == "KO"
        else:  # pragma: no cover - defensive against future pertpy schema drift
            adata.obs["mixscape_perturbed"] = (
                adata.obs.get("mixscape_class", pd.Series(index=adata.obs_names, dtype=str))
                .astype(str)
                .str.endswith(" KO")
            )
        logger.info(
            "Mixscape: %d KO cells, %d NP cells, %d control",
            int(adata.obs["mixscape_perturbed"].sum()),
            int((adata.obs["mixscape_class_global"] == "NP").sum())
            if "mixscape_class_global" in adata.obs
            else -1,
            int((adata.obs[perturbation_key].astype(str) == cfg.control_label).sum()),
        )
        return adata

    _fallback_mixscape(adata, cfg=cfg, perturbation_key=perturbation_key)
    return adata


def _ensure_pca(adata: ad.AnnData, *, layer: str | None) -> None:
    """Make sure ``adata.obsm['X_pca']`` exists, computing one on a log-normalized copy.

    We compute on a copy so the caller's ``X`` stays as raw counts (the DE
    module depends on raw counts being in ``X``).
    """
    if "X_pca" in adata.obsm:
        return
    tmp = adata.copy()
    if layer is None:
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
    sc.pp.scale(tmp, max_value=10)
    n_comps = min(30, min(tmp.shape) - 1)
    sc.tl.pca(tmp, n_comps=n_comps, random_state=0)
    adata.obsm["X_pca"] = tmp.obsm["X_pca"]


# ----------------------------------------------------------------------
# Fallback implementation (used only when pertpy is not installed).
# ----------------------------------------------------------------------


def _fallback_perturbation_signature(
    adata: ad.AnnData,
    *,
    cfg: PerturbationAnalysisConfig,
    perturbation_key: str,
    layer: str | None,
) -> None:
    """In-house perturbation signature: subtract mean of k-NN controls.

    Computes a PCA-space k-NN restricted to control cells, then for every
    cell subtracts the mean log-normalized expression of its nearest
    controls. Equivalent to pertpy's :meth:`Mixscape.perturbation_signature`
    in spirit but unvalidated against the published implementation —
    install pertpy for production runs.
    """
    work = adata.copy() if layer is None else None
    target = work if work is not None else adata
    src_layer = layer

    if src_layer is None:
        # Standard scanpy log-normalize on a copy, leaving caller's X untouched.
        sc.pp.normalize_total(target, target_sum=1e4)
        sc.pp.log1p(target)
        src_layer = None

    matrix = (
        target.layers[src_layer]
        if src_layer is not None and src_layer in target.layers
        else target.X
    )
    dense = matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)

    # Restrict to genes not on the exclude prefix list (e.g. ribosomal,
    # mitochondrial). These contribute large variance unrelated to KO biology.
    exclude_mask = np.zeros(target.n_vars, dtype=bool)
    for prefix in cfg.exclude_gene_prefixes:
        exclude_mask |= target.var_names.str.startswith(prefix)
    kept = ~exclude_mask

    # PCA on kept genes; lightweight implementation that doesn't mutate adata.
    centered = dense[:, kept] - dense[:, kept].mean(axis=0, keepdims=True)
    n_comp = min(30, min(centered.shape) - 1)
    # SVD-based PCA: u * s gives PC scores.
    u, s, _vt = np.linalg.svd(centered, full_matrices=False)
    pcs = u[:, :n_comp] * s[:n_comp]

    # k-NN restricted to controls.
    control_mask = (target.obs[perturbation_key].astype(str) == cfg.control_label).values
    if control_mask.sum() < cfg.n_neighbors:
        raise ValueError(
            f"Need at least n_neighbors={cfg.n_neighbors} control cells; "
            f"found {int(control_mask.sum())}"
        )
    ctrl_pcs = pcs[control_mask]

    # Brute force k-NN distance from each cell to control cells. n_cells is
    # typically <500k in practice — fine without a kdtree for the fallback.
    sq_norm_cells = (pcs**2).sum(axis=1, keepdims=True)
    sq_norm_ctrl = (ctrl_pcs**2).sum(axis=1)
    cross = pcs @ ctrl_pcs.T
    d2 = sq_norm_cells + sq_norm_ctrl[None, :] - 2 * cross
    nn_idx = np.argpartition(d2, kth=cfg.n_neighbors, axis=1)[:, : cfg.n_neighbors]

    # Mean control expression for each cell's NN set, subtracted from cell expression.
    expr = dense
    sig = np.empty_like(expr)
    ctrl_expr = expr[control_mask]
    for i in range(expr.shape[0]):
        sig[i] = expr[i] - ctrl_expr[nn_idx[i]].mean(axis=0)

    adata.layers["perturbation_signature"] = sig


def _fallback_mixscape(
    adata: ad.AnnData,
    *,
    cfg: PerturbationAnalysisConfig,
    perturbation_key: str,
) -> None:
    """Two-component GMM on the L2 norm of each cell's perturbation signature.

    For each perturbation independently, fit a 2-component Gaussian mixture
    on the signature magnitudes. The component with the larger mean is the
    KO arm; cells assigned to that component are flagged ``KO``, the others
    ``NP``. Control cells are passed through with the control label.
    """
    from sklearn.mixture import GaussianMixture

    sig = adata.layers["perturbation_signature"]
    sig_dense = sig.toarray() if hasattr(sig, "toarray") else np.asarray(sig)
    magnitudes = np.linalg.norm(sig_dense, axis=1)

    classes = pd.Series(index=adata.obs_names, dtype=object, name="mixscape_class")
    perturbed = pd.Series(index=adata.obs_names, dtype=bool, name="mixscape_perturbed")
    perts = adata.obs[perturbation_key].astype(str)
    for pert in perts.unique():
        mask = perts == pert
        if pert == cfg.control_label:
            classes[mask] = cfg.control_label
            perturbed[mask] = False
            continue
        sub = magnitudes[mask.values]
        if len(sub) < 5:
            # Too few cells to fit; mark them all NP.
            classes[mask] = f"{pert} NP"
            perturbed[mask] = False
            continue
        gmm = GaussianMixture(n_components=2, random_state=0, reg_covar=1e-3)
        gmm.fit(sub.reshape(-1, 1))
        labels = gmm.predict(sub.reshape(-1, 1))
        ko_component = int(np.argmax(gmm.means_.ravel()))
        is_ko = labels == ko_component
        classes[mask] = np.where(is_ko, f"{pert} KO", f"{pert} NP")
        perturbed[mask] = is_ko

    adata.obs["mixscape_class"] = classes.astype("category")
    adata.obs["mixscape_class_global"] = (
        classes.str.extract(r"(KO|NP)$", expand=False).fillna(cfg.control_label).astype("category")
    )
    adata.obs["mixscape_perturbed"] = perturbed.astype(bool)
    logger.info(
        "Fallback Mixscape: %d KO cells, %d NP cells across %d perturbations",
        int(perturbed.sum()),
        int((~perturbed & (perts != cfg.control_label)).sum()),
        perts.nunique() - 1,
    )
