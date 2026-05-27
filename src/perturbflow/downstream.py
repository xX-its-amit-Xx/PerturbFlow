"""Downstream interpretation: pathway scoring and cell-state effect maps.

Two questions this module answers:

1. **What pathways are perturbed?** Given DE statistics per perturbation,
   score gene-set collections (MSigDB Hallmarks, Reactome, PROGENy) via
   :mod:`decoupler-py` so each perturbation gets a per-pathway activity
   estimate plus a significance.
2. **Where in cell-state space does the perturbation push cells?** Project
   perturbed and control cells onto UMAP coordinates and quantify the
   centroid shift per perturbation. This is the figure that goes in slide 1
   of the manuscript — "here's the cell-state landscape and here's how
   each KO moves cells through it".
"""

from __future__ import annotations

import logging
from importlib.util import find_spec

import anndata as ad
import numpy as np
import pandas as pd

from perturbflow.config import DownstreamConfig

logger = logging.getLogger(__name__)


def score_pathways(
    de_results: dict[str, pd.DataFrame],
    *,
    config: DownstreamConfig | None = None,
    organism: str = "human",
) -> pd.DataFrame:
    """Score pathway activity per perturbation from DE statistics.

    Decoupler's ``run_ulm`` (univariate linear model) is the default — it's
    robust, fast, and works on a single statistic per gene per
    perturbation. We feed it the ``stat`` column (Wald statistic from
    DESeq2, or t-statistic from the fallback) because using log2FC alone
    underweights highly-variable genes; the test statistic is the right
    quantity for a t-distribution-based decoupler test.

    Returns a long-format DataFrame:
    ``perturbation | pathway | score | pvalue``.

    When decoupler is unavailable, returns an empty DataFrame and logs a
    warning. The pipeline does not fail because pathway scoring is an
    interpretive step, not a structural prerequisite for the report.
    """
    cfg = config or DownstreamConfig()
    if not cfg.enable_pathway_scoring:
        return pd.DataFrame(columns=["perturbation", "pathway", "score", "pvalue"])
    if not de_results:
        logger.warning("score_pathways: no DE results to score; returning empty frame")
        return pd.DataFrame(columns=["perturbation", "pathway", "score", "pvalue"])

    if find_spec("decoupler") is None:
        logger.warning(
            "decoupler is not installed; skipping pathway scoring. "
            "Install with `pip install perturbflow[pathways]` to enable it."
        )
        return pd.DataFrame(columns=["perturbation", "pathway", "score", "pvalue"])

    net = _load_pathway_network(cfg.pathway_net, organism=organism)
    if net.empty:
        logger.warning("Pathway network %r returned no entries; skipping", cfg.pathway_net)
        return pd.DataFrame(columns=["perturbation", "pathway", "score", "pvalue"])

    import decoupler as dc

    # Build (perturbation × gene) matrix of test statistics.
    wide = (
        pd.concat(
            {pert: df.set_index("gene")["stat"] for pert, df in de_results.items()},
            axis=1,
        )
        .fillna(0.0)
        .T
    )
    wide.index.name = "perturbation"

    method = cfg.pathway_method
    runner = getattr(dc, f"run_{method}", None)
    if runner is None:
        raise ValueError(
            f"Unknown decoupler method {method!r}; expected one of run_ulm/run_mlm/run_wsum"
        )
    estimate, pvals = runner(
        mat=wide,
        net=net,
        source="source",
        target="target",
        weight="weight",
        verbose=False,
    )
    long = (
        estimate.stack()
        .rename("score")
        .reset_index()
        .rename(columns={"level_0": "perturbation", "level_1": "pathway"})
    )
    long["pvalue"] = pvals.stack().reset_index(drop=True).values
    return long.sort_values(["perturbation", "pvalue"]).reset_index(drop=True)


def compute_cell_state_effects(
    adata: ad.AnnData,
    *,
    perturbation_key: str = "perturbation",
    control_label: str = "NT",
    embedding: str = "X_umap",
) -> pd.DataFrame:
    """Quantify how each perturbation shifts cell-state in an embedding.

    Computes per-perturbation:

    - ``centroid_shift`` — Euclidean distance between the perturbation's
      mean embedding coordinate and the control mean.
    - ``dispersion_ratio`` — perturbation's mean intra-group distance /
      control's mean intra-group distance. Values above 1 indicate the
      perturbation broadened cell-state heterogeneity.
    - ``n_cells`` — number of cells contributing.

    Requires an embedding to already be in ``adata.obsm[embedding]``. If
    it's missing we compute a minimal UMAP via Scanpy so the function does
    something useful out of the box (the pipeline's CLI does this for you).
    """
    if perturbation_key not in adata.obs.columns:
        raise KeyError(f"adata.obs[{perturbation_key!r}] not found")
    if embedding not in adata.obsm:
        raise KeyError(
            f"adata.obsm[{embedding!r}] missing — run sc.pp.neighbors + sc.tl.umap "
            "before compute_cell_state_effects, or use the CLI which does it for you."
        )

    coords = np.asarray(adata.obsm[embedding])
    perts = adata.obs[perturbation_key].astype(str).values
    if control_label not in set(perts):
        raise ValueError(f"control_label {control_label!r} not in {perturbation_key!r}")
    ctrl_mask = perts == control_label
    ctrl_centroid = coords[ctrl_mask].mean(axis=0)
    ctrl_disp = float(np.linalg.norm(coords[ctrl_mask] - ctrl_centroid, axis=1).mean())
    ctrl_disp = max(ctrl_disp, 1e-9)

    rows = []
    for pert in sorted(set(perts)):
        if pert == control_label:
            continue
        mask = perts == pert
        sub = coords[mask]
        centroid = sub.mean(axis=0)
        shift = float(np.linalg.norm(centroid - ctrl_centroid))
        disp = float(np.linalg.norm(sub - centroid, axis=1).mean())
        rows.append(
            {
                "perturbation": pert,
                "n_cells": int(mask.sum()),
                "centroid_shift": shift,
                "dispersion_ratio": disp / ctrl_disp,
            }
        )
    return pd.DataFrame(rows).sort_values("centroid_shift", ascending=False).reset_index(drop=True)


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------


def _load_pathway_network(name: str, *, organism: str) -> pd.DataFrame:
    """Resolve a pathway-network name into a decoupler-compatible long DataFrame.

    Supports the three commonly-used collections plus an escape hatch for
    a user-supplied TSV (path with a ``/`` or ``\\`` or ending in ``.tsv``/
    ``.csv``).
    """
    if find_spec("decoupler") is None:
        return pd.DataFrame()
    import decoupler as dc

    key = name.strip().lower()
    if key in {"hallmarks", "hallmark", "msigdb_hallmarks"}:
        net = dc.get_resource("MSigDB", organism=organism)
        if isinstance(net, pd.DataFrame) and "collection" in net.columns:
            net = net[net["collection"] == "hallmark"]
        if "geneset" in net.columns:
            net = net.rename(columns={"geneset": "source", "genesymbol": "target"})
        net["weight"] = 1.0
        return net[["source", "target", "weight"]].drop_duplicates()
    if key == "reactome":
        net = dc.get_resource("MSigDB", organism=organism)
        if isinstance(net, pd.DataFrame) and "collection" in net.columns:
            net = net[net["collection"].str.contains("reactome", case=False, na=False)]
        if "geneset" in net.columns:
            net = net.rename(columns={"geneset": "source", "genesymbol": "target"})
        net["weight"] = 1.0
        return net[["source", "target", "weight"]].drop_duplicates()
    if key == "progeny":
        return dc.get_progeny(top=500, organism=organism)
    # Treat anything else as a path.
    from pathlib import Path

    p = Path(name)
    if not p.exists():
        raise ValueError(f"Unknown pathway_net {name!r}: not a recognized name and not a file path")
    sep = "\t" if p.suffix.lower() in {".tsv", ".txt"} else ","
    df = pd.read_csv(p, sep=sep)
    needed = {"source", "target", "weight"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"Custom pathway network {p} missing columns {missing}")
    return df[["source", "target", "weight"]]
