"""PerturbFlow: opinionated Perturb-seq analysis pipeline.

Public surface
--------------
The top-level package re-exports the small set of functions that compose the
pipeline. Use these when building a notebook workflow; the CLI (``perturbflow
run``) wires them together for batch jobs.

>>> import perturbflow as pf
>>> adata = pf.read_10x_h5("filtered_feature_bc_matrix.h5")
>>> guides = pf.read_guide_calls("guide_calls.csv")
>>> adata = pf.assign_guides(adata, guides, guide_metadata=meta)
>>> adata = pf.run_mixscape(adata)
>>> de = pf.run_pseudobulk_de(adata)
"""

from __future__ import annotations

from perturbflow._version import __version__
from perturbflow.adapters import (
    from_replogle_2022_anndata,
    guide_metadata_from_cellranger_features,
    read_cellranger_protospacer_calls,
)
from perturbflow.de import run_pseudobulk_de
from perturbflow.downstream import compute_cell_state_effects, score_pathways
from perturbflow.guide_assignment import assign_guides
from perturbflow.io import (
    read_10x_h5,
    read_10x_mtx,
    read_guide_calls,
    read_guide_metadata,
    read_h5ad,
)
from perturbflow.perturbation_analysis import (
    compute_perturbation_signature,
    run_mixscape,
)
from perturbflow.qc import (
    per_cell_qc,
    per_guide_qc,
    per_perturbation_qc,
)
from perturbflow.report import write_html_report
from perturbflow.validation import (
    assert_raw_counts,
    multi_guide_concordance,
    warn_if_pseudoreplicates,
)

__all__ = [
    "__version__",
    "assert_raw_counts",
    "assign_guides",
    "compute_cell_state_effects",
    "compute_perturbation_signature",
    "from_replogle_2022_anndata",
    "guide_metadata_from_cellranger_features",
    "multi_guide_concordance",
    "per_cell_qc",
    "per_guide_qc",
    "per_perturbation_qc",
    "read_10x_h5",
    "read_10x_mtx",
    "read_cellranger_protospacer_calls",
    "read_guide_calls",
    "read_guide_metadata",
    "read_h5ad",
    "run_mixscape",
    "run_pseudobulk_de",
    "score_pathways",
    "warn_if_pseudoreplicates",
    "write_html_report",
]
