"""Strongly-typed configuration loader for PerturbFlow.

The shape mirrors ``workflow/config.yaml`` exactly. Each subsection is its own
``dataclass`` so that downstream modules can accept the narrow slice they need
(``GuideAssignmentConfig`` vs the full ``PerturbFlowConfig``) without picking
keys out of a dict.

The loader is intentionally strict: unknown keys raise rather than silently
ignore, so a typo in a long-running batch config fails fast.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass
from pathlib import Path
from typing import Any

import yaml


class ConfigError(ValueError):
    """Raised when the YAML config cannot be validated against the schema."""


@dataclass(frozen=True)
class RunConfig:
    name: str = "perturbflow_run"
    outdir: str = "perturbflow"
    seed: int = 0
    log_level: str = "INFO"


@dataclass(frozen=True)
class InputConfig:
    matrix_h5: str | None = None
    matrix_mtx_dir: str | None = None
    h5ad: str | None = None
    guide_calls: str = ""
    guide_metadata: str = ""
    sample_col: str | None = None
    n_pseudo_replicates: int = 3

    def matrix_source(self) -> tuple[str, str]:
        """Return (kind, path) for whichever input matrix is configured.

        Exactly one of ``matrix_h5``, ``matrix_mtx_dir``, ``h5ad`` must be set.
        """
        candidates = {
            "h5": self.matrix_h5,
            "mtx": self.matrix_mtx_dir,
            "h5ad": self.h5ad,
        }
        provided = {k: v for k, v in candidates.items() if v}
        if len(provided) != 1:
            raise ConfigError(
                "Exactly one of input.matrix_h5, input.matrix_mtx_dir, input.h5ad "
                f"must be set; got {list(provided)}"
            )
        kind, path = next(iter(provided.items()))
        assert path is not None
        return kind, path


@dataclass(frozen=True)
class QCConfig:
    min_genes_per_cell: int = 500
    max_pct_mito: float = 20.0
    min_cells_per_gene: int = 10
    mito_prefix: str = "MT-"


@dataclass(frozen=True)
class GuideAssignmentConfig:
    min_guide_umi: int = 5
    dominance_ratio: float = 2.0
    max_guides: int = 1
    drop_unassigned: bool = True


@dataclass(frozen=True)
class PerturbationAnalysisConfig:
    enable_mixscape: bool = True
    control_label: str = "NT"
    n_neighbors: int = 20
    mixscape_pval_cutoff: float = 5.0e-2
    exclude_gene_prefixes: tuple[str, ...] = ("RPS", "RPL", "MT-")


@dataclass(frozen=True)
class DEConfig:
    enable: bool = True
    min_replicates_per_group: int = 2
    min_cells_per_replicate: int = 10
    lfc_threshold: float = 1.0
    padj_threshold: float = 0.05
    use_mixscape_filter: bool = True


@dataclass(frozen=True)
class DownstreamConfig:
    enable_pathway_scoring: bool = True
    pathway_net: str = "hallmarks"
    pathway_method: str = "ulm"
    emit_umap_overlays: bool = True


@dataclass(frozen=True)
class ReportConfig:
    enable: bool = True
    bundle_html: bool = True


@dataclass(frozen=True)
class PerturbFlowConfig:
    run: RunConfig = field(default_factory=RunConfig)
    input: InputConfig = field(default_factory=InputConfig)
    qc: QCConfig = field(default_factory=QCConfig)
    guide_assignment: GuideAssignmentConfig = field(default_factory=GuideAssignmentConfig)
    perturbation_analysis: PerturbationAnalysisConfig = field(
        default_factory=PerturbationAnalysisConfig
    )
    de: DEConfig = field(default_factory=DEConfig)
    downstream: DownstreamConfig = field(default_factory=DownstreamConfig)
    report: ReportConfig = field(default_factory=ReportConfig)


def _build(cls: type, data: dict[str, Any] | None, *, path: str) -> Any:
    """Construct a dataclass instance, validating that no unknown keys slipped in."""
    if data is None:
        return cls()
    if not isinstance(data, dict):
        raise ConfigError(f"{path}: expected a mapping, got {type(data).__name__}")
    allowed = {f.name for f in fields(cls)}
    unknown = set(data) - allowed
    if unknown:
        raise ConfigError(f"{path}: unknown key(s) {sorted(unknown)} (allowed: {sorted(allowed)})")
    kwargs: dict[str, Any] = {}
    for f in fields(cls):
        if f.name not in data:
            continue
        value = data[f.name]
        if is_dataclass(f.type):
            kwargs[f.name] = _build(f.type, value, path=f"{path}.{f.name}")
        elif (
            f.name == "exclude_gene_prefixes" and value is not None and not isinstance(value, tuple)
        ):
            kwargs[f.name] = tuple(value)
        else:
            kwargs[f.name] = value
    return cls(**kwargs)


def load_config(path: str | Path) -> PerturbFlowConfig:
    """Load and validate a YAML config file."""
    p = Path(path)
    if not p.exists():
        raise ConfigError(f"Config file not found: {p}")
    with p.open() as fh:
        raw = yaml.safe_load(fh) or {}
    if not isinstance(raw, dict):
        raise ConfigError(f"Top-level config must be a mapping, got {type(raw).__name__}")
    allowed = {f.name for f in fields(PerturbFlowConfig)}
    unknown = set(raw) - allowed
    if unknown:
        raise ConfigError(f"Unknown top-level section(s) {sorted(unknown)}")
    sub: dict[str, Any] = {}
    for f in fields(PerturbFlowConfig):
        sub[f.name] = _build(f.type, raw.get(f.name), path=f.name)
    return PerturbFlowConfig(**sub)


def from_dict(data: dict[str, Any]) -> PerturbFlowConfig:
    """Build a config from an in-memory dict (used in tests and notebooks)."""
    allowed = {f.name for f in fields(PerturbFlowConfig)}
    unknown = set(data) - allowed
    if unknown:
        raise ConfigError(f"Unknown top-level section(s) {sorted(unknown)}")
    sub: dict[str, Any] = {}
    for f in fields(PerturbFlowConfig):
        sub[f.name] = _build(f.type, data.get(f.name), path=f.name)
    return PerturbFlowConfig(**sub)
