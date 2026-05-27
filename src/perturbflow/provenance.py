"""Reproducibility provenance: capture the run context so results are auditable.

What we record
--------------
- Package version and a content hash of the active config.
- Git revision and dirty flag (so a manuscript pipeline run can be tied to a
  specific commit).
- Pinned versions of the major scientific dependencies.
- Seed, run name, output directory.

The provenance dict is written to ``perturbflow_run.log`` alongside the
outputs and embedded in the HTML report.
"""

from __future__ import annotations

import hashlib
import json
import logging
import platform
import subprocess
import sys
from dataclasses import asdict
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from perturbflow._version import __version__
from perturbflow.config import PerturbFlowConfig

logger = logging.getLogger(__name__)


PACKAGES_TO_RECORD = (
    "anndata",
    "scanpy",
    "numpy",
    "pandas",
    "scipy",
    "scikit-learn",
    "matplotlib",
    "pertpy",
    "pydeseq2",
    "decoupler",
)


def collect(config: PerturbFlowConfig, *, config_path: str | Path | None) -> dict[str, Any]:
    """Build the provenance record for a run."""
    cfg_dict = _config_to_dict(config)
    cfg_hash = hashlib.sha256(
        json.dumps(cfg_dict, sort_keys=True, default=str).encode()
    ).hexdigest()[:12]

    return {
        "perturbflow_version": __version__,
        "config_path": str(config_path) if config_path else None,
        "config_hash": cfg_hash,
        "config": cfg_dict,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "git": _git_info(),
        "packages": _package_versions(),
    }


def write(provenance: dict[str, Any], path: str | Path) -> Path:
    """Write the provenance dict as JSON to ``path``."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(provenance, indent=2, default=str, sort_keys=True))
    logger.info("Wrote provenance: %s", p)
    return p


# ----------------------------------------------------------------------


def _config_to_dict(cfg: PerturbFlowConfig) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for f, value in asdict(cfg).items():
        out[f] = value
    return out


def _git_info() -> dict[str, str | bool]:
    try:
        rev = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL, text=True
        ).strip()
        dirty = bool(
            subprocess.check_output(
                ["git", "status", "--porcelain"], stderr=subprocess.DEVNULL, text=True
            ).strip()
        )
        return {"revision": rev, "dirty": dirty}
    except (subprocess.CalledProcessError, FileNotFoundError):
        return {"revision": "unknown", "dirty": False}


def _package_versions() -> dict[str, str]:
    out: dict[str, str] = {}
    for name in PACKAGES_TO_RECORD:
        try:
            out[name] = version(name)
        except PackageNotFoundError:
            out[name] = "not installed"
    return out
