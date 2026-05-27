"""Shared pytest fixtures.

We build the synthetic Perturb-seq dataset once per session and share it
across tests. Mutating tests take a ``.copy()``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd
import pytest

from tests.fixtures.make_synthetic import build_synthetic


@pytest.fixture(scope="session")
def synthetic_bundle() -> dict[str, Any]:
    """Build the (adata, guide_calls, guide_metadata) bundle once per session."""
    return build_synthetic(seed=0)


@pytest.fixture(scope="session")
def synthetic_adata(synthetic_bundle: dict[str, Any]) -> ad.AnnData:
    return synthetic_bundle["adata"]


@pytest.fixture(scope="session")
def synthetic_guide_calls(synthetic_bundle: dict[str, Any]) -> pd.DataFrame:
    return synthetic_bundle["guide_calls"]


@pytest.fixture(scope="session")
def synthetic_guide_metadata(synthetic_bundle: dict[str, Any]) -> pd.DataFrame:
    return synthetic_bundle["guide_metadata"]


@pytest.fixture
def tmp_run_dir(tmp_path: Path) -> Path:
    out = tmp_path / "perturbflow_run"
    out.mkdir()
    return out


def _have(mod: str) -> bool:
    from importlib.util import find_spec

    return find_spec(mod) is not None


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Skip tests marked requires_X when X isn't installed."""
    skip_map = {
        "requires_pertpy": pytest.mark.skip(reason="pertpy not installed"),
        "requires_pydeseq2": pytest.mark.skip(reason="pydeseq2 not installed"),
        "requires_decoupler": pytest.mark.skip(reason="decoupler not installed"),
    }
    have = {
        "requires_pertpy": _have("pertpy"),
        "requires_pydeseq2": _have("pydeseq2"),
        "requires_decoupler": _have("decoupler"),
    }
    for item in items:
        for mark_name, present in have.items():
            if mark_name in item.keywords and not present:
                item.add_marker(skip_map[mark_name])
