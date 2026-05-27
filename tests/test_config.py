from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from perturbflow.config import ConfigError, from_dict, load_config


def _write_yaml(path: Path, content: dict) -> Path:
    path.write_text(yaml.safe_dump(content))
    return path


def test_load_default_workflow_config(tmp_path: Path) -> None:
    cfg_path = _write_yaml(
        tmp_path / "c.yaml",
        {
            "run": {"name": "x", "outdir": "out", "seed": 7},
            "input": {"matrix_h5": "a.h5", "guide_calls": "g.csv", "guide_metadata": "m.csv"},
        },
    )
    cfg = load_config(cfg_path)
    assert cfg.run.name == "x"
    assert cfg.run.seed == 7
    assert cfg.input.matrix_h5 == "a.h5"
    # Defaults filled in for unspecified sections.
    assert cfg.qc.min_genes_per_cell == 500
    assert cfg.perturbation_analysis.control_label == "NT"


def test_unknown_section_rejected(tmp_path: Path) -> None:
    cfg_path = _write_yaml(tmp_path / "c.yaml", {"run": {}, "wat": {"x": 1}})
    with pytest.raises(ConfigError, match="wat"):
        load_config(cfg_path)


def test_unknown_key_rejected(tmp_path: Path) -> None:
    cfg_path = _write_yaml(tmp_path / "c.yaml", {"run": {"weird_key": 1}})
    with pytest.raises(ConfigError, match="weird_key"):
        load_config(cfg_path)


def test_matrix_source_validation() -> None:
    cfg = from_dict(
        {
            "input": {
                "matrix_h5": "a.h5",
                "matrix_mtx_dir": "b/",
                "guide_calls": "g.csv",
                "guide_metadata": "m.csv",
            }
        }
    )
    with pytest.raises(ConfigError):
        cfg.input.matrix_source()


def test_matrix_source_picks_kind() -> None:
    cfg = from_dict(
        {"input": {"h5ad": "x.h5ad", "guide_calls": "g.csv", "guide_metadata": "m.csv"}}
    )
    kind, path = cfg.input.matrix_source()
    assert kind == "h5ad"
    assert path == "x.h5ad"


def test_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(ConfigError, match="not found"):
        load_config(tmp_path / "nope.yaml")
