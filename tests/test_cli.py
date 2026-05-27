from __future__ import annotations

from pathlib import Path

import yaml
from click.testing import CliRunner

from perturbflow.cli import main


def test_version_command() -> None:
    result = CliRunner().invoke(main, ["version"])
    assert result.exit_code == 0
    assert result.output.strip()


def test_validate_config_passes(tmp_path: Path) -> None:
    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        yaml.safe_dump(
            {
                "input": {
                    "h5ad": "x.h5ad",
                    "guide_calls": "g.csv",
                    "guide_metadata": "m.csv",
                }
            }
        )
    )
    result = CliRunner().invoke(main, ["validate-config", "--config", str(cfg)])
    assert result.exit_code == 0


def test_validate_config_rejects_unknown_key(tmp_path: Path) -> None:
    cfg = tmp_path / "c.yaml"
    cfg.write_text(yaml.safe_dump({"run": {"unknown_key": 1}}))
    result = CliRunner().invoke(main, ["validate-config", "--config", str(cfg)])
    assert result.exit_code != 0


def test_help() -> None:
    result = CliRunner().invoke(main, ["--help"])
    assert result.exit_code == 0
    assert "Perturb-seq" in result.output
