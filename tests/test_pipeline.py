"""End-to-end pipeline integration test on the synthetic fixture."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from perturbflow.config import load_config
from perturbflow.pipeline import run


def _write_inputs(tmp_path: Path, synthetic_bundle: dict) -> Path:
    """Materialize the synthetic bundle to disk and write a YAML config."""
    adata = synthetic_bundle["adata"]
    h5ad_path = tmp_path / "raw.h5ad"
    adata.write_h5ad(h5ad_path)

    guide_path = tmp_path / "guide_calls.csv"
    synthetic_bundle["guide_calls"].to_csv(guide_path, index=False)
    meta_path = tmp_path / "guide_metadata.csv"
    synthetic_bundle["guide_metadata"].to_csv(meta_path, index=False)

    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(
        yaml.safe_dump(
            {
                "run": {
                    "name": "integration",
                    "outdir": str(tmp_path / "out"),
                    "seed": 0,
                    "log_level": "WARNING",
                },
                "input": {
                    "h5ad": str(h5ad_path),
                    "guide_calls": str(guide_path),
                    "guide_metadata": str(meta_path),
                    "n_pseudo_replicates": 3,
                },
                "qc": {"mito_prefix": "MT-", "min_genes_per_cell": 5},
                "guide_assignment": {"drop_unassigned": True},
                "perturbation_analysis": {
                    "enable_mixscape": True,
                    "control_label": "NT",
                    "n_neighbors": 15,
                },
                "de": {
                    "enable": True,
                    "min_replicates_per_group": 2,
                    "min_cells_per_replicate": 5,
                    "lfc_threshold": 0.5,
                    "padj_threshold": 0.1,
                    "use_mixscape_filter": True,
                },
                "downstream": {"enable_pathway_scoring": False},
                "report": {"enable": True, "bundle_html": True},
            }
        )
    )
    return cfg_path


@pytest.mark.integration
def test_end_to_end_pipeline(tmp_path: Path, synthetic_bundle: dict) -> None:
    cfg_path = _write_inputs(tmp_path, synthetic_bundle)
    cfg = load_config(cfg_path)
    artifacts = run(cfg, config_path=cfg_path)

    outdir = Path(cfg.run.outdir)
    assert (outdir / "report.html").exists()
    assert (outdir / "perturbflow.h5ad").exists()
    assert (outdir / "qc" / "per_cell.csv").exists()
    assert (outdir / "qc" / "per_guide.csv").exists()
    assert (outdir / "qc" / "per_perturbation.csv").exists()
    assert (outdir / "provenance.json").exists()

    # DE should have one CSV per KO perturbation.
    de_csvs = list((outdir / "de").glob("*.csv"))
    de_names = {p.stem for p in de_csvs}
    assert {"KOA", "KOB", "KOC", "KOD"}.issubset(de_names)

    # Pipeline returned artifacts should be consistent.
    assert artifacts.report_path is not None
    assert artifacts.report_path.exists()
    assert len(artifacts.de_results) >= 4
