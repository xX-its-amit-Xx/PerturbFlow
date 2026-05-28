"""Tests for the per-stage runners that the Snakemake DAG calls."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd
import yaml

from perturbflow.config import load_config
from perturbflow.stages import (
    STAGE_REGISTRY,
    stage_assign,
    stage_de,
    stage_embedding,
    stage_load,
    stage_mixscape,
    stage_pathways,
    stage_qc,
)


def _write_config(tmp_path: Path, synthetic_bundle: dict) -> Path:
    """Write inputs to disk and produce a valid config pointing at them."""
    h5ad_path = tmp_path / "raw.h5ad"
    synthetic_bundle["adata"].write_h5ad(h5ad_path)
    g_path = tmp_path / "guides.csv"
    synthetic_bundle["guide_calls"].to_csv(g_path, index=False)
    m_path = tmp_path / "meta.csv"
    synthetic_bundle["guide_metadata"].to_csv(m_path, index=False)
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(
        yaml.safe_dump(
            {
                "run": {"name": "stages", "outdir": str(tmp_path / "out"), "seed": 0},
                "input": {
                    "h5ad": str(h5ad_path),
                    "guide_calls": str(g_path),
                    "guide_metadata": str(m_path),
                    "n_pseudo_replicates": 3,
                },
                "qc": {"min_genes_per_cell": 5},
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
                "report": {"enable": False},
            }
        )
    )
    return cfg_path


def test_stage_registry_complete() -> None:
    assert set(STAGE_REGISTRY) == {
        "load",
        "assign",
        "mixscape",
        "embedding",
        "qc",
        "de",
        "pathways",
    }
    for fn in STAGE_REGISTRY.values():
        assert callable(fn)


def test_stage_load_writes_h5ad(tmp_path: Path, synthetic_bundle: dict) -> None:
    cfg_path = _write_config(tmp_path, synthetic_bundle)
    cfg = load_config(cfg_path)
    out = tmp_path / "01_raw.h5ad"
    stage_load(cfg, out)
    assert out.exists()
    a = ad.read_h5ad(out)
    assert a.n_obs == synthetic_bundle["adata"].n_obs


def test_stage_assign_writes_h5ad_and_per_guide(tmp_path: Path, synthetic_bundle: dict) -> None:
    cfg_path = _write_config(tmp_path, synthetic_bundle)
    cfg = load_config(cfg_path)
    raw = tmp_path / "01_raw.h5ad"
    stage_load(cfg, raw)
    assigned = tmp_path / "02_assigned.h5ad"
    per_guide = tmp_path / "out" / "qc" / "per_guide.csv"
    stage_assign(cfg, raw, assigned, per_guide)
    assert assigned.exists()
    assert per_guide.exists()
    a = ad.read_h5ad(assigned)
    assert "perturbation" in a.obs.columns
    pg = pd.read_csv(per_guide)
    assert {"guide_id", "n_cells", "is_control"}.issubset(pg.columns)


def test_stage_full_dag_in_order(tmp_path: Path, synthetic_bundle: dict) -> None:
    """Run every stage in order and check the final outputs exist."""
    cfg_path = _write_config(tmp_path, synthetic_bundle)
    cfg = load_config(cfg_path)
    outdir = tmp_path / "out"
    stages_dir = outdir / "stages"
    qc_dir = outdir / "qc"
    de_dir = outdir / "de"

    raw = stages_dir / "01_raw.h5ad"
    stage_load(cfg, raw)
    assigned = stages_dir / "02_assigned.h5ad"
    stage_assign(cfg, raw, assigned, qc_dir / "per_guide.csv")
    mix = stages_dir / "03_mixscape.h5ad"
    stage_mixscape(cfg, assigned, mix)
    emb = stages_dir / "04_embedded.h5ad"
    stage_embedding(cfg, mix, emb, qc_dir / "cell_state_effects.csv")
    stage_qc(cfg, emb, qc_dir / "per_cell.csv", qc_dir / "per_perturbation.csv")
    stage_de(cfg, emb, de_dir)
    stage_pathways(cfg, de_dir, de_dir / "pathway_scores.csv")

    for p in [
        raw,
        assigned,
        mix,
        emb,
        qc_dir / "per_guide.csv",
        qc_dir / "per_cell.csv",
        qc_dir / "per_perturbation.csv",
        qc_dir / "cell_state_effects.csv",
        de_dir / "pathway_scores.csv",
    ]:
        assert p.exists(), f"Missing stage output: {p}"

    # DE produced at least one perturbation table.
    de_csvs = [c for c in de_dir.glob("*.csv") if c.stem != "pathway_scores"]
    assert len(de_csvs) >= 1
