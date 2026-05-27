from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from perturbflow.report import write_html_report


def _fake_per_cell() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cell_barcode": [f"c{i}" for i in range(8)],
            "total_counts": [1000] * 8,
            "n_genes": [600] * 8,
            "pct_mito": [5.0] * 8,
            "guide": ["g1"] * 4 + ["g_nt"] * 4,
            "perturbation": ["KOA"] * 4 + ["NT"] * 4,
            "assignment_status": ["assigned"] * 8,
            "passes_qc": [True] * 8,
        }
    )


def _fake_per_perturbation() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "perturbation": ["KOA"],
            "n_cells": [4],
            "escape_fraction": [0.1],
            "on_target_log2fc": [-2.3],
        }
    )


def _fake_de() -> dict[str, pd.DataFrame]:
    return {
        "KOA": pd.DataFrame(
            {
                "gene": ["KOA", "G1", "G2"],
                "log2FoldChange": [-2.3, 1.5, -0.2],
                "padj": [1e-10, 1e-5, 0.5],
                "stat": [-15.0, 6.0, -0.4],
                "baseMean": [100.0, 50.0, 20.0],
                "significant": [True, True, False],
            }
        )
    }


def test_write_html_report_creates_file(tmp_path: Path) -> None:
    out = tmp_path / "report.html"
    written = write_html_report(
        out_path=out,
        run_name="unit_test",
        seed=0,
        version="0.0.0",
        git_rev="abc1234",
        per_cell=_fake_per_cell(),
        per_perturbation=_fake_per_perturbation(),
        de_results=_fake_de(),
        pathway_scores=pd.DataFrame(columns=["perturbation", "pathway", "score", "pvalue"]),
        provenance={"placeholder": True},
    )
    assert written == out
    html = out.read_text(encoding="utf-8")
    assert "PerturbFlow" in html
    assert "KOA" in html
    assert "Run summary" in html
    # Inline base64 plot present.
    assert "data:image/png;base64," in html


def test_write_html_report_with_umap(tmp_path: Path) -> None:
    out = tmp_path / "report.html"
    coords = np.random.default_rng(0).normal(size=(8, 2))
    labels = pd.Series(["KOA"] * 4 + ["NT"] * 4, name="perturbation")
    write_html_report(
        out_path=out,
        run_name="umap_test",
        seed=0,
        version="0.0.0",
        git_rev="abc1234",
        per_cell=_fake_per_cell(),
        per_perturbation=_fake_per_perturbation(),
        de_results=_fake_de(),
        pathway_scores=pd.DataFrame(columns=["perturbation", "pathway", "score", "pvalue"]),
        umap_coords=coords,
        umap_labels=labels,
    )
    assert "Cell-state landscape" in out.read_text(encoding="utf-8")
