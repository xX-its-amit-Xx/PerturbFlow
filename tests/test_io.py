from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from perturbflow.io import (
    read_guide_calls,
    read_guide_metadata,
    read_h5ad,
)


def test_read_guide_calls_validates_columns(tmp_path: Path) -> None:
    bad = tmp_path / "bad.csv"
    bad.write_text("a,b\n1,2\n")
    with pytest.raises(ValueError, match="missing required columns"):
        read_guide_calls(bad)


def test_read_guide_calls_round_trip(tmp_path: Path, synthetic_guide_calls: pd.DataFrame) -> None:
    p = tmp_path / "g.csv"
    synthetic_guide_calls.to_csv(p, index=False)
    df = read_guide_calls(p)
    assert list(df.columns) == ["cell_barcode", "guide_id", "umi_count"]
    assert df["umi_count"].dtype.kind in {"i", "u"}
    assert df.equals(
        synthetic_guide_calls.astype({"cell_barcode": str, "guide_id": str, "umi_count": int})
    )


def test_read_guide_calls_tsv(tmp_path: Path, synthetic_guide_calls: pd.DataFrame) -> None:
    p = tmp_path / "g.tsv"
    synthetic_guide_calls.to_csv(p, sep="\t", index=False)
    df = read_guide_calls(p)
    assert len(df) == len(synthetic_guide_calls)


def test_read_guide_metadata_requires_control(tmp_path: Path) -> None:
    p = tmp_path / "m.csv"
    pd.DataFrame({"guide_id": ["g1"], "target_gene": ["X"], "is_control": [False]}).to_csv(
        p, index=False
    )
    with pytest.raises(ValueError, match="no control guides"):
        read_guide_metadata(p)


def test_read_guide_metadata_rejects_duplicates(tmp_path: Path) -> None:
    p = tmp_path / "m.csv"
    pd.DataFrame(
        {
            "guide_id": ["g1", "g1", "nt1"],
            "target_gene": ["X", "X", ""],
            "is_control": [False, False, True],
        }
    ).to_csv(p, index=False)
    with pytest.raises(ValueError, match="duplicate"):
        read_guide_metadata(p)


def test_read_h5ad_round_trip(tmp_path: Path, synthetic_adata) -> None:
    p = tmp_path / "x.h5ad"
    synthetic_adata.write_h5ad(p)
    a = read_h5ad(p)
    assert a.shape == synthetic_adata.shape


def test_read_h5ad_missing() -> None:
    with pytest.raises(FileNotFoundError):
        read_h5ad("/no/such/file.h5ad")
