from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from perturbflow.adapters import (
    from_replogle_2022_anndata,
    guide_metadata_from_cellranger_features,
    read_cellranger_protospacer_calls,
)


def test_cellranger_protospacer_single_guide(tmp_path: Path) -> None:
    p = tmp_path / "protospacer_calls_per_cell.csv"
    p.write_text(
        "cell_barcode,num_features,feature_call,num_umis\n"
        "AAA-1,1,GENE_A_g1,42\n"
        "BBB-1,1,GENE_B_g1,15\n"
        "CCC-1,1,NT_g1,28\n"
    )
    df = read_cellranger_protospacer_calls(p)
    assert list(df.columns) == ["cell_barcode", "guide_id", "umi_count"]
    assert df.shape == (3, 3)
    assert df.set_index("cell_barcode").loc["AAA-1", "umi_count"] == 42


def test_cellranger_protospacer_multi_guide(tmp_path: Path) -> None:
    """Multi-feature cells must explode to one row per (cell, guide)."""
    p = tmp_path / "multi.csv"
    p.write_text(
        "cell_barcode,num_features,feature_call,num_umis\n"
        "AAA-1,2,GENE_A_g1|GENE_A_g2,23|17\n"
        "BBB-1,1,GENE_B_g1,30\n"
    )
    df = read_cellranger_protospacer_calls(p)
    assert df.shape == (3, 3)
    aaa = df[df["cell_barcode"] == "AAA-1"].sort_values("guide_id")
    assert list(aaa["guide_id"]) == ["GENE_A_g1", "GENE_A_g2"]
    assert list(aaa["umi_count"]) == [23, 17]


def test_cellranger_protospacer_drop_multifeature(tmp_path: Path) -> None:
    p = tmp_path / "multi.csv"
    p.write_text(
        "cell_barcode,num_features,feature_call,num_umis\n"
        "AAA-1,2,GENE_A_g1|GENE_A_g2,23|17\n"
        "BBB-1,1,GENE_B_g1,30\n"
    )
    df = read_cellranger_protospacer_calls(p, drop_multifeature=True)
    assert df["cell_barcode"].tolist() == ["BBB-1"]


def test_cellranger_protospacer_rejects_wrong_schema(tmp_path: Path) -> None:
    p = tmp_path / "wrong.csv"
    p.write_text("a,b,c\n1,2,3\n")
    with pytest.raises(ValueError, match="not a CellRanger"):
        read_cellranger_protospacer_calls(p)


def test_cellranger_protospacer_consistency_check(tmp_path: Path) -> None:
    """A row with 2 guides but only 1 UMI count should error, not silently drop."""
    p = tmp_path / "mismatch.csv"
    p.write_text(
        "cell_barcode,num_features,feature_call,num_umis\nAAA-1,2,GENE_A_g1|GENE_A_g2,23\n"
    )
    with pytest.raises(ValueError, match="Inconsistent feature_call"):
        read_cellranger_protospacer_calls(p)


def test_guide_metadata_from_cellranger_features(tmp_path: Path) -> None:
    p = tmp_path / "features.tsv"
    # CellRanger features.tsv schema for CRISPR: id, name, type, target_gene_id, target_gene_name
    p.write_text(
        "ENSG_dummy_A\tGENE_A\tGene Expression\t\t\n"
        "GENE_A_g1\tGENE_A_g1\tCRISPR Guide Capture\tENSG_GENE_A\tGENE_A\n"
        "GENE_A_g2\tGENE_A_g2\tCRISPR Guide Capture\tENSG_GENE_A\tGENE_A\n"
        "NT_g1\tNT_g1\tCRISPR Guide Capture\t\tnon-targeting\n"
    )
    meta = guide_metadata_from_cellranger_features(p)
    assert set(meta["guide_id"]) == {"GENE_A_g1", "GENE_A_g2", "NT_g1"}
    nt = meta[meta["guide_id"] == "NT_g1"]
    assert bool(nt["is_control"].iloc[0]) is True
    assert nt["target_gene"].iloc[0] == ""
    a1 = meta[meta["guide_id"] == "GENE_A_g1"]
    assert a1["target_gene"].iloc[0] == "GENE_A"
    assert bool(a1["is_control"].iloc[0]) is False


def test_guide_metadata_rejects_no_crispr_features(tmp_path: Path) -> None:
    p = tmp_path / "gex_only.tsv"
    p.write_text("ENSG_x\tGENE_X\tGene Expression\t\t\n")
    with pytest.raises(ValueError, match="no rows with feature_type"):
        guide_metadata_from_cellranger_features(p)


def test_from_replogle_2022_anndata() -> None:
    rng = np.random.default_rng(0)
    n = 30
    obs = pd.DataFrame(
        {"gene_target": ["GENE_A"] * 10 + ["GENE_B"] * 10 + ["non-targeting"] * 10},
        index=[f"cell_{i}" for i in range(n)],
    )
    X = rng.poisson(2, size=(n, 5)).astype(np.int32)
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(5)]))

    out = from_replogle_2022_anndata(adata)
    assert "perturbation" in out.obs
    assert (out.obs["perturbation"].astype(str).iloc[:10] == "GENE_A").all()
    assert (out.obs["perturbation"].astype(str).iloc[-10:] == "NT").all()
    assert "assignment_status" in out.obs
    assert (out.obs["assignment_status"].astype(str) == "assigned").all()
    assert bool(out.obs["is_control"].iloc[-1]) is True
    assert bool(out.obs["is_control"].iloc[0]) is False


def test_from_replogle_2022_missing_column() -> None:
    adata = ad.AnnData(
        X=np.zeros((3, 3), dtype=np.int32),
        obs=pd.DataFrame({"wrong": [1, 2, 3]}, index=["a", "b", "c"]),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )
    with pytest.raises(KeyError, match="gene_target"):
        from_replogle_2022_anndata(adata)
