# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Doublet detection** (`qc.detect_doublets`) wrapping
  `scanpy.pp.scrublet`. Adds `obs['doublet_score']` and
  `obs['predicted_doublet']`. Wired into the pipeline via
  `qc.detect_doublets` / `qc.drop_doublets` config flags (both off by
  default). Falls back to an explicit threshold (0.25) when
  scikit-image isn't available, so the hook works without the
  `perturbflow[doublets]` extra installed.
- **`perturbflow[doublets]` install extra** (`scrublet>=0.2.3` +
  scikit-image).
- **"Project status and roadmap" section** in README explicitly listing
  what's shipped, what's in progress (PyPI publish, Zenodo DOI), and
  what's TBD (Mixscale wrapper, batch correction integration,
  multi-factor DE design, Snakemake cluster profiles).

### Added

- **CellRanger adapter** (`adapters.read_cellranger_protospacer_calls`)
  for the 10x `crispr_analysis/protospacer_calls_per_cell.csv` schema.
  Pipe-separated multi-feature cells are expanded to one row per
  (cell, guide) so the downstream multi-guide classifier sees them.
- **CellRanger features parser** (`adapters.guide_metadata_from_cellranger_features`)
  for `features.tsv` with `CRISPR Guide Capture` feature type. Configurable
  non-targeting pattern.
- **Replogle 2022 adapter** (`adapters.from_replogle_2022_anndata`) that
  normalizes the figshare release's `gene_target` schema to PerturbFlow's
  canonical obs columns.
- **Raw-counts validation** (`validation.assert_raw_counts`). DE now
  fails loud if `adata.X` is log-normalized, instead of silently producing
  a wrong volcano.
- **Pseudo-replicate warning** (`validation.warn_if_pseudoreplicates`)
  emitted at the DE call site, not just buried in `docs/methodology.md`.
- **Multi-guide concordance** (`validation.multi_guide_concordance`):
  Pearson r of log2FC vectors across guides targeting the same gene.
- **Per-stage Snakemake DAG.** The previous Snakefile was a single-rule
  wrapper; the new one has eight independently re-runnable stages with
  cached intermediates and per-stage timing for `snakemake --report`.
- **`perturbflow stage <name>`** CLI subcommand for the new DAG.
- **Pertpy-parity test** that verifies our Mixscape wrapper produces the
  same classifications as directly-called pertpy (≥98% agreement).
- **Ground-truth recovery integration test** that runs the full
  pipeline on the seeded synthetic fixture and verifies the seeded
  downstream genes are recovered.
- **Real-data example notebooks**: `cellranger_to_perturbflow.ipynb` and
  `replogle2022_real_data.ipynb`. The latter includes a known-biology
  recovery check on five hand-picked perturbations.
- **`CITATION.cff`** for `cff-version: 1.2.0`-style citation metadata.

### Fixed

- `qc.per_perturbation_qc` was averaging the on-target log2FC over **all**
  carrying cells, including Mixscape-NP (escaped) cells. The KO-only
  value is now reported as `on_target_log2fc`; the all-cells value is
  kept alongside as `on_target_log2fc_all_cells` for transparency.
- `perturbation_analysis._ensure_pca` was log-normalizing on top of
  already-log-normalized X, silently producing meaningless neighbor
  graphs. It now detects already-normalized X and skips the renorm.
- `de.run_pseudobulk_de` now skips-with-warning the case where every
  cell of a perturbation gets Mixscape-classified as NP (empty
  treatment arm), instead of crashing pydeseq2.

### Changed

- `qc.per_perturbation_qc` adds two new columns: `n_ko_cells` and
  `on_target_log2fc_all_cells` (see "Fixed" above).
- Snakefile rewritten as a real DAG (breaking change for anyone who was
  manually invoking the old single-rule wrapper; the high-level
  `perturbflow run` is unchanged).

## [0.1.0] - 2026-05-27

Initial public release.

- End-to-end Perturb-seq pipeline: IO → guide assignment → Mixscape
  KO/NP classification → pseudobulk DE → pathway scoring → HTML report.
- CellRanger H5 / MTX / h5ad readers; opinionated three-rule guide
  assignment; pertpy Mixscape wrapper with a GMM fallback; pseudobulk DE
  via pydeseq2 with a Welch fallback.
- Docker (multi-stage mamba) + AWS Batch deployment.
- GitHub Actions CI on Python 3.12 with full pertpy + pydeseq2 +
  decoupler installation.
- Methodology, architecture, and 10-recipe cookbook docs.
