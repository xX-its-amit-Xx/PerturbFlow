# PerturbFlow

> Opinionated, reproducible Perturb-seq analysis: from CellRanger output to
> defensible biological calls.

**Perturb-seq** is a pooled CRISPR screen with a single-cell RNA-seq readout.
Each cell receives one (or a few) guide RNAs from a library, so after
sequencing you know *what* each cell was perturbed with and *how* it
responded transcriptionally. PerturbFlow takes you from the raw 10x matrices
and guide-call CSV to per-perturbation differential expression, escape
detection, pathway scores, and a self-contained HTML report you can attach
to a project page.

This pipeline is for working bench-adjacent computational biologists who
want a reproducible, opinionated path from raw outputs to figures — not
another thin wrapper around Scanpy. The opinions are in three places: how
guides are assigned to cells, how escaped (non-perturbed) cells are
separated from true knock-outs, and how pseudobulk DE is set up. Each of
those choices is documented under [docs/methodology.md](docs/methodology.md)
so you can defend it to a reviewer.

[![CI](https://github.com/xX-its-amit-Xx/PerturbFlow/actions/workflows/ci.yml/badge.svg)](https://github.com/xX-its-amit-Xx/PerturbFlow/actions/workflows/ci.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](pyproject.toml)

---

## Quickstart

### Docker (recommended for production runs)

```bash
docker build -t perturbflow:latest -f docker/Dockerfile .

docker run --rm \
  -v $(pwd)/data:/workspace/data:ro \
  -v $(pwd)/results:/workspace/results \
  perturbflow:latest run --config /workspace/data/run.yaml --outdir /workspace/results
```

### Local install (development)

```bash
pip install -e ".[all,dev]"

perturbflow validate-config --config workflow/config.yaml
perturbflow run --config workflow/config.yaml
```

### Snakemake

```bash
snakemake -p --configfile workflow/config.yaml --cores 8
```

### AWS Batch

See [aws/README.md](aws/README.md). One `aws batch submit-job` call after
registering [aws/batch_job_definition.json](aws/batch_job_definition.json).

---

## Methodology: a decision tree

```
                            How were guides captured?
                              /                 \
                  Direct capture (CROP-seq,      Cellular indexing
                  10x CRISPR Guide Capture)      (CITE-seq HTO style)
                            |                              |
                  expression-based                 hashtag-demultiplex
                  guide assignment                 upstream (e.g. hashsolo),
                  (this pipeline)                  then enter pipeline with
                                                   guide already in obs

                            |
                Run Mixscape escape detection?
                          /            \
                  Per-target signal is        Bulk knock-down efficacy
                  uncertain (e.g. low-MOI,    is known to be high (>90%
                  partial CRISPRi knockdown)  protein loss confirmed)
                          |                            |
                  YES — separate KO from NP    Skip: marginal benefit,
                  cells before DE              extra compute
```

PerturbFlow defaults to **expression-based guide assignment** with the
three-rule procedure described in [docs/methodology.md#guide-assignment](docs/methodology.md#guide-assignment).
Mixscape escape detection is **on by default** because in practice we see
escape rates of 10–40% even in well-validated CRISPRi screens.

---

## Outputs

A run produces the following layout under `config.run.outdir`:

```
perturbflow/
├── qc/
│   ├── per_cell.csv                # one row per cell: counts, mito%, guide, mixscape class
│   ├── per_guide.csv               # one row per library guide: n_cells, mean UMI, escape rate
│   ├── per_perturbation.csv        # one row per perturbation: n_cells, escape, on-target log2FC
│   └── cell_state_effects.csv      # UMAP centroid shift + dispersion ratio per perturbation
├── de/
│   ├── <PERTURBATION>.csv          # DESeq2 pseudobulk results
│   ├── <PERTURBATION>.parquet      # same, faster to reload
│   └── pathway_scores.csv          # decoupler pathway activity, long format
├── figures/                        # PNGs of QC plots, volcanos (also embedded in report)
├── perturbflow.h5ad                # final AnnData with mixscape_class etc. for notebook follow-up
├── provenance.json                 # full run context: git rev, package versions, config hash
└── report.html                     # single-file interactive report
```

The CSVs and the `.h5ad` are the data products; the HTML is the
shareable artifact.

---

## Why not just use [X]?

| Tool | What PerturbFlow does differently |
|---|---|
| **Mixscape in Seurat** | Same idea, Python-native, packaged as a pipeline with config-driven runs and Docker/Batch deployment. PerturbFlow uses pertpy's pure-Python Mixscape implementation under the hood, so the *statistics* are the same — what's different is the surrounding pipeline (assignment, pseudobulk DE, report). |
| **Mixscale** | Mixscale is a continuous-effect model; PerturbFlow uses the discrete KO/NP labels because downstream pseudobulk DE wants categorical groups. If you need a continuous effect-size per cell, run Mixscale first and feed its output back into PerturbFlow via the `h5ad` entrypoint. |
| **pertpy alone** | pertpy is a library — PerturbFlow is a pipeline. Use pertpy if you want to script your own analysis; use PerturbFlow if you want a reproducible, config-driven run that produces a report and an AnnData. |
| **Scanpy / a custom notebook** | PerturbFlow forces explicit, defensible choices on guide assignment thresholds, Mixscape filtering, and pseudobulk DE design that you'd otherwise have to reinvent every project. |

PerturbFlow is *not* a guide-counting tool — assume you've already run
CellRanger with CRISPR Guide Capture or the equivalent. It's also not a
batch-correction tool — if you have multiple donors, bring your own
integration upstream and pass an `h5ad`.

---

## Documentation

- [Methodology](docs/methodology.md) — the *why* behind every step
- [Cookbook](docs/cookbook.md) — recipes for common real-world questions
- [Architecture](docs/architecture.md) — module layout and extension points
- [CellRanger → PerturbFlow notebook](examples/cellranger_to_perturbflow.ipynb)
  — zero-hand-rolled-CSV path from CellRanger 7.x output
- [Replogle 2022 real-data notebook](examples/replogle2022_real_data.ipynb)
  — runs the pipeline on the actual public dataset with a known-biology
  recovery check
- [Synthetic fixture walkthrough](examples/replogle2022_walkthrough.ipynb)
  — runnable end-to-end on the built-in test fixture, no download required

## Validity and correctness

PerturbFlow ships validity guards that catch the most common failure
modes other Perturb-seq pipelines silently produce wrong numbers for:

- **Raw-counts guard** — DE refuses to run on log-normalized X. Without
  this, the volcano looks fine but every gene is wrong.
- **Pseudo-replicate warning** — loud warning at the DE call site that
  hash-bin pseudo-replicates produce anticonservative p-values.
- **Mixscape-KO-only on-target log2FC** — the "did the screen work?"
  metric excludes escaped cells (which dilute the signal). The
  all-cells value is reported alongside for transparency.
- **Multi-guide concordance** — when running per-guide DE, the
  `multi_guide_concordance` helper flags guides whose effects are not
  reproducible across the per-gene library (off-target signal).
- **Parity with pertpy** — a test in CI verifies our Mixscape wrapper
  produces the same classifications as directly-called pertpy, so
  upstream API drift never silently changes our outputs.
- **Ground-truth recovery** — a CI test runs the full pipeline on a
  synthetic fixture with known downstream genes and verifies the
  pipeline finds them.

See [docs/methodology.md](docs/methodology.md) for the full statistical
justifications.

---

## Citation

If PerturbFlow helps your project, please cite it. Until the Zenodo
record is minted, cite the GitHub release tag plus the specific commit
hash you used:

```
Shenoy A. PerturbFlow: an opinionated Perturb-seq analysis pipeline.
GitHub (2026). https://github.com/xX-its-amit-Xx/PerturbFlow
Commit: <git rev-parse HEAD>
```

A machine-readable citation block is in [CITATION.cff](CITATION.cff).
The Zenodo DOI is **TBD** — it will be minted on the first tagged
release (`v0.1.0`), at which point this section will be updated with the
DOI and a permanent archive link.

When PerturbFlow drives a published analysis, please also cite the upstream
methods it composes:

- Papalexi *et al.* "Characterizing the molecular regulation of inhibitory
  immune checkpoints with multimodal single-cell screens." *Nature Genetics*
  53, 322–331 (2021).  — Mixscape
- Love *et al.* "Moderated estimation of fold change and dispersion for
  RNA-seq data with DESeq2." *Genome Biology* 15, 550 (2014). — DESeq2
- Badia-i-Mompel *et al.* "decoupleR: ensemble of computational methods to
  infer biological activities from omics data." *Bioinformatics Advances*
  2, vbac016 (2022). — pathway scoring

---

## Project status and roadmap

PerturbFlow follows semantic versioning. The current version is **0.1.0**,
which means the public API is stable for routine use but reserves the
right to evolve in 0.2.0+.

### Done (shipped in current `main`)

- Full pipeline: IO → guide assignment → Mixscape → pseudobulk DE →
  pathway scoring → HTML report.
- CellRanger 7.x adapter, Replogle 2022 adapter.
- Validity guards: raw-counts hard-fail, pseudo-replicate warning,
  KO-only on-target metric, multi-guide concordance.
- Per-stage Snakemake DAG with cached intermediates.
- Doublet detection via `scanpy.pp.scrublet` (optional, off by default).
- pertpy parity test + ground-truth recovery integration test in CI.
- Docker (mamba multi-stage) + AWS Batch job definition.
- Three example notebooks: synthetic walkthrough, CellRanger path,
  real-data Replogle 2022 with known-biology recovery check.

### In progress

- **PyPI publishing.** Release workflow ([release.yml](.github/workflows/release.yml))
  is wired and tested; the first publish will fire when `v0.1.0` is
  tagged and the `PYPI_API_TOKEN` repo secret is in place. Until then,
  install with `pip install git+https://github.com/xX-its-amit-Xx/PerturbFlow`.
- **Zenodo DOI.** Will be minted on the `v0.1.0` GitHub release via the
  Zenodo–GitHub integration. The DOI placeholder in the citation block
  above will be replaced at that point.

### TBD (not started; open for contribution)

- Continuous-effect classifier (Mixscale wrapper). Useful for CRISPRa
  screens where the per-cell effect is dose-dependent.
- Built-in batch correction integration (Harmony / scVI passes that run
  before guide assignment). Currently you bring your own and pass an
  `h5ad`.
- Native multi-factor DE design (perturbation × timepoint,
  perturbation × donor). Currently you split inputs and join post-hoc.
- Cluster execution profiles for Snakemake (Slurm, AWS Batch, K8s).
  The DAG is profile-ready; the profiles themselves are TBD.

If any of these block your work, please [open an issue](https://github.com/xX-its-amit-Xx/PerturbFlow/issues)
so they can be prioritized.

## License

[GPL v3.0](LICENSE).
