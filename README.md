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
- [Example notebook](examples/replogle2022_walkthrough.ipynb) — end-to-end on
  public data

---

## Citation

If PerturbFlow helps your project, please cite it:

```
Shenoy A. PerturbFlow: an opinionated Perturb-seq analysis pipeline.
Zenodo (2026). DOI: 10.5281/zenodo.0000000  [placeholder — release pending]
```

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

## License

[GPL v3.0](LICENSE).
