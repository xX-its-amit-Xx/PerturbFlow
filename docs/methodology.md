# Methodology

This document explains the biological and statistical reasoning behind
each step of the pipeline. Read it once before you defend a PerturbFlow
run to a collaborator or a reviewer.

The pipeline is composed of five major stages: QC, guide assignment,
perturbation signal detection (Mixscape), pseudobulk DE, and downstream
interpretation (pathways + cell-state). Each is opinionated, and the
opinions are listed below.

---

## 1. QC

We compute the three standard single-cell QC metrics (total counts per
cell, genes detected per cell, % mitochondrial counts) but do **not**
apply a hard filter at this stage. The reason is that mitochondrial cutoffs
are tissue-specific (cardiomyocytes are 30%, T cells are 5%) and your
upstream pipeline almost certainly already filtered. We compute and report
them so the report has the numbers; thresholds are advisory.

If you want hard filtering, do it before feeding the AnnData to
PerturbFlow.

**On `mito_prefix`.** The default is `MT-` (human Ensembl convention). Set
it to `mt-` for mouse, `Mt` for Drosophila, `M` for some custom annotation,
or `""` to disable.

---

## 2. Guide assignment

Three sequential rules:

1. **UMI floor (`min_guide_umi`, default 5).** Guide capture libraries are
   intentionally sequenced shallow — they're a fraction of a percent of
   total reads. Below ~5 UMI per guide per cell, calls are dominated by
   ambient barcode swap. Lowering this rescues cells but at the cost of
   higher noise; raising it loses real cells.
2. **Dominance ratio (`dominance_ratio`, default 2.0).** If the top guide
   and runner-up are within a 2× ratio, the cell is **ambiguous**. This
   catches cells where two guides are near-tied — almost always
   contamination from library prep, not a real double infection.
3. **Multi-guide cap (`max_guides`, default 1).** Cells with more than one
   guide above the UMI floor are **multi-guide** (commonly doublets at
   MOI~0.3). The default rejects these; raise to 2 if your screen was
   designed for multi-perturbation per cell.

**The four-bucket status (`assigned`, `ambiguous`, `multi-guide`,
`unassigned`) is reported per-cell so you can audit assignment yield**, not
hidden behind a "passes_qc" boolean.

### When to use cellular indexing instead

If your screen used cellular indexing (CITE-seq HTOs marking each guide
clone, common in arrayed screens) rather than direct guide capture, run
your hashtag demultiplexer upstream (e.g.
[hashsolo](https://github.com/calico/solo)), write the `guide` column
into `adata.obs`, and feed the resulting `.h5ad` directly. PerturbFlow's
assignment module is bypassed when `guide` is already populated.

---

## 3. Mixscape: separating KO from escaped cells

### The problem

Even in well-validated CRISPRi/Cas9 screens, **10–40% of cells that carry a
guide RNA fail to perturb the target**. Reasons vary: incomplete editing,
in-frame indels, escape via the alternate allele, dCas9-KRAB depletion in
late G2/M, low protein turnover. These "escaped" cells dilute the DE
signal toward the null because they're indistinguishable from controls at
the protein level.

### The fix

[Mixscape (Papalexi et al. 2021)](https://www.nature.com/articles/s41588-021-00778-2)
models the per-cell perturbation signature as a two-component Gaussian
mixture (KO vs NP "non-perturbed") and classifies each cell. PerturbFlow
wraps [pertpy](https://pertpy.readthedocs.io)'s Python implementation.

The classification is then used to filter the DE arm: `use_mixscape_filter:
true` (default) drops NP cells from the treatment group when running
pseudobulk DE.

### What you give up

Mixscape requires enough cells per perturbation (>50 is comfortable, >100
is great) to fit a stable 2-component mixture. For sparse perturbations
(few cells per guide) the classification is unreliable and you may want to
disable Mixscape and accept the signal dilution. Set
`perturbation_analysis.enable_mixscape: false`.

For protein-level validation (a small panel of KOs with
flow-cytometry-confirmed knockdown), Mixscape often agrees ~85–90% with
the gold standard — see Figure 2 of Papalexi et al. for the validation
panel.

---

## 4. Pseudobulk DE

### Why pseudobulk and not single-cell DE

Single-cell DE methods (wilcoxon, MAST, etc.) treat each cell as an
independent replicate. They are not: cells from the same donor / clone /
biological replicate share most of their gene-expression variability. This
inflates p-values dramatically. [Squair et al. 2021](https://www.nature.com/articles/s41467-021-25960-2)
showed empirically that single-cell DE methods produce >10× more false
positives than pseudobulk methods on real benchmark datasets.

The fix is to **pool cells within a (perturbation × replicate) group, sum
their counts, and run DESeq2/edgeR on the resulting bulk-like table**.
PerturbFlow uses [pydeseq2](https://pydeseq2.readthedocs.io), a faithful
Python reimplementation of DESeq2.

### Real vs pseudo-replicates

If you have biological replicates (donors, time points, batches), set
`input.sample_col` to the column name. PerturbFlow will use those as the
replicate factor.

If you don't (a single-donor screen, common for early discovery work), the
pipeline synthesizes pseudo-replicates by hashing cell barcodes into
`n_pseudo_replicates` bins. **Pseudo-replicates do not buy you statistical
power that wasn't there.** They give DESeq2's dispersion estimator
something to chew on so it doesn't produce wildly conservative p-values,
but the design is still single-replicate per group. Treat the resulting
p-values as a ranking, not a calibrated test.

### Contrast

For each non-control perturbation we contrast `treatment` (KO cells from
that perturbation) vs `control` (all NT cells). Multi-factor designs
(perturbation × treatment, perturbation × time) are not currently
supported; if you need them, use `pydeseq2` directly on the pseudobulk
matrix that PerturbFlow writes to `de/pseudobulk_counts.parquet`.

---

## 5. Downstream interpretation

### Pathway scoring (decoupler-py)

For each perturbation, we feed the per-gene **test statistic** (Wald
statistic from DESeq2, or t-statistic from the fallback) into
[decoupler-py](https://decoupler-py.readthedocs.io)'s `run_ulm` against a
chosen gene-set collection. By default this is MSigDB Hallmarks (50 broad
biological processes); Reactome and PROGENy are also supported, and you
can pass a custom TSV.

**Why test statistic, not log2FC?** The statistic already accounts for
gene-level variance; log2FC alone overweights low-expression genes with
artifactually large fold-changes. Decoupler's `run_ulm` is a t-distributed
test, which is the right test to pair with a t-statistic.

### Cell-state effect map

For each perturbation we compute:

- **Centroid shift** — Euclidean distance in UMAP space between the
  perturbation's mean coordinate and the control mean.
- **Dispersion ratio** — perturbation's mean intra-group distance divided
  by the control's. Values >1 indicate the perturbation broadened
  cell-state heterogeneity (a common signature of differentiation
  perturbations).

These give you the figure-1 plot: "the cell-state landscape, and here's
how each KO moves cells through it". They are intentionally cheap to
compute; for the rigorous treatment (latent-space effects, generative
models) see [scGen](https://github.com/theislab/scgen) or
[CPA](https://github.com/theislab/cpa).

---

## 6. Reproducibility

Every run writes a `provenance.json` capturing:

- The PerturbFlow version and a 12-character SHA256 hash of the config
- Git revision and dirty flag
- Python version and platform
- Pinned versions of the major scientific dependencies

This file is also embedded in the HTML report. **Cite the config hash, not
the run date**, when referencing a specific analysis in a manuscript or
ticket.

Seeds: numpy and python `random` are seeded with `run.seed` (default 0);
torch is also seeded when present. UMAP is seeded via Scanpy. DESeq2 has
no random component.

---

## What PerturbFlow is not

- **Not a guide-counting tool.** Run CellRanger with CRISPR Guide Capture
  or the equivalent upstream and pass the guide-call CSV.
- **Not a batch-correction tool.** If you have donor batches, run Harmony,
  scVI, or BBKNN upstream and pass the resulting `.h5ad`.
- **Not a single-cell DE tool.** We use pseudobulk because we have to;
  if you want a cell-level statistic per gene per perturbation, run
  [Augur](https://github.com/neurorestore/Augur) on the same AnnData.
- **Not a continuous-effect classifier.** Mixscape labels are categorical
  (KO/NP). For per-cell continuous effect sizes use Mixscale or CINEMA-OT.
