# Cookbook

Recipes for real-world Perturb-seq questions. Each recipe is a short
config snippet plus a paragraph on when to use it.

If the recipe you need isn't here, open an issue —
[github.com/xX-its-amit-Xx/PerturbFlow/issues](https://github.com/xX-its-amit-Xx/PerturbFlow/issues).

---

## 1. Genome-scale CRISPRi screen (Replogle 2022 style)

**Question:** You have a K562 CRISPRi screen covering ~2,000 essential genes.
You want a defensible DE table per gene to feed into a downstream
gene–program clustering analysis.

**Why it's hard:** With ~2,000 perturbations and 30–50 cells per guide,
single-cell DE is wildly underpowered and noisy. Pseudobulk DE with
pseudo-replicates is the only viable approach at this scale, and you need
Mixscape to cut escape noise — at this depth, even a 10% escape rate is
enough to flatten log2FCs to within noise.

```yaml
run:
  name: replogle_2022_essential
  outdir: results/replogle_essential
  seed: 0

input:
  matrix_h5: data/replogle_essential.h5
  guide_calls: data/guide_calls.csv
  guide_metadata: data/guide_library.csv
  sample_col: null              # single donor — synthesize replicates
  n_pseudo_replicates: 3

guide_assignment:
  min_guide_umi: 5
  dominance_ratio: 3.0          # tighter — CRISPRi guides cross-contaminate
  max_guides: 1
  drop_unassigned: true

perturbation_analysis:
  enable_mixscape: true
  control_label: "non-targeting"
  n_neighbors: 25
  mixscape_pval_cutoff: 5.0e-2

de:
  enable: true
  min_replicates_per_group: 2
  min_cells_per_replicate: 5    # low — screens at this scale are sparse
  lfc_threshold: 0.5
  use_mixscape_filter: true

downstream:
  enable_pathway_scoring: true
  pathway_net: hallmarks
```

**Output to look at first:** `qc/per_perturbation.csv` — sort by
`escape_fraction` to find perturbations where Mixscape couldn't find a
signal. A high escape rate on a known-essential gene means the guide didn't
work or the cells are dying too fast to capture.

---

## 2. Small targeted CRISPR-KO screen (≤50 perturbations, multi-donor)

**Question:** You have a focused screen of 30 candidate tumor-suppressor
genes in CD8 T cells from 4 donors. You want donor-controlled DE so you
can distinguish perturbation effects from donor biology.

**Why it's hard:** Donor effects in primary cells are enormous (often >
perturbation effects). You need to model donors as the replicate factor,
not synthesize pseudo-replicates.

```yaml
run:
  name: cd8_tumor_suppressor_screen
  outdir: results/cd8_screen

input:
  h5ad: data/cd8_screen_postharmony.h5ad   # already batch-corrected
  guide_calls: data/guide_calls.csv
  guide_metadata: data/guide_library.csv
  sample_col: donor_id                      # the real biological replicate
  n_pseudo_replicates: 1                    # ignored when sample_col is set

guide_assignment:
  min_guide_umi: 3                          # primary T cells: lower guide capture
  dominance_ratio: 2.0
  max_guides: 1

perturbation_analysis:
  enable_mixscape: true
  control_label: "scrambled"
  n_neighbors: 20

de:
  enable: true
  min_replicates_per_group: 3               # need at least 3 donors per perturbation
  min_cells_per_replicate: 20
  use_mixscape_filter: true
  lfc_threshold: 1.0
  padj_threshold: 0.05
```

**Output to look at first:** the per-perturbation card in the HTML
report. If `n_cells` < 60 (3 donors × 20 cells/donor) for any perturbation,
flag it — DESeq2 will run but the effect estimates are noisy.

---

## 3. CRISPRa screen with continuous dose

**Question:** You ran a CRISPRa screen using dCas9-VPR with multiple guides
per gene of varying activation strength. You want to use only the strongest
activation cells (top quartile of target gene expression) for DE.

**Why it's hard:** CRISPRa effects vary continuously per cell. Mixscape's
KO/NP binarization is a poor fit; you want to filter by *level*, not by
discrete class. The recipe: pre-filter the AnnData to top-quartile cells
per perturbation outside PerturbFlow, then run with Mixscape disabled.

```python
# Pre-filter (outside PerturbFlow)
import scanpy as sc
import perturbflow as pf

adata = pf.read_10x_h5("crispr_a.h5")
guides = pf.read_guide_calls("guide_calls.csv")
meta = pf.read_guide_metadata("guide_library.csv")
adata = pf.assign_guides(adata, guides, meta)

# Keep only cells where the target gene is in the top quartile of its perturbation.
import numpy as np
keep_mask = np.zeros(adata.n_obs, dtype=bool)
for pert in adata.obs["perturbation"].unique():
    if pert in {"NT", "unassigned"}:
        keep_mask |= adata.obs["perturbation"] == pert  # always keep controls
        continue
    if pert not in adata.var_names:
        continue
    pert_mask = (adata.obs["perturbation"] == pert).values
    target_expr = adata[:, pert].X.toarray().ravel()
    cutoff = np.quantile(target_expr[pert_mask], 0.75)
    keep_mask |= pert_mask & (target_expr >= cutoff)

adata = adata[keep_mask].copy()
adata.write_h5ad("crispr_a_topquartile.h5ad")
```

```yaml
# config.yaml — point at the pre-filtered h5ad
input:
  h5ad: data/crispr_a_topquartile.h5ad
  guide_calls: data/guide_calls.csv
  guide_metadata: data/guide_library.csv

perturbation_analysis:
  enable_mixscape: false        # binary KO/NP doesn't fit continuous CRISPRa

de:
  enable: true
  use_mixscape_filter: false    # we already filtered
```

---

## 4. Re-running just the DE step on a previous run

**Question:** You ran the pipeline yesterday, the report looks fine, but
you realize you wanted `lfc_threshold: 0.5` instead of `1.0`. You don't
want to redo guide assignment + Mixscape (slow).

**Why it's easy:** PerturbFlow writes the final AnnData to
`perturbflow.h5ad`. Edit the config and re-run from that.

```yaml
input:
  h5ad: results/replogle_essential/perturbflow.h5ad   # last run's output
  guide_calls: data/guide_calls.csv
  guide_metadata: data/guide_library.csv

guide_assignment:
  drop_unassigned: false        # already dropped; this is a no-op

perturbation_analysis:
  enable_mixscape: false        # mixscape_class already in obs

de:
  enable: true
  lfc_threshold: 0.5            # the change you wanted
  use_mixscape_filter: true
```

The first three stages become near-no-ops because the AnnData already has
the right columns. **The pipeline does NOT detect this automatically** —
you have to set `enable_mixscape: false` yourself. The Snakemake
workflow is a better fit for fine-grained re-running (see `workflow/Snakefile`).

---

## 5. Screen with multiple time points

**Question:** You have a CRISPRi screen with samples collected at 24h, 48h,
and 72h. You want a DE table per perturbation per timepoint.

**Why it's hard:** PerturbFlow's current contrast model is binary
(treatment vs control). For per-timepoint DE, run the pipeline three times,
once per timepoint, and join the results.

```bash
# Split your h5ad upstream
python -c "
import anndata as ad
a = ad.read_h5ad('data/timeseries.h5ad')
for tp in ['24h', '48h', '72h']:
    a[a.obs['timepoint'] == tp].write_h5ad(f'data/{tp}.h5ad')
"

# One config per timepoint
for tp in 24h 48h 72h; do
  cp configs/template.yaml configs/${tp}.yaml
  sed -i "s|INPUT_H5AD|data/${tp}.h5ad|; s|RUN_NAME|screen_${tp}|; s|OUTDIR|results/${tp}|" \
    configs/${tp}.yaml
  perturbflow run --config configs/${tp}.yaml
done
```

Then post-process the three per-perturbation DE tables in a notebook to
build a (gene × perturbation × timepoint) tensor for trajectory analysis.

---

## 6. Sanity-checking a published Perturb-seq result

**Question:** A collaborator sent you a list of "high-confidence hits" from
a Perturb-seq paper. You want to re-analyze the public dataset to confirm
the hits actually replicate under defensible methodology.

**Recipe:**

```bash
# 1. Get the public data (Replogle 2022 GEO accession GSE177150 or figshare).
mkdir -p data/replogle_2022
# … download .h5 + guide call CSV + guide library …

# 2. Run with strict defaults.
perturbflow run --config configs/replogle_strict.yaml

# 3. Check overlap with reported hits.
python -c "
import pandas as pd
reported = set(open('collaborator_hits.txt').read().split())
de = pd.read_csv('results/replogle_strict/de/TARGET_GENE.csv')
ours = set(de[(de['padj'] < 0.05) & (de['log2FoldChange'].abs() > 1)]['gene'])
print('Concordance:', len(reported & ours) / len(reported))
print('Missing from our call:', reported - ours)
"
```

**What to look for:**

- A concordance below 60% with strict defaults usually means the paper
  used single-cell DE or skipped Mixscape.
- "Missing from our call" genes worth investigating manually — are they
  in the lower-padj tail (real but borderline) or completely absent?

---

## 7. Producing a single per-perturbation report card for a slide

**Question:** You want one summary PNG per perturbation that you can paste
into a presentation.

**Recipe:** The HTML report already contains a card per perturbation. To
extract them as standalone PNGs:

```python
import perturbflow as pf
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

de_dir = Path("results/screen/de")
for csv_path in de_dir.glob("*.csv"):
    pert = csv_path.stem
    df = pd.read_csv(csv_path)
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.scatter(df["log2FoldChange"], -df["padj"].clip(1e-300).apply(lambda x: -__import__("math").log10(x)), s=4, alpha=0.4)
    ax.set_title(f"{pert} volcano")
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10(padj)")
    fig.savefig(f"figures/{pert}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
```

---

## 8. Adding a custom QC gate (e.g., minimum on-target knockdown)

**Question:** You want to exclude any perturbation whose on-target log2FC
is weaker than -1.0 — those guides clearly didn't work and including
them in the report just adds noise.

**Recipe:** Filter the per-perturbation table downstream:

```python
import pandas as pd
qc = pd.read_csv("results/screen/qc/per_perturbation.csv")
worked = qc[qc["on_target_log2fc"] < -1.0]
print(f"{len(worked)}/{len(qc)} perturbations passed the knockdown gate")
worked.to_csv("results/screen/qc/per_perturbation_passing.csv", index=False)
```

To make this hard-fail on a future config, add a post-pipeline check to
your Snakefile:

```python
rule on_target_gate:
    input: "results/screen/qc/per_perturbation.csv"
    output: "results/screen/qc/passing_perturbations.csv"
    run:
        import pandas as pd
        df = pd.read_csv(input[0])
        passing = df[df["on_target_log2fc"] < -1.0]
        if len(passing) < 0.5 * len(df):
            raise ValueError("Less than 50% of perturbations passed on-target gate; aborting")
        passing.to_csv(output[0], index=False)
```

---

## 9. Comparing two Perturb-seq screens

**Question:** You ran the same library in two cell lines (K562 vs U2OS) and
want to compare effect sizes per perturbation.

**Recipe:** Run PerturbFlow twice, once per cell line, then merge the
per-gene log2FCs:

```python
import pandas as pd
from pathlib import Path

def load_all(run_dir: Path) -> pd.DataFrame:
    frames = []
    for csv in (run_dir / "de").glob("*.csv"):
        pert = csv.stem
        if pert == "pathway_scores":
            continue
        df = pd.read_csv(csv)[["gene", "log2FoldChange", "padj"]]
        df["perturbation"] = pert
        frames.append(df)
    return pd.concat(frames, ignore_index=True)

k562 = load_all(Path("results/k562")).rename(columns={"log2FoldChange": "lfc_k562", "padj": "padj_k562"})
u2os = load_all(Path("results/u2os")).rename(columns={"log2FoldChange": "lfc_u2os", "padj": "padj_u2os"})
both = k562.merge(u2os, on=["gene", "perturbation"], how="outer")

# Conserved hits: significant in both with same direction
conserved = both[
    (both["padj_k562"] < 0.05)
    & (both["padj_u2os"] < 0.05)
    & (both["lfc_k562"] * both["lfc_u2os"] > 0)
]
conserved.to_csv("conserved_hits.csv", index=False)
```

---

## 10. Profiling a slow run

**Question:** A 10k-cell screen is taking 20 minutes per run and you want
to know where the time is going.

**Recipe:**

```bash
# Drop into Python with cProfile
python -m cProfile -o profile.out -m perturbflow.cli run --config workflow/config.yaml
python -c "
import pstats
p = pstats.Stats('profile.out')
p.strip_dirs().sort_stats('cumulative').print_stats(30)
"
```

**Usual culprits, in rough rank order:**

1. UMAP — Scanpy's UMAP is single-threaded and dominant for >50k cells.
   Set `OMP_NUM_THREADS=1` and use `rapids-singlecell` if you have a GPU.
2. PCA in `_ensure_pca` — `sc.tl.pca` is reasonably fast but
   `sc.pp.scale` makes a dense copy that can OOM on large screens.
3. Mixscape — pertpy's `perturbation_signature` builds a k-NN index per
   call. For genome-scale screens (>10k cells, >2k perturbations) this
   dominates.

If you're below 20k cells, the time is almost certainly in DESeq2; raise
`min_cells_per_replicate` so fewer pseudobulks get tested.
