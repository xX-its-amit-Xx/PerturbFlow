# Architecture

A short reference for navigating the codebase and extending the pipeline.

## Module map

```
src/perturbflow/
├── __init__.py                # Public re-exports
├── _version.py                # Single-source version string
├── config.py                  # Strongly-typed YAML loader (frozen dataclasses)
├── io.py                      # 10x H5/MTX/h5ad readers, guide CSV schema validation
├── guide_assignment.py        # Three-rule per-cell guide call
├── perturbation_analysis.py   # pertpy Mixscape wrapper + GMM fallback
├── qc.py                      # Per-cell / per-guide / per-perturbation QC tables
├── de.py                      # Pseudobulk + pydeseq2 (Welch fallback)
├── downstream.py              # decoupler pathway scoring + UMAP centroid shifts
├── report.py                  # Jinja2 single-file HTML, base64-inlined PNGs
├── provenance.py              # Git rev / package versions / config hash
├── pipeline.py                # The composed `run()` (CLI + Snakemake call into it)
└── cli.py                     # `perturbflow {run,validate-config,version}` click app
```

## Pipeline data flow

```
config.yaml                                      [validated by config.py]
     │
     ▼
io.read_10x_h5 / read_10x_mtx / read_h5ad        ── AnnData
io.read_guide_calls                              ── long-format DataFrame
io.read_guide_metadata                           ── guide_id → target_gene map
     │
     ▼
guide_assignment.assign_guides                   ── AnnData with .obs[guide, perturbation, assignment_status]
     │
     ▼
qc.per_cell_qc / per_guide_qc                    ── DataFrames written to qc/
     │
     ▼
perturbation_analysis.compute_perturbation_signature   ── AnnData with .layers[X_pert]
perturbation_analysis.run_mixscape               ── AnnData with .obs[mixscape_class, mixscape_perturbed]
     │
     ▼
qc.per_perturbation_qc                           ── DataFrame (with escape_fraction, on_target_log2fc)
     │
     ▼
pipeline._ensure_embedding                       ── AnnData with .obsm[X_umap]
downstream.compute_cell_state_effects            ── DataFrame (centroid_shift, dispersion_ratio)
     │
     ▼
de.run_pseudobulk_de                             ── dict[perturbation -> DataFrame] written to de/
     │
     ▼
downstream.score_pathways                        ── DataFrame written to de/pathway_scores.csv
     │
     ▼
report.write_html_report                         ── perturbflow/report.html
```

## Extension points

### Adding a new pathway collection

`downstream._load_pathway_network` resolves named collections (`hallmarks`,
`reactome`, `progeny`) and falls through to file-path mode for anything
else. To add a new built-in:

```python
# downstream.py
if key == "kegg":
    net = dc.get_resource("KEGG", organism=organism)
    return net.rename(columns={"geneset": "source", "genesymbol": "target"})[
        ["source", "target", "weight"]
    ]
```

…and document the new value under `downstream.pathway_net` in
`workflow/config.yaml`.

### Replacing the DE backend

`de.run_pseudobulk_de` checks `find_spec("pydeseq2")` and falls back to a
Welch t-test on log-CPM. To plug in a different backend (edgeR via
rpy2, limma-voom, etc.):

```python
# de.py
def _run_edger(pb: ad.AnnData, *, pert: str, control: str) -> pd.DataFrame:
    ...
    return pd.DataFrame(
        {"gene", "baseMean", "log2FoldChange", "stat", "pvalue", "padj"}
    )

# In run_pseudobulk_de:
if cfg.backend == "edger":
    df = _run_edger(sub, pert=pert, control=control_label)
```

Add a `backend: str = "pydeseq2"` field to `DEConfig` and document under
`de.backend`.

### Replacing Mixscape

The wrapper in `perturbation_analysis.run_mixscape` only requires that
the result populate `adata.obs["mixscape_perturbed"]` (bool, True for
true-KO). Anything that satisfies that contract — `Mixscale`, `CINEMA-OT`,
a custom GMM — is a drop-in.

### Hooking custom QC plots into the report

The report template (`report._TEMPLATE`) iterates over `cards`, one per
perturbation. Each card carries an arbitrary dict so adding a new image is
two lines:

```python
# In write_html_report:
cards.append({..., "my_custom_plot": _render_my_thing(de_df, pert)})

# In _TEMPLATE:
{% if card.my_custom_plot %}
  <img src="data:image/png;base64,{{ card.my_custom_plot }}" alt="…">
{% endif %}
```

## Configuration philosophy

The config is one YAML file, validated against frozen dataclasses with
**strict unknown-key rejection**. The cost of a typo in a long-running
batch job is too high to silently ignore it. If you want to add a new
knob:

1. Add a field to the relevant `*Config` dataclass in `config.py`.
2. Default it sensibly (the cluster-friendly default, not the
   experiment-friendly default).
3. Document it in `workflow/config.yaml` with a comment explaining
   *what it does* and *why you'd change it*.
4. Use it in the module that reads it. If the module doesn't already
   take the config dataclass, add it as a kwarg.

## Reproducibility contract

Every run produces a `provenance.json` with:

- `perturbflow_version`
- `config_hash` (SHA256 of the resolved config, truncated to 12 chars)
- `git.revision`, `git.dirty`
- `python`, `platform`
- `packages` (versions of the scientific stack)

This file is also embedded into the HTML report. **Manuscript citation
practice: reference the run by `(perturbflow_version, config_hash,
git_revision)`** — those three uniquely identify the analysis even if the
output directory is lost.
