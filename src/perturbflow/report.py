"""Self-contained HTML report.

The report is one ``.html`` file with all images inline as base64 data
URIs. That format ships well: you can email it, drop it in Slack, attach
it to a Notion page; there are no broken-image issues because there are no
external assets.

The page layout: a header with the run summary, then a per-perturbation
card grid. Each card has the cells, escape fraction, on-target knockdown,
the top DE genes, and the top pathways.
"""

from __future__ import annotations

import base64
import io
import json
import logging
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from jinja2 import Template

logger = logging.getLogger(__name__)


_TEMPLATE = Template(
    """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>PerturbFlow report — {{ run_name }}</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
         margin: 0; padding: 0; color: #1f2933; background: #f7f9fc; }
  header { background: #102a43; color: #fff; padding: 1.5rem 2rem; }
  header h1 { margin: 0 0 0.25rem 0; font-size: 1.6rem; }
  header .meta { font-size: 0.85rem; opacity: 0.75; }
  main { padding: 1.5rem 2rem; max-width: 1200px; margin: 0 auto; }
  section { background: #fff; border-radius: 8px; padding: 1.25rem 1.5rem;
            margin-bottom: 1.5rem; box-shadow: 0 1px 3px rgba(16,42,67,0.08); }
  section h2 { margin: 0 0 1rem 0; font-size: 1.2rem; color: #102a43; }
  .stat-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
               gap: 0.75rem; }
  .stat { background: #f0f4f8; border-radius: 6px; padding: 0.75rem; }
  .stat .label { font-size: 0.75rem; text-transform: uppercase; color: #486581;
                 letter-spacing: 0.05em; }
  .stat .value { font-size: 1.4rem; font-weight: 600; color: #102a43; margin-top: 0.25rem; }
  table { width: 100%; border-collapse: collapse; font-size: 0.85rem; }
  th, td { text-align: left; padding: 0.4rem 0.6rem; border-bottom: 1px solid #e4e7eb; }
  th { background: #f0f4f8; }
  img { max-width: 100%; height: auto; border-radius: 4px; }
  .card-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(380px, 1fr));
               gap: 1rem; }
  .card { background: #f7f9fc; border: 1px solid #e4e7eb; border-radius: 6px; padding: 1rem; }
  .card h3 { margin: 0 0 0.5rem 0; font-size: 1.05rem; }
  .badge { display: inline-block; padding: 0.1rem 0.5rem; border-radius: 999px;
           font-size: 0.7rem; background: #d9e2ec; color: #243b53; margin-left: 0.5rem; }
  .badge.warn { background: #fcd9b6; color: #7c2d12; }
  pre { background: #f0f4f8; padding: 0.6rem; border-radius: 4px; overflow-x: auto;
        font-size: 0.78rem; }
  footer { text-align: center; padding: 1rem 2rem; color: #627d98; font-size: 0.8rem; }
</style>
</head>
<body>
<header>
  <h1>PerturbFlow — {{ run_name }}</h1>
  <div class="meta">
    perturbflow v{{ version }} · {{ generated_at }} · seed={{ seed }} · git={{ git_rev }}
  </div>
</header>
<main>
  <section>
    <h2>Run summary</h2>
    <div class="stat-grid">
      <div class="stat"><div class="label">Cells (input)</div><div class="value">{{ n_cells_in }}</div></div>
      <div class="stat"><div class="label">Cells (assigned)</div><div class="value">{{ n_cells_assigned }}</div></div>
      <div class="stat"><div class="label">Perturbations</div><div class="value">{{ n_perturbations }}</div></div>
      <div class="stat"><div class="label">Guides called</div><div class="value">{{ n_guides_called }}</div></div>
      <div class="stat"><div class="label">DE tests run</div><div class="value">{{ n_de_tests }}</div></div>
      <div class="stat"><div class="label">Mean escape rate</div><div class="value">{{ mean_escape }}</div></div>
    </div>
  </section>

  {% if assignment_plot %}
  <section>
    <h2>Guide assignment</h2>
    <img src="data:image/png;base64,{{ assignment_plot }}" alt="Guide assignment status">
  </section>
  {% endif %}

  {% if umap_plot %}
  <section>
    <h2>Cell-state landscape</h2>
    <img src="data:image/png;base64,{{ umap_plot }}" alt="UMAP coloured by perturbation">
  </section>
  {% endif %}

  <section>
    <h2>Per-perturbation cards</h2>
    <div class="card-grid">
    {% for card in cards %}
      <div class="card">
        <h3>{{ card.perturbation }}
          {% if card.escape_fraction is not none and card.escape_fraction > 0.4 %}
            <span class="badge warn">high escape</span>
          {% endif %}
        </h3>
        <div class="stat-grid" style="grid-template-columns: repeat(3, 1fr);">
          <div class="stat"><div class="label">Cells</div><div class="value">{{ card.n_cells }}</div></div>
          <div class="stat"><div class="label">Escape</div><div class="value">{{ card.escape_display }}</div></div>
          <div class="stat"><div class="label">On-target log2FC</div><div class="value">{{ card.on_target_display }}</div></div>
        </div>
        {% if card.volcano_plot %}
          <img src="data:image/png;base64,{{ card.volcano_plot }}" alt="Volcano for {{ card.perturbation }}">
        {% endif %}
        <h4 style="margin-bottom: 0.25rem;">Top DE genes</h4>
        {% if card.top_de_rows %}
        <table>
          <thead><tr><th>gene</th><th>log2FC</th><th>padj</th></tr></thead>
          <tbody>
          {% for row in card.top_de_rows %}
            <tr><td>{{ row.gene }}</td><td>{{ "%.2f"|format(row.log2FoldChange) }}</td><td>{{ "%.2e"|format(row.padj) }}</td></tr>
          {% endfor %}
          </tbody>
        </table>
        {% else %}
        <p style="color: #627d98; font-size: 0.85rem;">No DE results for this perturbation.</p>
        {% endif %}
        {% if card.top_pathways %}
        <h4 style="margin-bottom: 0.25rem;">Top pathways</h4>
        <table>
          <thead><tr><th>pathway</th><th>score</th><th>pvalue</th></tr></thead>
          <tbody>
          {% for row in card.top_pathways %}
            <tr><td>{{ row.pathway }}</td><td>{{ "%.2f"|format(row.score) }}</td><td>{{ "%.2e"|format(row.pvalue) }}</td></tr>
          {% endfor %}
          </tbody>
        </table>
        {% endif %}
      </div>
    {% endfor %}
    </div>
  </section>

  <section>
    <h2>Provenance</h2>
    <pre>{{ provenance_json }}</pre>
  </section>
</main>
<footer>
  Generated by <a href="https://github.com/xX-its-amit-Xx/PerturbFlow">PerturbFlow</a>.
</footer>
</body>
</html>
"""
)


def write_html_report(
    *,
    out_path: str | Path,
    run_name: str,
    seed: int,
    version: str,
    git_rev: str,
    per_cell: pd.DataFrame,
    per_perturbation: pd.DataFrame,
    de_results: dict[str, pd.DataFrame],
    pathway_scores: pd.DataFrame,
    cell_state: pd.DataFrame | None = None,
    umap_coords: np.ndarray | None = None,
    umap_labels: pd.Series | None = None,
    provenance: dict[str, Any] | None = None,
) -> Path:
    """Render a single-file HTML report and write it to ``out_path``."""
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    assignment_plot = _render_assignment_plot(per_cell)
    umap_plot = (
        _render_umap_plot(umap_coords, umap_labels)
        if umap_coords is not None and umap_labels is not None
        else None
    )

    cards = []
    for _, prow in per_perturbation.iterrows():
        pert = str(prow["perturbation"])
        de_df = de_results.get(pert)
        top_de_rows: list[dict[str, Any]] = []
        if de_df is not None and not de_df.empty:
            ranked = de_df.assign(_abs=de_df["log2FoldChange"].abs()).sort_values(
                ["padj", "_abs"], ascending=[True, False]
            )
            top_de_rows = ranked.head(10)[["gene", "log2FoldChange", "padj"]].to_dict(
                orient="records"
            )
        top_pathways: list[dict[str, Any]] = []
        if not pathway_scores.empty:
            sub = pathway_scores[pathway_scores["perturbation"] == pert].head(5)
            top_pathways = sub[["pathway", "score", "pvalue"]].to_dict(orient="records")
        cards.append(
            {
                "perturbation": pert,
                "n_cells": int(prow["n_cells"]),
                "escape_fraction": _coerce_float(prow.get("escape_fraction")),
                "escape_display": _fmt_pct(prow.get("escape_fraction")),
                "on_target_display": _fmt_lfc(prow.get("on_target_log2fc")),
                "top_de_rows": top_de_rows,
                "top_pathways": top_pathways,
                "volcano_plot": _render_volcano(de_df, pert) if de_df is not None else None,
            }
        )

    mean_escape = (
        f"{per_perturbation['escape_fraction'].dropna().mean() * 100:.1f}%"
        if "escape_fraction" in per_perturbation.columns
        and per_perturbation["escape_fraction"].notna().any()
        else "—"
    )

    prov_json = json.dumps(provenance or {}, indent=2, sort_keys=True, default=str)

    html = _TEMPLATE.render(
        run_name=run_name,
        seed=seed,
        version=version,
        git_rev=git_rev,
        generated_at=pd.Timestamp.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC"),
        n_cells_in=len(per_cell),
        n_cells_assigned=int(per_cell["assignment_status"].astype(str).eq("assigned").sum())
        if "assignment_status" in per_cell.columns
        else len(per_cell),
        n_perturbations=int(per_perturbation.shape[0]),
        n_guides_called=int(per_cell["guide"].dropna().nunique())
        if "guide" in per_cell.columns
        else 0,
        n_de_tests=len(de_results),
        mean_escape=mean_escape,
        assignment_plot=assignment_plot,
        umap_plot=umap_plot,
        cards=cards,
        provenance_json=prov_json,
    )
    out.write_text(html, encoding="utf-8")
    logger.info("Wrote HTML report: %s", out)
    return out


# ----------------------------------------------------------------------
# Plot helpers (all return base64-encoded PNG strings, or None on failure)
# ----------------------------------------------------------------------


def _to_data_uri(fig: plt.Figure) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _render_assignment_plot(per_cell: pd.DataFrame) -> str | None:
    if "assignment_status" not in per_cell.columns:
        return None
    counts = per_cell["assignment_status"].astype(str).value_counts()
    fig, ax = plt.subplots(figsize=(5, 3))
    colors = {
        "assigned": "#2e7d32",
        "ambiguous": "#f9a825",
        "multi-guide": "#c62828",
        "unassigned": "#9e9e9e",
    }
    ax.bar(counts.index, counts.values, color=[colors.get(k, "#607d8b") for k in counts.index])
    ax.set_ylabel("Cells")
    ax.set_title("Guide assignment status")
    for i, v in enumerate(counts.values):
        ax.text(i, v, str(int(v)), ha="center", va="bottom", fontsize=8)
    return _to_data_uri(fig)


def _render_umap_plot(coords: np.ndarray, labels: pd.Series) -> str | None:
    if coords.shape[0] != len(labels):
        return None
    fig, ax = plt.subplots(figsize=(6, 5))
    cats = labels.astype(str)
    uniq = sorted(cats.unique())
    cmap = plt.get_cmap("tab20", max(len(uniq), 1))
    for i, lbl in enumerate(uniq):
        m = (cats == lbl).values
        ax.scatter(coords[m, 0], coords[m, 1], s=4, alpha=0.7, label=lbl, color=cmap(i))
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("Cell-state landscape, coloured by perturbation")
    ax.legend(fontsize=7, loc="best", markerscale=2)
    return _to_data_uri(fig)


def _render_volcano(de_df: pd.DataFrame | None, pert: str) -> str | None:
    if de_df is None or de_df.empty:
        return None
    fig, ax = plt.subplots(figsize=(4.5, 3.5))
    pvals = de_df["padj"].clip(lower=1e-300)
    nlog10 = -np.log10(pvals)
    sig = de_df.get("significant", pd.Series(False, index=de_df.index))
    ax.scatter(de_df["log2FoldChange"], nlog10, s=6, alpha=0.4, c="#9e9e9e", label="ns")
    if sig.any():
        ax.scatter(
            de_df.loc[sig, "log2FoldChange"],
            nlog10[sig],
            s=8,
            alpha=0.85,
            c="#c62828",
            label="significant",
        )
    ax.set_xlabel("log2 fold-change")
    ax.set_ylabel("-log10(padj)")
    ax.set_title(f"Volcano: {pert}")
    ax.axhline(-np.log10(0.05), color="#90a4ae", linestyle="--", linewidth=0.7)
    ax.legend(fontsize=7)
    return _to_data_uri(fig)


def _coerce_float(v: object) -> float | None:
    try:
        f = float(v)  # type: ignore[arg-type]
        return None if np.isnan(f) else f
    except (TypeError, ValueError):
        return None


def _fmt_pct(v: object) -> str:
    f = _coerce_float(v)
    return f"{f * 100:.1f}%" if f is not None else "—"


def _fmt_lfc(v: object) -> str:
    f = _coerce_float(v)
    return f"{f:+.2f}" if f is not None else "—"
