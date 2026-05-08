"""Self-contained HTML report builder.

Driven by a small TOML inventory (`sagebench-report.toml`). For each
listed run we evaluate at three levels (ion / peptide / protein), then
render headline tables, calibration plots, coverage plots, and three
sets of overlap Venns. A configurable `[deep_dive]` section adds a
sagepy-TDC vs sage-native comparison at each level.
"""
from __future__ import annotations

import html
import sys
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover
    import tomli as tomllib

from . import ground_truth, metrics, plots, sage_output
from . import tdc as tdc_mod
from .ground_truth import LEVELS, Level

_DEFAULT_Q = (0.001, 0.005, 0.01, 0.05)
_ROC_POINTS = 250
_VENN_Q = 0.01


@dataclass
class _LevelView:
    """Per-level evaluation slice for one run."""

    sweep: pd.DataFrame
    roc: pd.DataFrame
    headline: dict
    keys_at_1pct: set
    n_truth: int
    n_fragmentable: int


@dataclass
class _Section:
    """All three level-views for one run."""

    name: str
    n_psms: int
    n_targets: int
    n_decoys: int
    sage_parquet: str
    truth_paths: list[str]
    levels: dict[Level, _LevelView] = field(default_factory=dict)


def _load_section(spec: dict) -> _Section:
    truths = [ground_truth.load(p) for p in spec["truth"]]
    truth = ground_truth.union(*truths)
    hits = sage_output.load(spec["sage_out"])

    section = _Section(
        name=spec["name"],
        n_psms=hits.n_psms,
        n_targets=hits.n_targets,
        n_decoys=hits.n_psms - hits.n_targets,
        sage_parquet=spec["sage_out"],
        truth_paths=list(spec["truth"]),
    )
    for level in LEVELS:
        sw = metrics.sweep(hits, truth, q_cutoffs=list(_DEFAULT_Q), level=level)
        roc = metrics.roc_points(hits, truth, n_points=_ROC_POINTS, level=level)
        h_row = sw[sw["q_cutoff"] == 0.01]
        headline = h_row.iloc[0].to_dict() if len(h_row) else {}
        keys = metrics.precursor_keys_at(hits, q_cutoff=_VENN_Q, level=level)
        section.levels[level] = _LevelView(
            sweep=sw,
            roc=roc,
            headline=headline,
            keys_at_1pct=keys,
            n_truth=truth.n_at(level),
            n_fragmentable=truth.n_fragmentable_at(level),
        )
    return section


# ---------------------------------------------------------------------------
# HTML helpers
# ---------------------------------------------------------------------------


def _summary_table(sections: list[_Section], level: Level) -> str:
    rows = []
    for s in sections:
        v = s.levels[level]
        h = v.headline
        rows.append(
            {
                "run": s.name,
                "n_truth": v.n_truth,
                "n_frag": v.n_fragmentable,
                f"{level}s @ q≤0.01": int(h.get("reported_keys", 0)),
                "TP": int(h.get("tp", 0)),
                "FP": int(h.get("fp", 0)),
                "true FDR": f"{h.get('true_fdr', 0.0) * 100:.2f} %",
                "TPR (full)": f"{h.get('tpr', 0.0) * 100:.2f} %",
                "TPR (frag)": f"{h.get('tpr_fragmentable', 0.0) * 100:.2f} %",
            }
        )
    return pd.DataFrame(rows).to_html(index=False, classes="summary",
                                       border=0, escape=False)


def _sweep_table(view: _LevelView) -> str:
    df = view.sweep.copy()
    df["true_fdr"] = df["true_fdr"].map(lambda v: f"{v * 100:.2f} %")
    df["tpr"] = df["tpr"].map(lambda v: f"{v * 100:.2f} %")
    df["tpr_fragmentable"] = df["tpr_fragmentable"].map(lambda v: f"{v * 100:.2f} %")
    df = df.rename(
        columns={
            "q_cutoff": "q ≤",
            "reported_psms": "PSMs",
            "reported_keys": "keys",
            "tpr": "TPR (full)",
            "tpr_fragmentable": "TPR (frag)",
        }
    )
    return df[["q ≤", "PSMs", "keys", "tp", "fp",
               "true_fdr", "TPR (full)", "TPR (frag)"]].to_html(
        index=False, classes="sweep", border=0, escape=False
    )


def _format_notes(raw: str | None) -> str:
    if not raw:
        return ""
    parts = [p.strip() for p in raw.split("\n\n") if p.strip()]
    return "\n".join(
        f"<p>{html.escape(p).replace(chr(10), '<br>')}</p>" for p in parts
    )


def _level_block(level: Level, sections: list[_Section]) -> str:
    """Headline + plots + venns for a single level."""
    summary = _summary_table(sections, level)

    cal_curves = [plots.RunCurve(name=s.name, roc=s.levels[level].roc) for s in sections]
    cal_svg = plots.calibration(cal_curves, log=True)
    cov_svg = plots.coverage(cal_curves, y="tp")
    frag_svg = plots.coverage(cal_curves, y="tpr_fragmentable")

    return (
        f"<details open><summary><strong>Level: {level}</strong></summary>"
        f"<h3>Headline (q ≤ 0.01)</h3>"
        f"{summary}"
        f"<div class='plots'>"
        f"<div class='plot'>{cal_svg}</div>"
        f"<div class='plot'>{cov_svg}</div>"
        f"</div>"
        f"<div class='plots' style='grid-template-columns: 1fr;'>"
        f"<div class='plot'>{frag_svg}</div>"
        f"</div>"
        f"</details>"
    )


def _venn_block(level: Level, sections: list[_Section], specs: list[dict]) -> str:
    by_name = {s.name: s for s in sections}
    blocks = []
    for spec in specs:
        run_names = spec["runs"]
        if len(run_names) != 3:
            raise ValueError(f"venn '{spec.get('title','?')}' needs exactly 3 runs")
        try:
            triples = [(name, by_name[name].levels[level].keys_at_1pct)
                       for name in run_names]
        except KeyError as e:
            raise ValueError(f"venn references unknown run name: {e}")
        svg = plots.venn3_overlap(
            triples,
            title=spec.get("title", "overlap"),
            subtitle=spec.get("subtitle"),
        )
        blocks.append(f"<div class='plot'>{svg}</div>")
    return f"<div class='plots'>{''.join(blocks)}</div>"


def _deep_dive_block(cfg: dict, sections: list[_Section]) -> str:
    dd = cfg.get("deep_dive")
    if not dd:
        return ""
    by_name = {s.name: s for s in sections}
    if dd["run_name"] not in by_name:
        raise ValueError(f"deep_dive.run_name {dd['run_name']!r} not in [[runs]]")
    s = by_name[dd["run_name"]]

    truth = ground_truth.union(*[ground_truth.load(p) for p in s.truth_paths])

    variants = dd.get("variants") or [
        ("hyperscore", "psm"),
        ("hyperscore", "peptide_psm_peptide"),
        ("sage_discriminant_score", "psm"),
        ("sage_discriminant_score", "peptide_psm_peptide"),
    ]

    pieces = [
        f"<h2>TDC variant analysis · {html.escape(s.name)}</h2>"
        f"<p class='meta'>For each level, sage's native q-value is "
        f"compared against sagepy's TDC implementation re-run on the raw "
        f"hyperscore (no ML rescoring) and on sage's LDA discriminant "
        f"(post-rescoring). If the LDA is the source of the q-value drift, "
        f"hyperscore curves should sit on / below y = x while LDA curves "
        f"should lift above it.</p>"
    ]

    for level in LEVELS:
        # Build TDC ROC + surviving sets per (score, method) at this level.
        # At the protein level the sagepy Rust binding rejects
        # `picked_protein`; we use `peptide_peptide_only` over a
        # match_idx = accession instead (effectively picked protein).
        if level == "protein":
            level_variants = [
                ("hyperscore", "peptide_peptide_only"),
                ("sage_discriminant_score", "peptide_peptide_only"),
            ]
        else:
            level_variants = list(variants)

        tdc_rocs = []
        for score, method in level_variants:
            try:
                tdc_rocs.append(tdc_mod.run(s.sage_parquet, truth, score=score,
                                             method=method, level=level))
            except Exception as exc:  # pragma: no cover
                pieces.append(
                    f"<p class='meta'><em>{level}/{score}/{method} skipped: "
                    f"{html.escape(str(exc))}</em></p>"
                )
        native_curve = plots.RunCurve(name="sage native", roc=s.levels[level].roc)
        tdc_svg = plots.tdc_calibration(native_curve, tdc_rocs, log=True)

        # Three-way Venn: sage native vs hyperscore-TDC vs LDA-TDC at this level.
        try:
            if level == "protein":
                hs_method = "peptide_peptide_only"
            elif level == "peptide":
                hs_method = "peptide_psm_peptide"
            else:
                hs_method = "psm"
            hs_set = tdc_mod.surviving_keys(s.sage_parquet, "hyperscore", hs_method,
                                            level=level, q_cutoff=_VENN_Q)
            ld_set = tdc_mod.surviving_keys(s.sage_parquet, "sage_discriminant_score",
                                            hs_method, level=level, q_cutoff=_VENN_Q)
            venn_svg = plots.venn3_overlap(
                [
                    ("sage native", s.levels[level].keys_at_1pct),
                    ("hyperscore TDC", hs_set),
                    ("LDA TDC", ld_set),
                ],
                title=f"{level} overlap @ q ≤ 0.01",
                subtitle="exclusive lobes show estimator-specific calls",
            )
        except Exception as exc:  # pragma: no cover
            venn_svg = f"<p class='meta'>Venn unavailable: {html.escape(str(exc))}</p>"

        # q ≤ 0.01 table per level.
        rows = []
        nh = s.levels[level].headline
        rows.append(
            {
                "estimator": "sage native q",
                "keys @ q≤0.01": int(nh.get("reported_keys", 0)),
                "TP": int(nh.get("tp", 0)),
                "FP": int(nh.get("fp", 0)),
                "true FDR": f"{nh.get('true_fdr', 0.0) * 100:.2f} %",
                "TPR": f"{nh.get('tpr', 0.0) * 100:.2f} %",
            }
        )
        for tdc in tdc_rocs:
            sub = tdc.df[tdc.df["q_cutoff"] <= 0.01]
            if len(sub):
                last = sub.iloc[-1]
                rows.append(
                    {
                        "estimator": tdc.label,
                        "keys @ q≤0.01": int(last["reported_keys"]),
                        "TP": int(last["tp"]),
                        "FP": int(last["fp"]),
                        "true FDR": f"{last['true_fdr'] * 100:.2f} %",
                        "TPR": f"{last['tpr'] * 100:.2f} %",
                    }
                )
        tdc_table = pd.DataFrame(rows).to_html(index=False, classes="summary",
                                                border=0, escape=False)

        pieces.append(
            f"<details open><summary><strong>Level: {level}</strong></summary>"
            f"<div class='plots'>"
            f"<div class='plot'>{tdc_svg}</div>"
            f"<div class='plot'>{venn_svg}</div>"
            f"</div>"
            f"{tdc_table}"
            f"</details>"
        )

    return "\n".join(pieces)


_HTML_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
          max-width: 1100px; margin: 2rem auto; padding: 0 1rem; color: #222; }}
  h1 {{ font-size: 1.6rem; }}
  h2 {{ font-size: 1.2rem; margin-top: 2.2rem; border-bottom: 1px solid #eee; padding-bottom: 0.3rem; }}
  h3 {{ font-size: 1.0rem; margin-top: 1.5rem; }}
  table {{ border-collapse: collapse; margin: 0.5rem 0; font-size: 0.88rem; }}
  th, td {{ padding: 4px 12px; text-align: right; border-bottom: 1px solid #eee; }}
  th:first-child, td:first-child {{ text-align: left; }}
  table.summary th {{ background: #f6f8fa; }}
  .meta {{ color: #666; font-size: 0.85rem; }}
  .plots {{ display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; margin: 1rem 0; }}
  .plot {{ border: 1px solid #eee; padding: 0.5rem; background: #fff; }}
  .notes {{ background: #f9f9f9; padding: 0.7rem 1rem; border-left: 3px solid #ddd; }}
  details {{ margin: 0.5rem 0; padding: 0.6rem 0.9rem; border: 1px solid #eee; background: #fafbfc; }}
  details > summary {{ cursor: pointer; color: #333; font-size: 0.95rem; }}
  details[open] > summary {{ margin-bottom: 0.5rem; }}
  @media (max-width: 760px) {{ .plots {{ grid-template-columns: 1fr; }} }}
</style>
</head>
<body>

<h1>{title}</h1>
<p class="meta">{subtitle}</p>

<div class="notes">{notes_html}</div>

<h2>Level: ion (precursor)</h2>
{ion_block}

<h2>Level: peptide</h2>
{peptide_block}

<h2>Level: protein</h2>
{protein_block}

<h2>Overlap Venns @ q ≤ 0.01</h2>
<p class="meta">Each circle is the set of level-keys reported by one
run at 1 % FDR. Big intersections = consistent identifications across
runs; large exclusive lobes = a run is finding (or missing) things the
others aren't. Three sets of Venns — one per level.</p>
<h3>ion-level overlap</h3>
{ion_venns}
<h3>peptide-level overlap</h3>
{peptide_venns}
<h3>protein-level overlap</h3>
{protein_venns}

{deep_dive_html}

<h2>Per-run sweep tables</h2>
{per_run_html}

</body>
</html>
"""


def build(config_path: Path | str, out_path: Path | str) -> Path:
    cfg = tomllib.loads(Path(config_path).read_text())
    sections = [_load_section(spec) for spec in cfg["runs"]]
    venn_specs = cfg.get("venns") or []

    blocks = {level: _level_block(level, sections) for level in LEVELS}
    venns = (
        {level: _venn_block(level, sections, venn_specs) for level in LEVELS}
        if venn_specs
        else {level: "<p class='meta'>No Venns configured.</p>" for level in LEVELS}
    )

    per_run = []
    for s in sections:
        sub_blocks = []
        for level in LEVELS:
            v = s.levels[level]
            sub_blocks.append(
                f"<h4>{level}</h4>"
                f"<p class='meta'>truth: {v.n_truth:,} entries</p>"
                f"{_sweep_table(v)}"
            )
        per_run.append(
            f"<details><summary><strong>{html.escape(s.name)}</strong> · "
            f"{s.n_psms:,} PSMs ({s.n_targets:,} target / {s.n_decoys:,} decoy)"
            f"</summary>{''.join(sub_blocks)}</details>"
        )

    out = _HTML_TEMPLATE.format(
        title=html.escape(cfg.get("title", "SAGEBench report")),
        subtitle=html.escape(cfg.get("subtitle", "")),
        notes_html=_format_notes(cfg.get("notes")),
        ion_block=blocks["ion"],
        peptide_block=blocks["peptide"],
        protein_block=blocks["protein"],
        ion_venns=venns["ion"],
        peptide_venns=venns["peptide"],
        protein_venns=venns["protein"],
        deep_dive_html=_deep_dive_block(cfg, sections),
        per_run_html="\n".join(per_run),
    )
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(out, encoding="utf-8")
    return out_path
