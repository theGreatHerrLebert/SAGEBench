"""Inline-SVG matplotlib plots for the SAGEBench HTML report.

Two plot kinds:

- `calibration(...)` — true FDR (x) vs reported q-value (y) for each run.
  A perfectly calibrated search engine sits on y = x. Above the
  diagonal = optimistic (claims tighter FDR than reality); below =
  conservative.

- `coverage(...)` — true FDR (x) vs cumulative TP precursors (y).
  Pareto-style "how many real IDs at what real FDR" view.

Both functions return a string of SVG markup, ready to drop into HTML.
"""
from __future__ import annotations

import io
from dataclasses import dataclass

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

try:
    from matplotlib_venn import venn3, venn3_circles  # noqa: E402
    _HAS_VENN = True
except ImportError:  # pragma: no cover
    _HAS_VENN = False


@dataclass
class RunCurve:
    """ROC-style curve for a single run: q-cutoffs paired with the true
    FDR, TP precursor count, and TPR observed at each cutoff.
    """

    name: str
    roc: pd.DataFrame  # columns: q_cutoff, true_fdr, tp, tpr (from metrics.roc_points)


def _save_svg(fig: plt.Figure) -> str:
    buf = io.StringIO()
    fig.savefig(buf, format="svg", bbox_inches="tight")
    plt.close(fig)
    raw = buf.getvalue()
    # Strip the XML declaration and the DOCTYPE so the SVG is safe to
    # embed inline inside <body>.
    if raw.startswith("<?xml"):
        raw = raw.split("?>", 1)[1].lstrip()
    if raw.startswith("<!DOCTYPE"):
        raw = raw.split(">", 1)[1].lstrip()
    return raw


def calibration(curves: list[RunCurve], log: bool = True) -> str:
    """Calibration plot: true FDR (x) vs reported q-value (y).

    A run plotted *above* y=x is over-confident (sage says 1 % FDR
    while reality is 4 %). Below y=x is conservative.
    """
    fig, ax = plt.subplots(figsize=(5.5, 5.0), constrained_layout=True)

    for c in curves:
        df = c.roc.sort_values("q_cutoff")
        ax.plot(df["true_fdr"], df["q_cutoff"], label=c.name, lw=1.5)

    lo, hi = (1e-4, 0.5) if log else (0.0, 0.2)
    ax.plot([lo, hi], [lo, hi], color="0.4", ls="--", lw=1.0, label="y = x (calibrated)")

    if log:
        ax.set_xscale("log")
        ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("true FDR  (= FP / (TP + FP) at the precursor level)")
    ax.set_ylabel("reported q-value cutoff")
    ax.set_title("Sage q-value calibration vs simulated ground truth")
    ax.grid(alpha=0.3, which="both")
    ax.legend(loc="lower right", fontsize=8, frameon=False)

    return _save_svg(fig)


def tdc_calibration(
    sage_native: RunCurve,
    tdc_curves: "list",
    log: bool = True,
) -> str:
    """Calibration plot for a single dataset across multiple FDR estimators.

    `sage_native` is the baseline (sage's reported peptide_q). `tdc_curves`
    is a list of `sagebench.tdc.TDCRoc` objects — one per (score, method)
    combination.
    """
    fig, ax = plt.subplots(figsize=(5.5, 5.0), constrained_layout=True)

    df = sage_native.roc.sort_values("q_cutoff")
    ax.plot(df["true_fdr"], df["q_cutoff"], label="sage peptide_q (native)",
            lw=2.0, color="black")

    palette = ["#1f77b4", "#aec7e8", "#d62728", "#ff9896"]
    for tdc, color in zip(tdc_curves, palette, strict=False):
        df = tdc.df.sort_values("q_cutoff")
        ax.plot(df["true_fdr"], df["q_cutoff"], label=tdc.label, lw=1.4, color=color)

    lo, hi = (1e-4, 0.5) if log else (0.0, 0.2)
    ax.plot([lo, hi], [lo, hi], color="0.4", ls="--", lw=1.0, label="y = x (calibrated)")

    if log:
        ax.set_xscale("log")
        ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("true FDR (precursor level, vs TimSim truth)")
    ax.set_ylabel("reported q-value")
    ax.set_title("TDC variant comparison: where does the q-value drift come from?")
    ax.grid(alpha=0.3, which="both")
    ax.legend(loc="lower right", fontsize=8, frameon=False)

    return _save_svg(fig)


def venn3_overlap(
    sets: list[tuple[str, set]],
    title: str,
    subtitle: str | None = None,
) -> str:
    """3-way Venn over precursor (peptide, charge) sets.

    `sets` must be exactly three (label, set-of-keys) pairs. Each set is
    typically the (peptide, charge) precursors a given run/estimator
    reports at q ≤ 0.01.
    """
    if not _HAS_VENN:
        return f"<p class='meta'><em>matplotlib_venn not installed — Venn skipped: {title}</em></p>"
    if len(sets) != 3:
        raise ValueError(f"venn3_overlap needs exactly 3 sets, got {len(sets)}")

    labels = [s[0] for s in sets]
    sets_only = [s[1] for s in sets]

    fig, ax = plt.subplots(figsize=(5.5, 5.0), constrained_layout=True)
    venn3(sets_only, set_labels=labels, ax=ax,
          set_colors=("#1f77b4", "#2ca02c", "#d62728"), alpha=0.55)
    venn3_circles(sets_only, ax=ax, lw=0.7, color="0.4")
    ax.set_title(title, fontsize=11)
    if subtitle:
        ax.text(0.5, -0.06, subtitle, transform=ax.transAxes,
                ha="center", va="top", fontsize=9, color="#555")
    return _save_svg(fig)


def coverage(curves: list[RunCurve], y: str = "tp") -> str:
    """Coverage plot: true FDR (x) vs TP keys / TPR / fragmentable-TPR.

    Args:
        y: one of "tp" (count), "tpr" (vs full truth),
           "tpr_fragmentable" (vs DDA-selected subset).
    """
    fig, ax = plt.subplots(figsize=(5.5, 5.0), constrained_layout=True)

    for c in curves:
        df = c.roc.sort_values("true_fdr")
        if y in df.columns:
            ax.plot(df["true_fdr"], df[y], label=c.name, lw=1.5)

    ax.axvline(0.01, color="0.4", ls="--", lw=1.0)
    ax.set_xlabel("true FDR")
    ax.set_ylabel({
        "tp": "true-positive keys",
        "tpr": "TPR vs full truth",
        "tpr_fragmentable": "TPR vs fragmentable subset",
    }[y])
    ax.set_title({
        "tp": "True-positive keys vs true FDR",
        "tpr": "TPR (full truth) vs true FDR",
        "tpr_fragmentable": "TPR (fragmentable subset) vs true FDR",
    }[y])
    ax.set_xlim(0, 0.1)
    if y in ("tpr", "tpr_fragmentable"):
        ax.set_ylim(0, 1.0)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=8, frameon=False)

    return _save_svg(fig)
