"""True-FDR / TPR metrics at three levels (ion, peptide, protein).

Sage emits per-level q-values (`spectrum_q`, `peptide_q`, `protein_q`).
For each level we threshold on the appropriate q column, deduplicate
sage's target hits to one row per level-key, and join against the
matching ground-truth set (precursor pairs, peptide forms, or
protein accessions).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

from .ground_truth import GroundTruth, Level
from .sage_output import SageHits


# Map level → (sage q column, dedup keys).
_LEVEL_SPEC: dict[Level, tuple[str, list[str]]] = {
    "ion":     ("peptide_q", ["peptide", "charge"]),
    "peptide": ("peptide_q", ["peptide"]),
    "protein": ("protein_q", ["protein_accession"]),
}


@dataclass(frozen=True)
class EvalResult:
    """Counts at a single q-value cutoff, aggregated to the chosen level."""

    level: Level
    q_cutoff: float
    reported_psms: int
    reported_keys: int  # distinct level-keys reported (precursors / peptides / proteins)
    tp: int
    fp: int
    n_truth: int
    n_fragmentable: int
    true_fdr: float
    tpr: float
    tpr_fragmentable: float  # TP / |fragmentable|; the "honest" recall denominator

    # Backwards-compat alias used by older code paths / tests.
    @property
    def reported_precursors(self) -> int:
        return self.reported_keys

    def as_dict(self) -> dict:
        return {
            "level": self.level,
            "q_cutoff": self.q_cutoff,
            "reported_psms": self.reported_psms,
            "reported_keys": self.reported_keys,
            "tp": self.tp,
            "fp": self.fp,
            "n_truth": self.n_truth,
            "n_fragmentable": self.n_fragmentable,
            "true_fdr": self.true_fdr,
            "tpr": self.tpr,
            "tpr_fragmentable": self.tpr_fragmentable,
        }


def _q_column(hits: SageHits, level: Level) -> str:
    pref, _ = _LEVEL_SPEC[level]
    if pref in hits.df.columns:
        return pref
    # Fall back to the default "q_value" column the loader sets up
    # (which is peptide_q in practice). Protein-level then degrades to
    # peptide-FDR thresholding — flag that to the caller via the column
    # name they get back.
    return "q_value"


def _annotate(hits: SageHits, truth: GroundTruth, level: Level) -> pd.DataFrame:
    """Filter to target PSMs, attach a `tp` column based on level keys."""
    df = hits.df[hits.df["is_decoy"] == 0].copy()

    truth_keys = truth.keys_at(level)
    if level == "ion":
        keys = list(zip(df["peptide"], df["charge"], strict=True))
    elif level == "peptide":
        keys = list(df["peptide"])
    elif level == "protein":
        keys = list(df["protein_accession"])
    else:
        raise ValueError(level)
    df["tp"] = [k in truth_keys for k in keys]
    df["q_for_level"] = df[_q_column(hits, level)]
    return df


def evaluate(
    hits: SageHits,
    truth: GroundTruth,
    q_cutoff: float = 0.01,
    level: Level = "ion",
) -> EvalResult:
    """Score Sage hits at a single q-value cutoff, aggregated to `level`."""
    annotated = _annotate(hits, truth, level)
    if annotated["q_for_level"].isna().all():
        raise ValueError(f"Sage hits have no q-values for level {level!r}")

    reported = annotated[annotated["q_for_level"] <= q_cutoff]
    _, dedup_cols = _LEVEL_SPEC[level]
    keys = reported.groupby(dedup_cols, as_index=False)["tp"].any()
    tp = int(keys["tp"].sum())
    fp = int(len(keys) - tp)
    n_truth = truth.n_at(level)
    n_frag = truth.n_fragmentable_at(level)
    true_fdr = fp / (tp + fp) if (tp + fp) > 0 else 0.0
    tpr = tp / n_truth if n_truth > 0 else 0.0
    tpr_frag = tp / n_frag if n_frag > 0 else 0.0

    return EvalResult(
        level=level,
        q_cutoff=q_cutoff,
        reported_psms=len(reported),
        reported_keys=len(keys),
        tp=tp,
        fp=fp,
        n_truth=n_truth,
        n_fragmentable=n_frag,
        true_fdr=true_fdr,
        tpr=tpr,
        tpr_fragmentable=tpr_frag,
    )


def sweep(
    hits: SageHits,
    truth: GroundTruth,
    q_cutoffs: list[float] | None = None,
    level: Level = "ion",
) -> pd.DataFrame:
    if q_cutoffs is None:
        q_cutoffs = [0.001, 0.005, 0.01, 0.05]
    rows = [evaluate(hits, truth, q, level=level).as_dict() for q in q_cutoffs]
    df = pd.DataFrame(rows)
    # Backwards-compat alias column.
    df["reported_precursors"] = df["reported_keys"]
    return df


def roc_points(
    hits: SageHits,
    truth: GroundTruth,
    n_points: int = 200,
    level: Level = "ion",
) -> pd.DataFrame:
    annotated = _annotate(hits, truth, level)
    if annotated["q_for_level"].isna().all():
        raise ValueError(f"Sage hits have no q-values for level {level!r}")

    qs = np.linspace(0.0, 1.0, n_points + 1)[1:]
    cutoffs = np.quantile(annotated["q_for_level"], qs)
    n_truth = truth.n_at(level)
    n_frag = truth.n_fragmentable_at(level)
    _, dedup_cols = _LEVEL_SPEC[level]
    rows = []
    for q in cutoffs:
        reported = annotated[annotated["q_for_level"] <= q]
        keys = reported.groupby(dedup_cols, as_index=False)["tp"].any()
        tp = int(keys["tp"].sum())
        fp = int(len(keys) - tp)
        rows.append(
            {
                "q_cutoff": float(q),
                "reported_psms": len(reported),
                "reported_keys": len(keys),
                "reported_precursors": len(keys),  # alias
                "tp": tp,
                "fp": fp,
                "true_fdr": fp / (tp + fp) if (tp + fp) > 0 else 0.0,
                "tpr": tp / n_truth if n_truth > 0 else 0.0,
                "tpr_fragmentable": tp / n_frag if n_frag > 0 else 0.0,
            }
        )
    return pd.DataFrame(rows)


def precursor_keys_at(
    hits: SageHits,
    q_cutoff: float = 0.01,
    level: Level = "ion",
    q_col_override: str | None = None,
) -> set:
    """Set of level-keys for target PSMs at q ≤ q_cutoff (for Venns)."""
    df = hits.df[hits.df["is_decoy"] == 0]
    q_col = q_col_override or _q_column(hits, level)
    df = df[df[q_col] <= q_cutoff]
    if level == "ion":
        return set(zip(df["peptide"], df["charge"], strict=True))
    if level == "peptide":
        return set(df["peptide"])
    if level == "protein":
        return set(df["protein_accession"])
    raise ValueError(level)
