"""Bypass sage's built-in q-values and re-rescore via sagepy's TDC.

This is the analysis that probes the q-value drift on simulated DDA:
is sage's reported FDR off because the LDA discriminant overfits, or
is the calibration error already present in the raw hyperscore?

We expose a small helper that, given a sage parquet and a ground-truth
set, runs `sagepy.qfdr.tdc.target_decoy_competition_pandas` with the
chosen score column and TDC method, then computes the same true-FDR
metrics our harness reports for sage's native q-values.

Three levels are supported via `level=`:

    ion     → match_idx = peptide+charge, default method = psm
    peptide → match_idx = peptide,        default method = peptide_psm_peptide
    protein → match_idx = accession,      default method = picked_protein
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from sagepy.qfdr.tdc import target_decoy_competition_pandas

from . import sage_output
from .ground_truth import GroundTruth, Level, _extract_accessions


VALID_METHODS = (
    "psm",
    "peptide_psm_only",
    "peptide_peptide_only",
    "peptide_psm_peptide",
)
VALID_SCORES = ("hyperscore", "sage_discriminant_score")
# At the protein level we set match_idx = accession and reuse
# `peptide_peptide_only` (which competes by match_idx) — this is
# effectively picked-protein TDC. The picked_protein method exists in
# the sagepy Python shim but the Rust binding rejects it as of v0.6.x.
DEFAULT_METHOD: dict[Level, str] = {
    "ion":     "psm",
    "peptide": "peptide_psm_peptide",
    "protein": "peptide_peptide_only",
}


def _tdc_input(parquet_path: str, level: Level = "ion") -> pd.DataFrame:
    """Build sagepy TDC's required-column frame from a sage parquet.

    Sagepy needs `spec_idx`, `match_idx`, `match_identity_candidates`,
    `decoy`, plus the score column. The shape of `match_idx` depends on
    the level we're analysing — peptide+charge for ion, peptide for
    peptide, accession for protein.
    """
    df = pq.read_table(parquet_path).to_pandas()

    pep_unimod = df["peptide"].astype(str).map(sage_output._peaks_to_unimod)
    charge = df["charge"].astype(int)

    if level == "ion":
        match_idx = pep_unimod + ":" + charge.astype(str)
    elif level == "peptide":
        match_idx = pep_unimod
    elif level == "protein":
        # First parsed accession; for decoys sage prefixes "rev_" so
        # target / decoy compete correctly when both share the same
        # underlying protein.
        def first_acc(s: str) -> str:
            accs = _extract_accessions(s)
            if accs:
                return accs[0]
            return s.split(";")[0].strip() or s
        match_idx = df["proteins"].astype(str).map(first_acc)
    else:
        raise ValueError(f"unknown level: {level!r}")

    out = pd.DataFrame(
        {
            "spec_idx": df["psm_id"].astype(str),
            "match_idx": match_idx.astype(str),
            "match_identity_candidates": [[m] for m in match_idx.astype(str)],
            "decoy": df["is_decoy"].astype(bool),
            "hyperscore": df["hyperscore"].astype(float),
            "sage_discriminant_score": df["sage_discriminant_score"].astype(float),
            "peptide": pep_unimod,
            "charge": charge,
        }
    )
    return out


@dataclass(frozen=True)
class TDCRoc:
    """ROC table for one (level, score, method) combination."""

    level: Level
    score: str
    method: str
    df: pd.DataFrame  # cols: q_cutoff, reported_keys, tp, fp, true_fdr, tpr

    @property
    def label(self) -> str:
        method_nice = {
            "psm": "PSM-level",
            "peptide_psm_only": "peptide (PSM-only)",
            "peptide_peptide_only": "peptide (direct)",
            "peptide_psm_peptide": "peptide (double)",
            "picked_peptide": "picked peptide",
            "picked_protein": "picked protein",
        }[self.method]
        score_nice = {
            "hyperscore": "hyperscore",
            "sage_discriminant_score": "LDA disc.",
        }[self.score]
        return f"{score_nice} · {method_nice}"


def run(
    sage_parquet: str,
    truth: GroundTruth,
    score: str = "hyperscore",
    method: Optional[str] = None,
    n_points: int = 250,
    level: Level = "ion",
) -> TDCRoc:
    """Run sagepy TDC and emit a true-FDR ROC table at the chosen level."""
    if score not in VALID_SCORES:
        raise ValueError(f"score must be one of {VALID_SCORES}, got {score!r}")
    if method is None:
        method = DEFAULT_METHOD[level]
    if method not in VALID_METHODS:
        raise ValueError(f"method must be one of {VALID_METHODS}, got {method!r}")

    tdc_in = _tdc_input(sage_parquet, level=level)
    tdc_out = target_decoy_competition_pandas(tdc_in, method=method, score=score)

    # Recover keys from match_idx for the join.
    if level == "ion":
        parts = tdc_out["match_idx"].str.rsplit(":", n=1, expand=True)
        peps = parts[0].map(sage_output._peaks_to_unimod)
        charges = parts[1].astype(int)
        targets = tdc_out[~tdc_out["decoy"].astype(bool)].copy()
        targets = targets.assign(peptide=peps[~tdc_out["decoy"].astype(bool)].values,
                                 charge=charges[~tdc_out["decoy"].astype(bool)].values)
    elif level == "peptide":
        peps = tdc_out["match_idx"].map(sage_output._peaks_to_unimod)
        targets = tdc_out[~tdc_out["decoy"].astype(bool)].copy()
        targets = targets.assign(peptide=peps[~tdc_out["decoy"].astype(bool)].values)
    elif level == "protein":
        targets = tdc_out[~tdc_out["decoy"].astype(bool)].copy()
        targets["protein"] = targets["match_idx"]
    else:
        raise ValueError(level)

    truth_keys = truth.keys_at(level)
    if level == "ion":
        targets["tp"] = [
            (p, c) in truth_keys for p, c in zip(targets["peptide"], targets["charge"])
        ]
        dedup = ["peptide", "charge"]
    elif level == "peptide":
        targets["tp"] = targets["peptide"].isin(truth_keys)
        dedup = ["peptide"]
    else:
        targets["tp"] = targets["protein"].isin(truth_keys)
        dedup = ["protein"]

    qs = np.linspace(0.0, 1.0, n_points + 1)[1:]
    if len(targets) == 0:
        cutoffs: list[float] = []
    else:
        cutoffs = sorted(set(np.quantile(targets["q_value"], qs).tolist()))

    rows = []
    n_truth = truth.n_at(level)
    for q in cutoffs:
        reported = targets[targets["q_value"] <= q]
        keys = reported.groupby(dedup, as_index=False)["tp"].any()
        tp = int(keys["tp"].sum())
        fp = int(len(keys) - tp)
        rows.append(
            {
                "q_cutoff": float(q),
                "reported_keys": len(keys),
                "reported_precursors": len(keys),
                "tp": tp,
                "fp": fp,
                "true_fdr": fp / (tp + fp) if (tp + fp) > 0 else 0.0,
                "tpr": tp / n_truth if n_truth > 0 else 0.0,
            }
        )
    return TDCRoc(level=level, score=score, method=method, df=pd.DataFrame(rows))


def surviving_keys(
    sage_parquet: str,
    score: str,
    method: str,
    level: Level,
    q_cutoff: float = 0.01,
) -> set:
    """Set of level-keys retained by TDC at q ≤ cutoff (for Venns)."""
    tdc_in = _tdc_input(sage_parquet, level=level)
    tdc_out = target_decoy_competition_pandas(tdc_in, method=method, score=score)
    keep = tdc_out[(~tdc_out["decoy"].astype(bool)) & (tdc_out["q_value"] <= q_cutoff)]
    if level == "ion":
        parts = keep["match_idx"].str.rsplit(":", n=1, expand=True)
        peps = parts[0].map(sage_output._peaks_to_unimod)
        charges = parts[1].astype(int)
        return set(zip(peps, charges, strict=True))
    if level == "peptide":
        return set(keep["match_idx"].map(sage_output._peaks_to_unimod))
    if level == "protein":
        return set(keep["match_idx"])
    raise ValueError(level)
