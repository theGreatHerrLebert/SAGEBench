"""Parse Sage's PSM parquet output into a normalized DataFrame.

We pull only the columns the harness needs:

    peptide          - peptide sequence with inline UniMod tags
    charge           - precursor charge
    proteins         - semicolon-joined protein accessions
    sage_discriminant_score / spectrum_q
    is_decoy         - 1 if decoy hit, 0 otherwise

Sage's column names have shifted across versions (`peptide` vs.
`stripped_peptide`, `is_decoy` vs. `label`); the loader is tolerant.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


# Sage emits modifications as bracketed +mass (e.g. "M[+15.9949]" or
# "C[+57.0215]"). Map the few masses we actually use back to UniMod
# tags so we can join against TimSim's UniMod-encoded sequences.
# We match a couple of common rounding variants per mass.
_MASS_TO_UNIMOD = {
    "57.0215":  "[UNIMOD:4]",    # Carbamidomethyl(C)
    "57.02146": "[UNIMOD:4]",
    "15.9949":  "[UNIMOD:35]",   # Oxidation(M)
    "15.99491": "[UNIMOD:35]",
    "42.0106":  "[UNIMOD:1]",    # Acetyl(N-term)
    "42.01057": "[UNIMOD:1]",
    "79.9663":  "[UNIMOD:21]",   # Phospho
    "79.96633": "[UNIMOD:21]",
}

# Match "[+15.9949]" — sage's bracketed mass-delta form.
_MASS_RE = re.compile(r"\[\+(\d+\.\d+)\]")
# Sage decorates terminal mods with a hyphen ("[+42.01]-MGAV..." or
# "...PEPTIDE-[+42.01]"); TimSim does not. Drop the hyphen post-translation.
_NTERM_HYPHEN_RE = re.compile(r"^(\[UNIMOD:\d+\])-")
_CTERM_HYPHEN_RE = re.compile(r"-(\[UNIMOD:\d+\])$")


def _peaks_to_unimod(seq: str) -> str:
    """Translate Sage-style ``[+mass]`` tags to inline UniMod tags and
    normalise terminal-modification punctuation to TimSim's convention."""

    def repl(m: re.Match[str]) -> str:
        mass = m.group(1)
        return _MASS_TO_UNIMOD.get(mass, m.group(0))

    seq = _MASS_RE.sub(repl, seq)
    seq = _NTERM_HYPHEN_RE.sub(r"\1", seq)
    seq = _CTERM_HYPHEN_RE.sub(r"\1", seq)
    return seq


@dataclass(frozen=True)
class SageHits:
    df: pd.DataFrame  # peptide, charge, proteins, score, q_value, is_decoy

    @property
    def n_psms(self) -> int:
        return len(self.df)

    @property
    def n_targets(self) -> int:
        return int((self.df["is_decoy"] == 0).sum())


def load(parquet_path: Path | str) -> SageHits:
    """Load a Sage `results.sage.parquet` (or equivalent) PSM table."""
    parquet_path = Path(parquet_path)
    df = pd.read_parquet(parquet_path)
    cols = {c.lower(): c for c in df.columns}

    pep_col = cols.get("peptide") or cols.get("stripped_peptide")
    if pep_col is None:
        raise ValueError(f"no peptide column in {parquet_path} (cols: {list(df.columns)})")

    charge_col = cols.get("charge") or cols.get("precursor_charge")
    if charge_col is None:
        raise ValueError(f"no charge column in {parquet_path}")

    proteins_col = cols.get("proteins") or cols.get("protein")
    score_col = (
        cols.get("sage_discriminant_score")
        or cols.get("discriminant_score")
        or cols.get("score")
    )
    # Default to peptide_q so a `q ≤ 0.01` cutoff means 1% peptide-FDR,
    # which is the convention TimSim's ground truth is defined against.
    # spectrum_q (PSM-level) and protein_q (protein-level) are kept on
    # the loaded frame so callers can override via the q_column kw.
    q_col = (
        cols.get("peptide_q")
        or cols.get("spectrum_q")
        or cols.get("q_value")
    )

    if "is_decoy" in cols:
        decoy = df[cols["is_decoy"]].astype(int)
    elif "label" in cols:
        decoy = (df[cols["label"]] != 1).astype(int)  # sage label: 1=target, -1=decoy
    else:
        raise ValueError(f"no decoy column in {parquet_path}")

    out = pd.DataFrame(
        {
            "peptide": df[pep_col].astype(str).map(_peaks_to_unimod),
            "charge": df[charge_col].astype(int),
            "proteins": df[proteins_col].astype(str) if proteins_col else "",
            "score": df[score_col].astype(float) if score_col else 0.0,
            "q_value": df[q_col].astype(float) if q_col else None,
            "is_decoy": decoy.astype(int),
        }
    )
    # Primary protein accession for protein-level joins (first non-decoy
    # accession parsed from the proteins field). Decoys get the rev_-tagged
    # form retained so target/decoy stay distinct in TDC.
    from .ground_truth import _extract_accessions
    out["protein_accession"] = out["proteins"].map(
        lambda s: (_extract_accessions(s) or [s])[0] if isinstance(s, str) else ""
    )

    # Carry the unfiltered q-columns through so callers can rescore at a
    # different FDR level (e.g. q_column="spectrum_q" for PSM-level).
    for level in ("spectrum_q", "peptide_q", "protein_q", "protein_group_q"):
        if level in cols:
            out[level] = df[cols[level]].astype(float)
    return SageHits(df=out.reset_index(drop=True))
