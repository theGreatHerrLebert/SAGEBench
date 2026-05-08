"""Load TimSim ground truth from a synthetic_data.db file.

Three views over the same simulated peptides:

- `ion`     : (sequence_unimod, charge) precursors  (one row per distinct ion)
- `peptide` : sequence_unimod                       (one row per modified-peptide form)
- `protein` : protein accession                     (one row per UniProt accession)

Decoys (`peptides.decoy = 1`) are excluded — they exist in the seed
but are not real signal.
"""
from __future__ import annotations

import re
import sqlite3
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import pandas as pd

Level = Literal["ion", "peptide", "protein"]
LEVELS: tuple[Level, ...] = ("ion", "peptide", "protein")

# Match `sp|P12345|...` or `tr|Q9...|...` headers; capture the accession.
_ACCESSION_RE = re.compile(r"\b(?:sp|tr)\|([A-Z0-9-]+)\|")


def _extract_accessions(protein_field: str) -> list[str]:
    """Pull accession IDs out of TimSim/sage-style protein descriptors.

    TimSim records the full FASTA header (`sp|P12345|GENE_HUMAN ...`)
    and may join multiple accessions with `:` for shared peptides; sage
    joins with `;`. We accept both. Falls back to the raw string when no
    structured accession is found (e.g. trivial test-only inputs).
    """
    if not isinstance(protein_field, str) or not protein_field:
        return []
    matches = _ACCESSION_RE.findall(protein_field)
    if matches:
        return matches
    # No sp|...| headers — split on either separator and return the
    # bare strings (for cases like "P12345;Q67890" or test data).
    parts = [p.strip() for p in re.split(r"[;:]", protein_field) if p.strip()]
    return parts or []


@dataclass(frozen=True)
class GroundTruth:
    """Three-level view over a TimSim synthetic_data.db's target peptides.

    Each level carries two denominators:

    - **full truth** — every entry the simulator injected
    - **fragmentable subset** — entries that were actually picked by DDA
      topN (i.e. appear in `pasef_meta` via the `precursors` table). For
      DDA the fragmentable subset is much smaller than the full truth;
      reporting TPR over both is honest about what's identifiable in
      principle vs what got an MS2 spectrum.
    """

    ions: pd.DataFrame      # cols: sequence, charge, peptide_id, protein, intensity_proxy
    peptides: pd.DataFrame  # cols: sequence, peptide_id, protein
    proteins: pd.DataFrame  # cols: accession (one per row, deduped)

    fragmentable_ion_keys: set = field(default_factory=set)
    fragmentable_peptide_keys: set = field(default_factory=set)
    fragmentable_protein_keys: set = field(default_factory=set)

    # ----- counts -----

    @property
    def n_ions(self) -> int:
        return len(self.ions)

    @property
    def n_peptides(self) -> int:
        return len(self.peptides)

    @property
    def n_proteins(self) -> int:
        return len(self.proteins)

    @property
    def n_targets(self) -> int:  # backwards compat
        return self.n_ions

    # ----- key sets -----

    @property
    def ion_keys(self) -> set[tuple[str, int]]:
        return set(zip(self.ions["sequence"], self.ions["charge"], strict=True))

    @property
    def peptide_keys(self) -> set[str]:
        return set(self.peptides["sequence"])

    @property
    def protein_keys(self) -> set[str]:
        return set(self.proteins["accession"])

    @property
    def keys(self) -> set[tuple[str, int]]:  # backwards compat
        return self.ion_keys

    def keys_at(self, level: Level) -> set:
        if level == "ion":
            return self.ion_keys
        if level == "peptide":
            return self.peptide_keys
        if level == "protein":
            return self.protein_keys
        raise ValueError(f"unknown level: {level!r}")

    def n_at(self, level: Level) -> int:
        if level == "ion":
            return self.n_ions
        if level == "peptide":
            return self.n_peptides
        if level == "protein":
            return self.n_proteins
        raise ValueError(f"unknown level: {level!r}")

    def fragmentable_at(self, level: Level) -> set:
        if level == "ion":
            return self.fragmentable_ion_keys
        if level == "peptide":
            return self.fragmentable_peptide_keys
        if level == "protein":
            return self.fragmentable_protein_keys
        raise ValueError(level)

    def n_fragmentable_at(self, level: Level) -> int:
        return len(self.fragmentable_at(level))


def _fragmentable_ion_keys(
    cx: sqlite3.Connection,
    ions: pd.DataFrame,
    *,
    mz_tol_da: float = 0.005,
) -> set[tuple[str, int]]:
    """Ions whose (m/z, charge) match a precursor that fired in pasef_meta.

    pasef_meta.precursor → precursors.id; we keep precursors that actually
    fired (= every row, since precursors is only populated for selected
    ones in TimSim). For each (charge, monoisotopic_mz) we collect all
    truth ions within `mz_tol_da` — that's the chimeric expansion of the
    isolation window.
    """
    pre = pd.read_sql_query(
        """
        SELECT DISTINCT p.id, p.charge, p.monoisotopic_mz
        FROM precursors p
        JOIN pasef_meta pm ON pm.precursor = p.id
        """,
        cx,
    )
    if len(pre) == 0 or len(ions) == 0:
        return set()

    # Bucket ions by (charge, mz_bucket) for fast membership test.
    bucket_size = mz_tol_da
    ions = ions.copy()
    ions["mz_bucket"] = (ions["mz"] / bucket_size).round().astype(int)
    by_key: dict[tuple[int, int], list[tuple[float, str]]] = {}
    for c, b, m, s in zip(ions["charge"], ions["mz_bucket"], ions["mz"],
                          ions["sequence"]):
        by_key.setdefault((int(c), int(b)), []).append((float(m), str(s)))

    frag: set[tuple[str, int]] = set()
    for c, mz in zip(pre["charge"], pre["monoisotopic_mz"]):
        bucket = int(round(mz / bucket_size))
        for b in (bucket - 1, bucket, bucket + 1):
            for cand_mz, cand_seq in by_key.get((int(c), b), ()):
                if abs(cand_mz - mz) <= mz_tol_da:
                    frag.add((cand_seq, int(c)))
    return frag


def load(db_path: Path | str) -> GroundTruth:
    """Read peptides + ions from a TimSim synthetic_data.db."""
    db_path = Path(db_path)
    if not db_path.exists():
        raise FileNotFoundError(db_path)

    with sqlite3.connect(db_path) as cx:
        ions_raw = pd.read_sql_query(
            "SELECT peptide_id, sequence, charge, mz, relative_abundance FROM ions",
            cx,
        )
        peptides_raw = pd.read_sql_query(
            """
            SELECT peptide_id, sequence, protein, events
            FROM peptides
            WHERE decoy = 0
            """,
            cx,
        )

        # Ion view: one row per (sequence, charge) precursor.
        ions = ions_raw.merge(peptides_raw, on="peptide_id", how="inner",
                              suffixes=("", "_pep"))
        ions["intensity_proxy"] = ions["events"] * ions["relative_abundance"]
        ions_view = (
            ions[["sequence", "charge", "peptide_id", "protein", "intensity_proxy"]]
            .drop_duplicates(subset=["sequence", "charge"])
            .reset_index(drop=True)
        )

        # Compute the fragmentable ion set BEFORE we drop the m/z column.
        frag_ions = _fragmentable_ion_keys(cx, ions[["sequence", "charge", "mz"]])

    # Peptide view: collapse charges. One row per modified-peptide form.
    peptides = (
        peptides_raw[["peptide_id", "sequence", "protein"]]
        .drop_duplicates(subset=["sequence"])
        .reset_index(drop=True)
    )

    # Protein view: explode TimSim's protein descriptor into accessions.
    accs: set[str] = set()
    pep_to_accs: dict[str, set[str]] = {}
    for seq, raw in zip(peptides_raw["sequence"], peptides_raw["protein"]):
        if not isinstance(raw, str):
            continue
        local = set(_extract_accessions(raw))
        accs.update(local)
        pep_to_accs.setdefault(seq, set()).update(local)
    proteins = pd.DataFrame({"accession": sorted(accs)})

    # Roll fragmentable ions up to peptide / protein.
    frag_peps = {seq for seq, _ in frag_ions}
    frag_prots: set[str] = set()
    for seq in frag_peps:
        frag_prots.update(pep_to_accs.get(seq, ()))

    return GroundTruth(
        ions=ions_view,
        peptides=peptides,
        proteins=proteins,
        fragmentable_ion_keys=frag_ions,
        fragmentable_peptide_keys=frag_peps,
        fragmentable_protein_keys=frag_prots,
    )


def union(*truths: GroundTruth) -> GroundTruth:
    """Combine ground truths from multiple replicates (multi-file runs)."""
    if not truths:
        raise ValueError("union(...) requires at least one GroundTruth")
    ions = pd.concat([t.ions for t in truths], ignore_index=True) \
        .drop_duplicates(subset=["sequence", "charge"]).reset_index(drop=True)
    peptides = pd.concat([t.peptides for t in truths], ignore_index=True) \
        .drop_duplicates(subset=["sequence"]).reset_index(drop=True)
    proteins = pd.concat([t.proteins for t in truths], ignore_index=True) \
        .drop_duplicates(subset=["accession"]).reset_index(drop=True)
    frag_ions = set().union(*(t.fragmentable_ion_keys for t in truths))
    frag_peps = set().union(*(t.fragmentable_peptide_keys for t in truths))
    frag_prots = set().union(*(t.fragmentable_protein_keys for t in truths))
    return GroundTruth(
        ions=ions, peptides=peptides, proteins=proteins,
        fragmentable_ion_keys=frag_ions,
        fragmentable_peptide_keys=frag_peps,
        fragmentable_protein_keys=frag_prots,
    )
