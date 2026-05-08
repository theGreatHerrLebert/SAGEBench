"""Extract a small TimSim from_findings seed CSV from an existing
synthetic_data.db.

Used to build seeds/hela_ci_smoke.csv from the canonical HeLa 150K G07M
run (SUBMISSION/MBR/data/150K_G07M/TIMSIM-DDA-HELA-150K-G07M-001/
synthetic_data.db). The output matches the 6-column from_findings
contract:

    protein, sequence, intensity, charge, rt, im

`rt` and `im` are simply omitted (the columns aren't emitted at all)
so TimSim re-simulates them — that keeps the smoke independent of the
source gradient and reference dataset's IM range. We deliberately
do NOT emit empty-string rt/im columns: TimSim's `_read_findings`
runs `dropna()` across every optional column it sees, so a present-
but-empty column wipes every row.

Run:
    python scripts/extract_smoke_seed.py \
        --src SUBMISSION/MBR/data/150K_G07M/TIMSIM-DDA-HELA-150K-G07M-001/synthetic_data.db \
        --out seeds/hela_ci_smoke.csv \
        --n 1500 \
        --seed 41
"""
from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path

import pandas as pd


def extract_seed(src: Path, n: int, seed: int) -> pd.DataFrame:
    with sqlite3.connect(src) as cx:
        ions = pd.read_sql_query(
            """
            SELECT i.peptide_id, i.sequence, i.charge, i.relative_abundance
            FROM ions AS i
            WHERE i.charge BETWEEN 2 AND 4
            """,
            cx,
        )
        peptides = pd.read_sql_query(
            "SELECT peptide_id, protein, events FROM peptides WHERE decoy = 0",
            cx,
        )

    merged = ions.merge(peptides, on="peptide_id", how="inner")
    # Per-peptide intensity proxy: events * relative_abundance, scaled so
    # the median is around 1e5 (in the range TimSim's auto-scaler is
    # comfortable with for intensity_multiplier ~ 1.0).
    merged["intensity"] = merged["events"] * merged["relative_abundance"]
    median = merged["intensity"].median()
    if median > 0:
        merged["intensity"] = merged["intensity"] * (1e5 / median)

    # Pick top-charge per peptide so each peptide contributes one row,
    # then sample down to n.
    one_per_pep = (
        merged.sort_values("intensity", ascending=False)
        .drop_duplicates("peptide_id")
        .reset_index(drop=True)
    )
    sampled = one_per_pep.sample(n=min(n, len(one_per_pep)), random_state=seed)
    sampled = sampled.sort_values("peptide_id").reset_index(drop=True)

    out = pd.DataFrame(
        {
            "protein": sampled["protein"].fillna(""),
            "sequence": sampled["sequence"],
            "intensity": sampled["intensity"].astype(float).round(2),
            "charge": sampled["charge"].astype(int),
        }
    )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--n", type=int, default=1500)
    ap.add_argument("--seed", type=int, default=41)
    args = ap.parse_args()

    out = extract_seed(args.src, args.n, args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    n_mod = out["sequence"].str.contains(r"\[UNIMOD:\d+\]", regex=True).sum()
    print(
        f"wrote {args.out} ({len(out):,} rows; "
        f"{out['sequence'].nunique():,} unique sequences; "
        f"{n_mod:,} carry at least one UniMod tag)"
    )


if __name__ == "__main__":
    main()
