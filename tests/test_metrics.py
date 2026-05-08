"""Round-trip test for the SAGEBench harness:

1. Pull a small ground-truth set out of the canonical HeLa 150K
   synthetic_data.db.
2. Fabricate a fake Sage parquet with a known mix of TPs and FPs.
3. Confirm `metrics.evaluate(...)` recovers the expected counts.
"""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd
import pytest

from sagebench import ground_truth, metrics, sage_output


HELA_TRUTH_DB = Path(
    "/scratch/timsim-demo/SUBMISSION/MBR/data/150K_G07M"
    "/TIMSIM-DDA-HELA-150K-G07M-001/synthetic_data.db"
)


@pytest.fixture(scope="module")
def truth() -> ground_truth.GroundTruth:
    if not HELA_TRUTH_DB.exists():
        pytest.skip(f"missing {HELA_TRUTH_DB}")
    return ground_truth.load(HELA_TRUTH_DB)


def _truth_keys_sample(truth: ground_truth.GroundTruth, k: int) -> list[tuple[str, int]]:
    return list(truth.keys)[:k]


def test_round_trip_fdr_tpr(tmp_path, truth):
    """Sage hits = 80 real + 20 fake → expect TP=80, FP=20, true_fdr=0.2."""
    real = _truth_keys_sample(truth, 80)
    fake = [(f"FAKEPEP{i}KR", 2) for i in range(20)]

    rows = []
    for seq, charge in real + fake:
        rows.append(
            {
                "peptide": seq,
                "charge": charge,
                "proteins": "P00000",
                "sage_discriminant_score": 1.0,
                "spectrum_q": 0.005,  # all under 0.01
                "is_decoy": 0,
            }
        )
    # A handful of decoys (these never count as TP or FP at q<=0.01,
    # since they're filtered out before the join.)
    for i in range(5):
        rows.append(
            {
                "peptide": f"DECOYSEQ{i}",
                "charge": 2,
                "proteins": "rev_P00000",
                "sage_discriminant_score": -0.5,
                "spectrum_q": 0.5,
                "is_decoy": 1,
            }
        )

    pq = tmp_path / "fake_sage.parquet"
    pd.DataFrame(rows).to_parquet(pq)

    hits = sage_output.load(pq)
    assert hits.n_psms == 105
    assert hits.n_targets == 100

    result = metrics.evaluate(hits, truth, q_cutoff=0.01)
    assert result.reported_psms == 100
    assert result.reported_precursors == 100
    assert result.tp == 80
    assert result.fp == 20
    assert abs(result.true_fdr - 0.20) < 1e-9
    assert result.tpr == pytest.approx(80 / truth.n_targets, abs=1e-12)


def test_sweep_returns_one_row_per_cutoff(tmp_path, truth):
    real = _truth_keys_sample(truth, 50)
    rows = [
        {
            "peptide": seq,
            "charge": charge,
            "proteins": "P00000",
            "sage_discriminant_score": 1.0,
            "spectrum_q": 0.001,
            "is_decoy": 0,
        }
        for seq, charge in real
    ]
    pq = tmp_path / "small_sage.parquet"
    pd.DataFrame(rows).to_parquet(pq)

    hits = sage_output.load(pq)
    table = metrics.sweep(hits, truth, q_cutoffs=[0.001, 0.01, 0.05])
    assert list(table["q_cutoff"]) == [0.001, 0.01, 0.05]
    assert (table["reported_psms"] == 50).all()
    assert (table["reported_precursors"] == 50).all()
    assert (table["true_fdr"] == 0.0).all()


def test_mass_to_unimod_translation(tmp_path, truth):
    """Sage emits 'C+57.02146' style mods; harness must translate to UniMod."""
    seq, charge = next(
        (s, c) for s, c in truth.keys if "[UNIMOD:4]" in s
    )
    sage_seq = seq.replace("[UNIMOD:4]", "[+57.0215]")
    rows = [
        {
            "peptide": sage_seq,
            "charge": charge,
            "proteins": "P00000",
            "sage_discriminant_score": 1.0,
            "spectrum_q": 0.005,
            "is_decoy": 0,
        }
    ]
    pq = tmp_path / "mod_sage.parquet"
    pd.DataFrame(rows).to_parquet(pq)

    hits = sage_output.load(pq)
    assert hits.df.iloc[0]["peptide"] == seq

    result = metrics.evaluate(hits, truth, q_cutoff=0.01)
    assert result.tp == 1
    assert result.fp == 0
