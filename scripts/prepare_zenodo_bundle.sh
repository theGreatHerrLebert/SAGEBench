#!/usr/bin/env bash
# Bundle the SAGEBench fixtures + eval datasets for upload to Zenodo.
#
# Produces:
#
#   <out_dir>/sagebench-ci-smoke.tar.gz           (~470 MB; both 5-min .d files,
#                                                  seed CSV, configs, regen script,
#                                                  sage configs — drop-in CI fixture)
#   <out_dir>/sagebench-hela-150k-g30m.tar.gz     (~6 GB; HeLa 150K eval rep001)
#   <out_dir>/sagebench-hla-10k-g40.tar.gz        (~2.5 GB; HLA 10K, all 3 reps)
#   <out_dir>/sagebench-hla-100k-g3600.tar.gz     (~9 GB; HLA 100K, all 3 reps)
#   <out_dir>/sagebench-results.tar.gz            (~few MB; runs/REPORT.html,
#                                                  RESULTS.md, eval CSVs, draft
#                                                  issue comment — read this first)
#
# Usage:
#   bash scripts/prepare_zenodo_bundle.sh <out_dir>
#
# Skip individual bundles via SKIP, e.g.:
#   SKIP=hla-100k bash scripts/prepare_zenodo_bundle.sh /scratch/zenodo
set -euo pipefail

REPO_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
OUT_DIR="${1:-zenodo-bundle}"
SKIP="${SKIP:-}"
mkdir -p "$OUT_DIR"

# Heavy data bundles (eval sets) — symlinked .d trees.
declare -A BUNDLES=(
  [sagebench-hela-150k-g30m]="data/refs/hela-150k-g30m configs/hela-150k-g30m.toml"
  [sagebench-hla-10k-g40]="data/refs/hla-10k configs/hla-10k-g40.toml"
  [sagebench-hla-100k-g3600]="data/refs/hla-100k configs/hla-100k-g3600.toml"
)

cd "$REPO_ROOT"

# ---------------------------------------------------------------------------
# 1. CI smoke fixture — small, self-contained, the most useful Zenodo entry
#    for sage's CI. Both replicates + everything needed to regenerate.
# ---------------------------------------------------------------------------
if [[ "$SKIP" != *"ci-smoke"* ]] && [[ -d data/ci-smoke ]]; then
  smoke_stage="$OUT_DIR/_staging/sagebench-ci-smoke"
  rm -rf "$smoke_stage"
  mkdir -p "$smoke_stage"

  echo "==> staging sagebench-ci-smoke"
  # The two simulated runs, each with synthetic_data.db alongside.
  for rep in SAGEBENCH-CI-HELA-SMOKE-001 SAGEBENCH-CI-HELA-SMOKE-002; do
    if [[ -d "data/ci-smoke/$rep" ]]; then
      cp -r "data/ci-smoke/$rep" "$smoke_stage/$rep"
    fi
  done
  # Seed + configs + regen script + sage configs the user might want.
  cp seeds/hela_ci_smoke.csv "$smoke_stage/"
  cp configs/hela-ci-smoke.toml "$smoke_stage/"
  cp configs/modifications.toml "$smoke_stage/"
  cp configs/sage-smoke.json "$smoke_stage/"
  cp scripts/extract_smoke_seed.py "$smoke_stage/"
  cp scripts/generate_hela_ci_smoke.sh "$smoke_stage/"

  cat > "$smoke_stage/README.md" <<'EOF'
# SAGEBench — HeLa CI smoke fixture

Two simulated TIMS-TOF DDA `.d` files (5-min gradient, 1 500 HeLa
peptides each) plus everything needed to regenerate them and to score
any DDA search engine against the ground truth recorded in
`*/synthetic_data.db`.

## Why it exists

This fixture reproduces the multi-input timsTOF code path that
panicked in [lazear/sage#228](https://github.com/lazear/sage/issues/228)
on first invocation — sage's `processing files 0 .. 1` log line
canonicalize-then-NotFound. With the read_tdf fix applied
(see `github/ISSUE_228_COMMENT.md` in the SAGEBench repo) sage
0.15.0-beta.2 runs both files clean in ~5 seconds.

## Layout

```
sagebench-ci-smoke/
├── README.md
├── SAGEBENCH-CI-HELA-SMOKE-001/
│   ├── synthetic_data.db                        (ground truth)
│   └── SAGEBENCH-CI-HELA-SMOKE-001.d/
│       ├── analysis.tdf
│       └── analysis.tdf_bin
├── SAGEBENCH-CI-HELA-SMOKE-002/
│   └── ...
├── hela_ci_smoke.csv                            (seed peptides)
├── hela-ci-smoke.toml                           (TimSim config)
├── modifications.toml
├── sage-smoke.json                              (sage search config)
├── extract_smoke_seed.py                        (rebuild the seed CSV)
└── generate_hela_ci_smoke.sh                    (rebuild the .d files)
```

## Run sage on it

```bash
sage --parquet --output_directory ./out sage-smoke.json \
    SAGEBENCH-CI-HELA-SMOKE-001/SAGEBENCH-CI-HELA-SMOKE-001.d \
    SAGEBENCH-CI-HELA-SMOKE-002/SAGEBENCH-CI-HELA-SMOKE-002.d
```

Adjust `sage-smoke.json`'s `database.fasta` to a local UniProt human
reviewed canonical FASTA before running.

## Score against ground truth

Use the SAGEBench `sagebench` Python harness (peptide-level FDR /
TPR with multi-replicate ground-truth union):

```bash
python -m sagebench eval \
    --sage-out out/results.sage.parquet \
    --truth   SAGEBENCH-CI-HELA-SMOKE-001/synthetic_data.db \
              SAGEBENCH-CI-HELA-SMOKE-002/synthetic_data.db
```
EOF

  archive="$OUT_DIR/sagebench-ci-smoke.tar.gz"
  echo "==> compressing → $archive"
  tar -C "$OUT_DIR/_staging" -czf "$archive" sagebench-ci-smoke
  rm -rf "$smoke_stage"
fi

# ---------------------------------------------------------------------------
# 2. Heavy eval bundles — symlinked .d trees + provenance.
# ---------------------------------------------------------------------------
for name in "${!BUNDLES[@]}"; do
  if [[ "$SKIP" == *"${name#sagebench-}"* ]]; then
    echo "skip $name (SKIP=$SKIP)"
    continue
  fi
  read -r src_dir cfg <<<"${BUNDLES[$name]}"

  if [[ ! -e "$src_dir" ]]; then
    echo "skip $name — missing $src_dir (set up data/refs/ symlinks first)" >&2
    continue
  fi

  staging="$OUT_DIR/_staging/$name"
  rm -rf "$staging"
  mkdir -p "$staging"

  # Resolve symlink so tar follows it (-h would too, but we want the
  # original directory name preserved inside the archive).
  resolved=$(readlink -f "$src_dir")
  rep_basename=$(basename "$resolved")

  echo "==> staging $name from $resolved"
  cp -L -r "$resolved" "$staging/$rep_basename"
  cp "$cfg" "$staging/config.toml"
  cp "configs/modifications.toml" "$staging/modifications.toml"

  cat > "$staging/README.md" <<EOF
# $name

Simulated TIMS-TOF DDA dataset generated by TimSim, packaged for the
SAGEBench harness (https://github.com/<user>/SAGEBench — adjust before
publishing).

Contents:

- \`config.toml\` — TimSim configuration that produced the run.
- \`modifications.toml\` — variable + fixed modifications.
- \`$rep_basename/\` — the simulated \`.d\` directory.
- \`$rep_basename/synthetic_data.db\` — ground truth: every
  (sequence_unimod, charge) precursor that was injected.

Use SAGEBench's harness to compute true FDR / TPR for any DDA search
engine that can read this \`.d\`:

\`\`\`
python -m sagebench eval \\
    --sage-out <your-sage-output>.parquet \\
    --truth   $rep_basename/synthetic_data.db
\`\`\`
EOF

  archive="$OUT_DIR/$name.tar.gz"
  echo "==> compressing → $archive"
  tar -C "$staging" -czf "$archive" .
  rm -rf "$staging"
done

# ---------------------------------------------------------------------------
# 3. Results-only bundle — small, lets people see findings without
#    pulling the multi-GB .d files.
# ---------------------------------------------------------------------------
if [[ "$SKIP" != *"results"* ]]; then
  results_stage="$OUT_DIR/_staging/sagebench-results"
  rm -rf "$results_stage"
  mkdir -p "$results_stage/runs"

  echo "==> staging sagebench-results"
  cp runs/REPORT.html "$results_stage/runs/" 2>/dev/null || true
  cp runs/RESULTS.md  "$results_stage/runs/" 2>/dev/null || true
  cp sagebench-report.toml "$results_stage/" 2>/dev/null || true
  for d in runs/*/eval.csv; do
    [[ -f "$d" ]] || continue
    name=$(basename "$(dirname "$d")")
    mkdir -p "$results_stage/runs/$name"
    cp "$d" "$results_stage/runs/$name/eval.csv"
  done
  if [[ -d github ]]; then
    cp -r github "$results_stage/"
  fi

  cat > "$results_stage/README.md" <<'EOF'
# SAGEBench — results bundle

Small text-only deposit summarising what we found running sage
0.15.0-beta.2 against the simulated fixtures. Pair with the
`sagebench-*-eval.tar.gz` deposits for the underlying `.d` files
and ground-truth DBs.

- `runs/REPORT.html` — self-contained interactive report
  (calibration plots, Venns, per-level sweeps, TDC variant analysis).
  Open in any browser; all plots are inline SVG.
- `runs/RESULTS.md` — short prose summary of the same findings.
- `runs/<dataset>/eval.csv` — per-cutoff sweep CSVs (q ≤ 0.001 / 0.005
  / 0.01 / 0.05 with reported keys, TP, FP, true FDR, TPR vs full and
  vs fragmentable subset).
- `github/ISSUE_228_COMMENT.md` — draft issue comment summarising the
  whole thing for `lazear/sage#228`.
EOF

  archive="$OUT_DIR/sagebench-results.tar.gz"
  echo "==> compressing → $archive"
  tar -C "$OUT_DIR/_staging" -czf "$archive" sagebench-results
  rm -rf "$results_stage"
fi

rm -rf "$OUT_DIR/_staging"
echo
echo "Bundles in $OUT_DIR/"
ls -lh "$OUT_DIR"
