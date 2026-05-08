#!/usr/bin/env bash
# Generate the SAGEBench HeLa CI smoke fixture: 2 replicates of a
# 1 500-peptide DDA run with a 5-min gradient, derived from
# configs/hela-ci-smoke.toml (REP_PLACEHOLDER substituted per rep).
#
# Outputs land at:
#   data/ci-smoke/SAGEBENCH-CI-HELA-SMOKE-001.d/
#   data/ci-smoke/SAGEBENCH-CI-HELA-SMOKE-002.d/
#
# Usage:
#   bash scripts/generate_hela_ci_smoke.sh
#   N=4 bash scripts/generate_hela_ci_smoke.sh         # more replicates
#   DRY_RUN=1 bash scripts/generate_hela_ci_smoke.sh   # plan only
#   TIMSIM=/path/to/timsim ...                          # explicit binary
#
# Requires a TimSim build with `from_findings` support (rustims
# feature/simulate-from-findings, commit 03398133+). On this host:
#   TIMSIM=/scratch/TMAlign/MHCbench/.venv/bin/timsim bash $0
set -euo pipefail

REPO_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
N="${N:-2}"
DRY_RUN="${DRY_RUN:-0}"
TIMSIM="${TIMSIM:-timsim}"

cd "$REPO_ROOT"

BASE_CFG="configs/hela-ci-smoke.toml"
SEED_CSV="seeds/hela_ci_smoke.csv"
OUT_DIR="data/ci-smoke"
RENDERED_DIR="$OUT_DIR/_rendered_configs"

if [[ ! -f "$BASE_CFG" ]]; then
  echo "Missing $BASE_CFG" >&2
  exit 1
fi

if [[ ! -f "$SEED_CSV" ]]; then
  echo "Missing $SEED_CSV — run scripts/extract_smoke_seed.py first." >&2
  exit 1
fi

if [[ "$DRY_RUN" != "1" ]] && ! command -v "$TIMSIM" >/dev/null 2>&1; then
  cat >&2 <<EOF
ERROR: '$TIMSIM' is not on PATH.

The smoke fixture needs a TimSim build with the from_findings code
path (rustims feature/simulate-from-findings, commit 03398133+).
Activate a venv that ships it, or set TIMSIM=/abs/path/to/timsim.
EOF
  exit 2
fi

mkdir -p "$OUT_DIR" "$RENDERED_DIR"

for i in $(seq 1 "$N"); do
  rep=$(printf '%03d' "$i")
  rendered="$RENDERED_DIR/hela-ci-smoke-rep$rep.toml"
  sed "s/REP_PLACEHOLDER/$rep/g" "$BASE_CFG" > "$rendered"
  echo
  echo "==> rep$rep — $TIMSIM $rendered"
  if [[ "$DRY_RUN" == "1" ]]; then
    echo "    (dry-run; skipping)"
  else
    "$TIMSIM" "$rendered"
  fi
done

echo
echo "CI smoke fixtures under $OUT_DIR/"
