#!/usr/bin/env bash
# Run Sage against one or more simulated .d directories and emit a
# parquet result table in the chosen output dir.
#
# Usage:
#   bash scripts/run_sage.sh <out_dir> <input.d> [<input.d> ...]
#
# Env:
#   SAGE          Path to the sage binary (default: `sage` on PATH).
#   SAGE_CONFIG   Path to the sage JSON config (default: configs/sage.json).
#   SAGE_FASTA    Path to the FASTA used for the search (default: pulled
#                 from SAGE_CONFIG).
#
# The sage-config-on-disk is left untouched; we only override --output_directory.
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "usage: $0 <out_dir> <input.d> [<input.d> ...]" >&2
  exit 64
fi

OUT_DIR="$1"; shift
INPUTS=("$@")

SAGE="${SAGE:-sage}"
SAGE_CONFIG="${SAGE_CONFIG:-configs/sage-smoke.json}"

if ! command -v "$SAGE" >/dev/null 2>&1; then
  echo "ERROR: '$SAGE' is not on PATH (set SAGE=/abs/path/to/sage)." >&2
  exit 2
fi

if [[ ! -f "$SAGE_CONFIG" ]]; then
  echo "ERROR: missing sage config at $SAGE_CONFIG (set SAGE_CONFIG=...)." >&2
  exit 2
fi

mkdir -p "$OUT_DIR"

echo "==> sage:        $($SAGE --version 2>&1 | head -1)"
echo "==> config:      $SAGE_CONFIG"
echo "==> output:      $OUT_DIR"
echo "==> inputs:      ${#INPUTS[@]} file(s)"
for f in "${INPUTS[@]}"; do
  echo "                 - $f"
done
echo

SAGE_FLAGS=("--parquet" "--annotate-matches")
[[ "${SAGE_NO_ANNOTATE:-0}" == "1" ]] && SAGE_FLAGS=("--parquet")

# Sage takes positional .d/.mzML inputs after --output_directory.
"$SAGE" "${SAGE_FLAGS[@]}" --output_directory "$OUT_DIR" "$SAGE_CONFIG" "${INPUTS[@]}"

echo
echo "==> done. Parquet at $OUT_DIR/results.sage.parquet (or similar)."
