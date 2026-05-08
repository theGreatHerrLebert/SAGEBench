# Skeleton PR — Sage Bruker CI smoke

Draft of a PR that wires the SAGEBench smoke fixture into `lazear/sage`'s
CI. Not yet adapted to upstream's actual workflow layout — open the
real repo, look at how their existing CI invokes sage, and shape this
to match.

## Proposed changes

1. Add a CI workflow (`.github/workflows/bruker-smoke.yml`) that:
   - downloads the SAGEBench smoke fixture (release asset / Zenodo /
     LFS — whichever the maintainer prefers);
   - builds sage in release mode;
   - runs sage against both `.d` directories with `--parquet`;
   - asserts the parquet has > 0 rows (or, more strictly, runs the
     SAGEBench harness and asserts true FDR ≤ 5%).

2. Add a small fixture-fetch script under `tests/data/` that pulls the
   smoke `.d` files and verifies their checksums. Sage's existing
   tests don't need it; this is just for the new CI job.

3. (Optional) Document the fixture in `README.md` so external users
   know how to reproduce the timsTOF integration check.

## Why this catches lazear/sage#228

The canonicalize panic surfaces specifically when sage iterates more
than one input (`processing files 0 .. 1`). The smoke fixture is
generated with two replicate `.d` directories deliberately so the
multi-input path is hit on every CI run.

## Files to adapt

- `workflows/bruker-smoke.yml` — pasted as a starting point. Sage's
  existing CI may already define a Linux job with a `sage` build
  step; if so, fold this in there rather than duplicating.
- The fixture URL is a placeholder. Replace with the chosen
  distribution mechanism once decided.
