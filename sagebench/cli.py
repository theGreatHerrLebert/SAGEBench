"""`python -m sagebench` — score Sage output against TimSim ground truth."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import ground_truth as gt
from . import report as report_mod
from . import sage_output, metrics


def _cmd_eval(args: argparse.Namespace) -> int:
    truths = [gt.load(p) for p in args.truth]
    truth = gt.union(*truths)
    hits = sage_output.load(args.sage_out)

    print(f"ground truth:    {truth.n_targets:>8,} (sequence, charge) pairs")
    print(f"sage PSMs total: {hits.n_psms:>8,} ({hits.n_targets:,} target / "
          f"{hits.n_psms - hits.n_targets:,} decoy)")
    print()

    table = metrics.sweep(hits, truth, q_cutoffs=args.q)
    pretty = table.copy()
    pretty["true_fdr"] = pretty["true_fdr"].map(lambda x: f"{x:.4f}")
    pretty["tpr"] = pretty["tpr"].map(lambda x: f"{x:.4f}")
    print(pretty.to_string(index=False))

    if args.csv is not None:
        # Re-evaluate with floats to preserve precision in the CSV.
        out = metrics.sweep(hits, truth, q_cutoffs=args.q)
        out.to_csv(args.csv, index=False)
        print(f"\nwrote {args.csv}")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="sagebench")
    sub = ap.add_subparsers(dest="cmd", required=True)

    ev = sub.add_parser(
        "eval",
        help="Compute true FDR / TPR for sage hits against TimSim ground truth.",
    )
    ev.add_argument(
        "--sage-out",
        required=True,
        type=Path,
        help="Sage parquet (results.sage.parquet)",
    )
    ev.add_argument(
        "--truth",
        required=True,
        type=Path,
        nargs="+",
        help="One or more synthetic_data.db files (one per simulated .d)",
    )
    ev.add_argument(
        "--q",
        type=float,
        nargs="+",
        default=[0.001, 0.005, 0.01, 0.05],
        help="q-value cutoffs to report",
    )
    ev.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="optional: write the report table to CSV",
    )
    ev.set_defaults(func=_cmd_eval)

    rp = sub.add_parser(
        "report",
        help="Build a self-contained HTML report from a TOML run inventory.",
    )
    rp.add_argument(
        "--config",
        required=True,
        type=Path,
        help="TOML inventory (see sagebench-report.toml)",
    )
    rp.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Output HTML path",
    )
    rp.set_defaults(func=_cmd_report)

    args = ap.parse_args(argv)
    return args.func(args)


def _cmd_report(args: argparse.Namespace) -> int:
    out = report_mod.build(args.config, args.out)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
