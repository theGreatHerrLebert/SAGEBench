"""Mokapot feature ablation against TimSim ground truth.

Trains mokapot/linear on sage's PSM parquet with the full feature set,
ranks features by absolute LinearSVC weight (averaged across folds),
then iteratively drops the top-k highest-weighted features and re-trains.
For each ablation prints peptide-level true FDR / TPR at q ≤ 0.01 so
we can see whether a single feature is responsible for the calibration
drift, or whether it's diffuse.

Run:
    PYTHONPATH=. python scripts/ablate_features.py \\
        --sage-parquet runs/hela-150k-g30m/results.sage.parquet \\
        --truth        data/refs/hela-150k-g30m/synthetic_data.db \\
        --max-drop 8
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# mokapot 0.10.0 still uses np.float_, removed in numpy 2.0. Mirror
# sagepy-rescore's compatibility shim before importing mokapot.
if not hasattr(np, "float_"):
    np.float_ = np.float64

import mokapot  # noqa: E402

# Reuse sagepy-rescore's feature builder + mokapot wrapper.
sys.path.insert(0, "/scratch/sagepy/sagepy-rescore/src")
from sagepy_rescore import features as feats  # type: ignore  # noqa: E402
from sagepy_rescore import mokapot_runner  # type: ignore  # noqa: E402

from sagebench import ground_truth, sage_output  # noqa: E402


def _build_dataset(sage_parquet: Path, feature_columns: list[str]) -> mokapot.LinearPsmDataset:
    df = pd.read_parquet(sage_parquet)
    df = feats.normalize_sage_cli_columns(df)
    bundle = feats.build_feature_frame(df, feature_columns=feature_columns)
    feat_df = bundle["frame"]
    used = bundle["feature_columns"]

    ds = mokapot.LinearPsmDataset(
        psms=feat_df,
        target_column=bundle["target_column"],
        spectrum_columns=bundle["spectrum_columns"],
        peptide_column=bundle["peptide_column"],
        feature_columns=used,
    )
    return ds, used


def _peptide_eval(out_files: list[str], gt) -> dict:
    """Score mokapot peptide TSV against truth at q ≤ 0.01 (peptide level)."""
    pep_path = next(p for p in out_files if p.endswith("mokapot.peptides.txt"))
    peptides = pd.read_csv(pep_path, sep="\t")
    peptides["peptide_unimod"] = peptides["sequence_modified"].astype(str).map(
        sage_output._peaks_to_unimod
    )
    target = peptides[peptides["target"] == True]
    sel = target[target["mokapot q-value"] <= 0.01]
    keys = set(sel["peptide_unimod"])
    tp = len(keys & gt.peptide_keys)
    fp = len(keys) - tp
    return {
        "keys": len(keys),
        "tp": tp,
        "fp": fp,
        "true_fdr": fp / max(1, tp + fp),
        "tpr_full": tp / max(1, gt.n_peptides),
        "tpr_frag": tp / max(1, gt.n_fragmentable_at("peptide")),
    }


def _linear_weights(models) -> dict[str, float]:
    """Mean |weight| per feature across mokapot's CV-fold models.

    mokapot's LinearSVC sits in a sklearn Pipeline whose final step is
    the actual estimator. We pull `coef_` from there.
    """
    weights = {}
    for m in models:
        est = m.estimator
        if hasattr(est, "named_steps"):
            est = list(est.named_steps.values())[-1]
        if hasattr(est, "coef_"):
            coef = np.asarray(est.coef_).ravel()
            for name, w in zip(m.features, coef, strict=True):
                weights.setdefault(name, []).append(abs(float(w)))
    return {k: float(np.mean(v)) for k, v in weights.items()}


def _train(features_to_use: list[str], sage_parquet: Path,
           train_fdr: float = 0.01, output_dir: Path | None = None):
    df = pd.read_parquet(sage_parquet)
    df = feats.normalize_sage_cli_columns(df)
    bundle = feats.build_feature_frame(df, feature_columns=features_to_use)
    out = mokapot_runner.rescore(
        bundle,
        output_dir=str(output_dir) if output_dir else "/tmp/_mokapot_ablate",
        folds=3,
        train_fdr=train_fdr,
        model_kind="linear",
        verbose=False,
    )
    return out["results"], out["models"], bundle["feature_columns"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sage-parquet", required=True, type=Path)
    ap.add_argument("--truth", required=True, type=Path)
    ap.add_argument("--max-drop", type=int, default=6,
                    help="how many of the top features to drop, one at a time")
    args = ap.parse_args()

    gt = ground_truth.load(args.truth)
    print(f"truth peptides: {gt.n_peptides:,}  fragmentable: {gt.n_fragmentable_at('peptide'):,}")
    print()

    # Pass 0: full feature set, get weights for ranking.
    base_features = list(feats._DEFAULT_FEATURES)
    print(f"[round 0] training on {len(base_features)} features (full set)…")
    out_dir0 = Path("/tmp/_mokapot_ablate_round0")
    results, models, used = _train(base_features, args.sage_parquet, output_dir=out_dir0)
    weights = _linear_weights(models)
    out_files0 = list(map(str, out_dir0.glob("mokapot*.txt")))
    metric = _peptide_eval(out_files0, gt)

    ranking = sorted(weights.items(), key=lambda kv: kv[1], reverse=True)
    print(f"  q≤0.01 peptides: keys={metric['keys']:5d}  TP={metric['tp']:5d}  "
          f"FP={metric['fp']:4d}  true_FDR={metric['true_fdr']*100:5.2f}%  "
          f"TPR_frag={metric['tpr_frag']*100:5.2f}%")
    print("  feature ranking by |weight|:")
    for name, w in ranking:
        print(f"    {w:.3f}  {name}")
    print()

    # Drop top-k progressively.
    rows = [{"dropped": "(none)", "n_features": len(used),
             "keys": metric["keys"], "tp": metric["tp"], "fp": metric["fp"],
             "true_fdr": metric["true_fdr"], "tpr_frag": metric["tpr_frag"]}]
    keep = list(used)
    for k in range(1, args.max_drop + 1):
        feat_to_drop = ranking[k - 1][0]
        if feat_to_drop not in keep:
            print(f"[round {k}] {feat_to_drop} already absent, skipping")
            continue
        keep = [f for f in keep if f != feat_to_drop]
        print(f"[round {k}] dropping {feat_to_drop} (|w|={ranking[k-1][1]:.3f}) "
              f"→ {len(keep)} features")
        out_dir_k = Path(f"/tmp/_mokapot_ablate_round{k}")
        results, models, used = _train(keep, args.sage_parquet, output_dir=out_dir_k)
        m = _peptide_eval(list(map(str, out_dir_k.glob("mokapot*.txt"))), gt)
        rows.append({"dropped": feat_to_drop, "n_features": len(used),
                     "keys": m["keys"], "tp": m["tp"], "fp": m["fp"],
                     "true_fdr": m["true_fdr"], "tpr_frag": m["tpr_frag"]})
        print(f"  q≤0.01 peptides: keys={m['keys']:5d}  TP={m['tp']:5d}  "
              f"FP={m['fp']:4d}  true_FDR={m['true_fdr']*100:5.2f}%  "
              f"TPR_frag={m['tpr_frag']*100:5.2f}%")

    print()
    print("=== ablation summary ===")
    summary = pd.DataFrame(rows)
    summary["true_fdr"] = summary["true_fdr"].map(lambda v: f"{v*100:.2f} %")
    summary["tpr_frag"] = summary["tpr_frag"].map(lambda v: f"{v*100:.2f} %")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
