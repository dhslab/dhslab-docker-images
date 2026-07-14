#!/usr/bin/env python3
"""
recall_vs_vaf.py — the honest limit, measured.

From the ECS⋈WGS training table, compute how often the WGS shape score recovers an
ECS-confirmed edit as a function of the ECS (error-corrected) VAF. This is THE
deliverable that keeps the WGS-only promise credible: it names the VAF above which
WGS-only detection is trustworthy, rather than implying WGS sees everything ECS sees.

WGS "detected" = shape score >= --hi (the LIKELY-EDIT threshold) AND the locus was
evaluable (not INSUFFICIENT COVERAGE). ECS edits below the WGS depth floor simply have
no supporting WGS reads and will show up as the low-VAF recall shortfall — by design.
"""
import argparse
import numpy as np
import pandas as pd

BINS = [0.0, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.01]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--training", required=True, help="training.tsv from join_training_table.py")
    ap.add_argument("--hi", type=float, default=0.60, help="WGS score >= HI counts as detected")
    ap.add_argument("--target-recall", type=float, default=0.80,
                    help="report the VAF floor where binned recall first reaches this")
    ap.add_argument("--out-metrics", default="recall_vs_vaf.csv")
    ap.add_argument("--out-curve", default="recall_vs_vaf.png")
    args = ap.parse_args()

    df = pd.read_csv(args.training, sep="\t")
    pos = df[df["label"] == 1].copy()
    if len(pos) == 0:
        print("no ECS-positive hotspots in training table; nothing to measure")
        pd.DataFrame(columns=["vaf_bin", "n", "recall"]).to_csv(args.out_metrics, index=False)
        return

    pos["score"] = pd.to_numeric(pos["score"], errors="coerce")
    pos["detected"] = (pos["score"] >= args.hi).fillna(False).astype(int)
    pos["vaf_bin"] = pd.cut(pd.to_numeric(pos["ecs_if"], errors="coerce"), bins=BINS,
                            right=False)

    g = (pos.groupby("vaf_bin", observed=True)
            .agg(n=("detected", "size"), n_detected=("detected", "sum"))
            .reset_index())
    g["recall"] = g["n_detected"] / g["n"]
    g.to_csv(args.out_metrics, index=False)

    overall = pos["detected"].mean()
    floor = None
    for _, r in g.iterrows():
        if r["n"] >= 1 and r["recall"] >= args.target_recall:
            floor = r["vaf_bin"].left
            break
    print(f"overall WGS recall of ECS-confirmed edits (score>={args.hi}): {overall:.2f} "
          f"over {len(pos)} sites")
    print(f"VAF floor for >= {args.target_recall:.0%} binned recall: "
          f"{'>%.3f' % floor if floor is not None else 'not reached in these bins'}")
    print(g.to_string(index=False))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        centers = [iv.left for iv in g["vaf_bin"]]
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot(centers, g["recall"], "o-", color="#2b6cb0")
        ax.axhline(args.target_recall, ls="--", color="#a0aec0",
                   label=f"target recall {args.target_recall:.0%}")
        ax.set_xscale("symlog", linthresh=0.005)
        ax.set_xlabel("ECS error-corrected VAF (lower bin edge)")
        ax.set_ylabel(f"WGS recall (score ≥ {args.hi})")
        ax.set_ylim(-0.02, 1.02)
        ax.set_title("WGS-only recovery of ECS-confirmed edits vs VAF")
        for _, r in g.iterrows():
            ax.annotate(f"n={r['n']}", (r["vaf_bin"].left, r["recall"]),
                        textcoords="offset points", xytext=(0, 6), fontsize=8, ha="center")
        ax.legend()
        fig.tight_layout()
        fig.savefig(args.out_curve, dpi=130)
        print(f"wrote {args.out_curve}")
    except Exception as e:  # plotting is a nicety; metrics CSV is the source of truth
        print(f"(curve not rendered: {e})")


if __name__ == "__main__":
    main()
