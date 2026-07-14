#!/usr/bin/env python3
"""
reconcile_offtarget_report.py — the joined off-target deliverable.

Annotates the genome-wide, PoN-filtered WGS worklist with:
  - is_hotspot     : the candidate falls within --pad bp of a predicted (ECS-panel) site
  - ecs_confirmed  : that hotspot has ECS error-corrected evidence (ecs_if > 0)
  - ecs_if         : the ECS VAF at the matched hotspot (NaN if none)
so a reviewer sees at a glance whether a homology-free WGS hit is a known predicted
site (ECS-backed) or a novel, unpredicted candidate. Writes the annotated worklist and
a short summary.
"""
import argparse
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--worklist", required=True, help="PoN-filtered genome-wide worklist CSV")
    ap.add_argument("--truth", help="ecs_hotspot_truth.csv (guide,chrom,start,ecs_if,...)")
    ap.add_argument("--pad", type=int, default=25, help="bp window to match a candidate to a hotspot")
    ap.add_argument("--verdict-col", default="verdict_pon",
                    help="which verdict column to summarize (falls back to 'verdict')")
    ap.add_argument("--out", default="offtarget_report.csv")
    args = ap.parse_args()

    wl = pd.read_csv(args.worklist)
    vcol = args.verdict_col if args.verdict_col in wl.columns else "verdict"

    wl["is_hotspot"] = 0
    wl["ecs_confirmed"] = 0
    wl["ecs_if"] = np.nan

    if args.truth:
        truth = pd.read_csv(args.truth)
        truth["chrom"] = truth["chrom"].astype(str)
        # index hotspots per chrom for a padded nearest-site match
        by_chrom = {c: t.sort_values("start") for c, t in truth.groupby("chrom")}
        for i, r in wl.iterrows():
            t = by_chrom.get(str(r["chrom"]))
            if t is None:
                continue
            d = (t["start"] - int(r["start"])).abs()
            j = d.idxmin() if len(d) else None
            if j is not None and d.loc[j] <= args.pad:
                wl.at[i, "is_hotspot"] = 1
                ecs_if = float(t.loc[j, "ecs_if"])
                wl.at[i, "ecs_if"] = ecs_if
                wl.at[i, "ecs_confirmed"] = int(ecs_if > 0)

    wl.to_csv(args.out, index=False)

    n = len(wl)
    likely = wl[wl[vcol].astype(str).str.contains("LIKELY EDIT", na=False)]
    print(f"wrote {args.out} ({n} candidates)")
    print(f"  LIKELY EDIT ({vcol}): {len(likely)}")
    print(f"  on a predicted hotspot: {int(wl['is_hotspot'].sum())}  "
          f"(ECS-confirmed: {int(wl['ecs_confirmed'].sum())})")
    novel = likely[likely["is_hotspot"] == 0]
    print(f"  LIKELY EDIT & NOT on any predicted hotspot (novel candidates): {len(novel)}")
    if len(novel):
        cols = [c for c in ["rank", "sample", "chrom", "start", "alt", "score", vcol]
                if c in novel.columns]
        print(novel[cols].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
