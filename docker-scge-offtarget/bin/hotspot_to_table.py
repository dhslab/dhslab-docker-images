#!/usr/bin/env python3
"""
hotspot_to_table.py — build a score.py-ready loci table to score WGS CRAMs at the
ECS-defined hotspot panel, and emit the matching ECS truth.

The ECS `*.offtarget_analysis.tsv` files ARE the hotspot panel: each row carries a
predicted site (chrom,start,end,target_info,is_target) plus the ECS error-corrected
edit fraction (`indel_fraction`) = ground truth VAF at that site. This script:
  1. collapses those TSVs to the unique hotspot panel PER GUIDE,
  2. cross-joins each guide's hotspots with that guide's WGS sample(s),
  3. writes a score.py table whose `sample_name` = the WGS sample (so score.py looks up
     the WGS CRAM). `indel_fraction` is forced to 1.0 so score.py's stage-1 gate keeps
     every hotspot — the real WGS signal is recomputed from the CRAM, this is a
     gate-passer only.
Separately writes the ECS truth (per guide/site VAF) for the training-table and
recall-vs-VAF joins, keyed on (guide, chrom, start).

Guide<->sample mapping comes from the unified samplesheet (columns: sample,datatype,guide).
ECS table filenames are `<ecs_sample>.offtarget_analysis.tsv`.
"""
import os, sys, argparse, glob
import pandas as pd

HOTSPOT_COLS = ["chrom", "start", "end", "is_target", "target_info"]


def load_samplesheet(path):
    df = pd.read_csv(path)
    df.columns = [c.strip().lower() for c in df.columns]
    need = {"sample", "datatype", "guide"}
    if not need.issubset(df.columns):
        sys.exit(f"samplesheet must have columns {need}; got {list(df.columns)}")
    df["datatype"] = df["datatype"].str.lower().str.strip()
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ecs-tables", nargs="+", required=True,
                    help="ECS <sample>.offtarget_analysis.tsv files (the hotspot panel + truth)")
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--edit-threshold", type=float, default=0.0,
                    help="ECS indel_fraction strictly above this = a true edit (label=1)")
    ap.add_argument("--out-table", default="wgs_hotspot_input_table.csv")
    ap.add_argument("--out-truth", default="ecs_hotspot_truth.csv")
    args = ap.parse_args()

    ss = load_samplesheet(args.samplesheet)
    ecs_sample_guide = dict(zip(ss.loc[ss.datatype == "ecs", "sample"],
                                ss.loc[ss.datatype == "ecs", "guide"]))
    wgs_by_guide = (ss[ss.datatype == "wgs"].groupby("guide")["sample"].apply(list).to_dict())

    tables = []
    for f in args.ecs_tables:
        ecs_sample = os.path.basename(f).replace(".offtarget_analysis.tsv", "")
        guide = ecs_sample_guide.get(ecs_sample)
        if guide is None:
            print(f"  WARN: {ecs_sample} not in samplesheet as ecs; skipping {f}", flush=True)
            continue
        t = pd.read_csv(f, sep="\t")
        t["guide"] = guide
        t["ecs_if"] = pd.to_numeric(t["indel_fraction"], errors="coerce").fillna(0.0)
        tables.append(t)
    if not tables:
        sys.exit("no ECS tables matched samplesheet ecs samples")
    ecs = pd.concat(tables, ignore_index=True)

    # ECS truth per (guide, site): take the strongest ECS evidence across replicates
    truth = (ecs.groupby(["guide"] + HOTSPOT_COLS, as_index=False, dropna=False)["ecs_if"]
                .max())
    truth["ecs_is_edit"] = (truth["ecs_if"] > args.edit_threshold).astype(int)
    truth.to_csv(args.out_truth, index=False)

    # score.py input: unique hotspot panel per guide x that guide's WGS sample(s)
    panel = ecs[["guide"] + HOTSPOT_COLS].drop_duplicates()
    rows = []
    for _, h in panel.iterrows():
        for w in wgs_by_guide.get(h["guide"], []):
            rows.append({"sample_name": w, "chrom": h["chrom"], "start": int(h["start"]),
                         "end": int(h["end"]), "indel_fraction": 1.0,
                         "control_indel_fraction": 0.0, "is_target": h["is_target"],
                         "target_info": h["target_info"], "guide": h["guide"]})
    out = pd.DataFrame(rows)
    out.to_csv(args.out_table, index=False)
    print(f"wrote {args.out_table}: {len(out)} (wgs_sample x hotspot) rows "
          f"across {out['guide'].nunique() if len(out) else 0} guides")
    print(f"wrote {args.out_truth}: {len(truth)} (guide x hotspot) truth rows "
          f"({int(truth['ecs_is_edit'].sum())} ECS-positive)")


if __name__ == "__main__":
    main()
