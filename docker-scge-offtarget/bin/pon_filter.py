#!/usr/bin/env python3
"""
pon_filter.py — Panel-of-Normals post-filter for the genome-wide worklist.

The matched-normal check (score.py control_check) subtracts germline using ONE
donor's normal, which here is ~1/5 the tumor depth. At the VAF we care about a
single shallow normal cannot reliably see a germline/mosaic indel, so the
homology-blind LIKELY-EDIT tail is contaminated by germline that the matched
normal missed. Pooling ALL 25 unedited normals restores the depth: an indel that
is truly germline/mosaic will show up in SOME donor's normal (and in the pooled
pileup) even if the one matched normal was too shallow to call it.

Rule (textbook PoN): a candidate LIKELY-EDIT indel that appears at the observed
position in ANY donor's normal, or in the pooled cross-donor pileup, is
germline/mosaic/recurrent-artifact -> demote. A real Cas9 off-target is somatic
and donor-private: it is absent from every unedited normal.

Only re-checks sites that already survived (LIKELY EDIT / POSSIBLE); pooling is a
cheap per-site scan (~sites x 25 normals), not a genome-wide re-run.
"""
import os, sys, argparse
from collections import Counter
import numpy as np
import pandas as pd
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from features import read_records, WIN
import score as S

REF = S.REF
# a normal read "supports the indel" if it carries an in-window indel within this
# many bp of the candidate's observed modal position (matches the tumor caller's
# +/-2 bp concordance window, widened by 1 for shallow-normal position jitter)
POS_TOL = 3
PON_PRESENT_IF = 0.03   # per-normal indel_frac >= this (with >= MIN_SPAN) = present
PON_MIN_SPAN = 4        # per-normal spanning floor to count that normal as evaluable
PON_POOL_IF = 0.02      # pooled cross-donor indel_frac >= this = germline/artifact


def normals_from_map(cram_map):
    """sample -> matched-normal CRAM path (drops the tumor, keeps <base>.cram)."""
    crams = S.cram_index(cram_map)
    out = {}
    for s, tp in crams.items():
        np_ = S.normal_cram_path(tp)
        if np_:
            out[s] = np_
    return out


def normal_support(bam, chrom, pos):
    """(n_indel_near_pos, spanning) in one normal at the candidate position.

    n_indel = reads whose largest in-window indel sits within POS_TOL bp of pos;
    spanning = reads covering pos. Returns (0, 0) if the locus can't be fetched
    (e.g. contig/ref decode failure) so a bad normal never masks a real edit."""
    try:
        rr = read_records(bam, str(chrom), int(pos) - 1, int(pos) + 1)
    except OSError:
        return 0, 0
    if rr is None:
        return 0, 0
    records = rr[0]
    spanning = sum(1 for r in records if r["spans"])
    n_indel = sum(1 for r in records
                  if r["indel"] is not None and abs(r["indel"][0] - int(pos)) <= POS_TOL)
    return n_indel, spanning


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--worklist", default="wgs_offtarget_worklist_genomewide.csv")
    ap.add_argument("--cram-list",
                    default=f"{os.path.dirname(os.path.abspath(__file__))}/wgs_cram_map.tsv")
    ap.add_argument("--ref", default=REF)
    ap.add_argument("--verdicts", default="LIKELY EDIT,POSSIBLE — review",
                    help="comma-list of verdicts to re-check against the PoN")
    ap.add_argument("--out", default="wgs_offtarget_worklist_pon.csv")
    ap.add_argument("--top", type=int, default=60)
    args = ap.parse_args()

    df = pd.read_csv(args.worklist)
    normals = normals_from_map(args.cram_list)
    print(f"Panel of Normals: {len(normals)} unedited normals", flush=True)

    check_verds = set(v.strip() for v in args.verdicts.split(","))
    todo = df[df["verdict"].isin(check_verds)].copy()
    print(f"re-checking {len(todo)} candidates ({', '.join(sorted(check_verds))}) "
          f"against the PoN\n", flush=True)

    # open all normals once
    bams = {s: pysam.AlignmentFile(p, "rc", reference_filename=args.ref)
            for s, p in normals.items()}

    pon = {}  # rank -> dict of pon metrics
    for i, (_, r) in enumerate(todo.iterrows(), 1):
        chrom = str(r["chrom"])
        pos = int(r["modal_pos"]) if pd.notna(r.get("modal_pos")) else int(r["start"])
        own = r["sample"]
        n_hit = 0            # normals (excluding own, already checked) with the indel
        pooled_i = pooled_s = 0
        worst_if = 0.0
        for s, bam in bams.items():
            ni, sp = normal_support(bam, chrom, pos)
            if s == own:      # own matched normal already vetted in stage-1 ctrl_if
                continue
            if sp >= PON_MIN_SPAN:
                pooled_i += ni
                pooled_s += sp
                fi = ni / sp
                worst_if = max(worst_if, fi)
                if fi >= PON_PRESENT_IF:
                    n_hit += 1
        pool_if = (pooled_i / pooled_s) if pooled_s else 0.0
        pon[r["rank"]] = dict(pon_n_hit=n_hit, pon_pool_if=round(pool_if, 4),
                              pon_worst_if=round(worst_if, 3), pon_span=pooled_s)
        if i % 40 == 0:
            print(f"  ...{i}/{len(todo)}", flush=True)
    for b in bams.values():
        b.close()

    # merge pon metrics; sites not re-checked get NaN
    for col in ("pon_n_hit", "pon_pool_if", "pon_worst_if", "pon_span"):
        df[col] = df["rank"].map(lambda k: pon.get(k, {}).get(col, np.nan))

    # PoN verdict: demote survivors that show germline in the cross-donor panel
    def revised(r):
        v = r["verdict"]
        if v not in check_verds:
            return v
        if pd.isna(r["pon_n_hit"]):
            return v
        if r["pon_n_hit"] >= 1 or r["pon_pool_if"] >= PON_POOL_IF:
            return "GERMLINE/ARTIFACT (in PoN)"
        return v
    df["verdict_pon"] = df.apply(revised, axis=1)

    df.to_csv(args.out, index=False)

    # report
    surv = df[df["verdict_pon"].isin(check_verds)]
    demoted = df[(df["verdict"].isin(check_verds)) &
                 (df["verdict_pon"] == "GERMLINE/ARTIFACT (in PoN)")]
    le0 = df[df["verdict"] == "LIKELY EDIT"]
    le1 = df[df["verdict_pon"] == "LIKELY EDIT"]
    print(f"\nwrote {args.out}")
    print(f"LIKELY EDIT: {len(le0)} -> {len(le1)}  "
          f"(PoN demoted {len(le0) - len(le1)} as germline/artifact)")
    print(f"total survivors ({'/'.join(sorted(check_verds))}): "
          f"{len(todo)} -> {len(surv)}")

    # positive control must survive
    pc = df[df.get("truth") == True] if "truth" in df.columns else df.iloc[0:0]
    for _, k in pc.iterrows():
        print(f"POSITIVE CONTROL {k['chrom']}:{k['start']} ({k['sample']}): "
              f"verdict_pon='{k['verdict_pon']}'  pon_n_hit={k['pon_n_hit']} "
              f"pon_pool_if={k['pon_pool_if']}")

    cols = ["rank", "sample", "chrom", "start", "alt", "dragen_af", "min_mm",
            "is_target", "indel_frac", "ctrl_if", "pon_n_hit", "pon_pool_if",
            "pon_worst_if", "score", "verdict_pon"]
    cols = [c for c in cols if c in df.columns]
    surv_ot = surv[surv["is_target"] == 0] if "is_target" in surv.columns else surv
    print(f"\n== PoN survivors, off-target, top {args.top} ==")
    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(surv_ot.sort_values("score", ascending=False)[cols]
              .head(args.top).to_string(index=False))
    if "min_mm" in surv_ot.columns and len(surv_ot):
        print("\nsurvivor off-target min_mm dist:",
              surv_ot["min_mm"].value_counts().to_dict())


if __name__ == "__main__":
    main()
