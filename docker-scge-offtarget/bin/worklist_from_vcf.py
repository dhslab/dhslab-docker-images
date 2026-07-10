#!/usr/bin/env python3
"""
worklist_from_vcf.py — GENOME-WIDE, HOMOLOGY-FREE off-target worklist.

Candidate source = DRAGEN somatic tumor/normal small-variant VCFs
(`<sample>.hard-filtered.vcf.gz`, tumor=edited, normal=<donor>-unedited). These are
genome-wide indels already normal-subtracted + systematic-noise-filtered by DRAGEN —
so candidates come from edited-vs-control evidence, NOT a guide-homology list. Guide
homology (min_mm / is_target) is attached only as an ANNOTATION by cross-referencing
the Cas-OFFinder table, so homology-blind real off-targets can't be excluded.

Per PASS somatic indel: walk the tumor CRAM -> count-free pileup shape -> shape
ranker; recompute in matched normal at the observed indel pos (cross-check on top of
DRAGEN's somatic call); cross-sample recurrence; verdict. Optionally auto-render an
edited-vs-normal pileup PNG per LIKELY EDIT.

Reuses score.py (control_check, normal_cram_path, add_recurrence, verdict, HI/MID)
and features.py. Known PLCB2 off-target (chr12:32679408) = built-in positive control.
"""
import os, sys, argparse
import numpy as np
import pandas as pd
import pysam
import joblib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from features import read_records, features_from_records, MODEL_FEATURES, parse_target_info, check_sklearn_version
import score as S

REF = S.REF
MODEL = S.MODEL
STD_CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}


def vcf_for_tumor(tumor_cram):
    """<dir>/<BASE>_tumor.cram -> <dir>/<BASE>.hard-filtered.vcf.gz"""
    if tumor_cram.endswith("_tumor.cram"):
        base = tumor_cram[:-len("_tumor.cram")]
        v = base + ".hard-filtered.vcf.gz"
        return v if os.path.exists(v) else None
    return None


def tumor_sample_id(vcf):
    """DRAGEN somatic: the non-'*-unedited' sample column is the tumor/edited one."""
    v = pysam.VariantFile(vcf)
    ids = list(v.header.samples)
    v.close()
    tum = [s for s in ids if not s.endswith("-unedited")]
    return tum[-1] if tum else (ids[-1] if ids else None)


def is_indel(rec):
    return any(len(a) != len(rec.ref) for a in (rec.alts or []))


def build_homology_index(table_path, sheet=None, pad=25):
    """chrom -> sorted list of (pos, min_mm, is_target) from the Cas-OFFinder table,
    for annotating VCF indels with the nearest predicted-homology site (within pad)."""
    df = S.load_table(table_path, sheet=sheet)
    idx = {}
    for _, r in df.iterrows():
        idx.setdefault(str(r["chrom"]), []).append(
            (int(r["start"]), int(r["_min_mm"]), int(bool(r["_tgt"]))))
    for c in idx:
        idx[c].sort()
    return idx, pad


def annotate_homology(idx_pad, chrom, pos):
    """Nearest table site within pad -> (min_mm, is_target, dist) or (99, 0, None)."""
    if idx_pad is None:
        return 99, 0, None
    idx, pad = idx_pad
    lst = idx.get(str(chrom))
    if not lst:
        return 99, 0, None
    import bisect
    starts = [x[0] for x in lst]
    i = bisect.bisect_left(starts, pos)
    best = None
    for j in (i - 1, i, i + 1):
        if 0 <= j < len(lst):
            d = abs(lst[j][0] - pos)
            if d <= pad and (best is None or d < best[2]):
                best = (lst[j][1], lst[j][2], d)
    return best if best else (99, 0, None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cram-list", default=f"{os.path.dirname(os.path.abspath(__file__))}/wgs_cram_map.tsv")
    ap.add_argument("--ref", default=REF)
    ap.add_argument("--model", default=MODEL)
    ap.add_argument("--homology-table",
                    default="/storage2/fs1/dspencer/Active/clinseq/projects/scge/"
                            "Manual Indel Review/cart_wgs/cart_wgs_merged.xlsx",
                    help="Cas-OFFinder table for min_mm/is_target ANNOTATION only")
    ap.add_argument("--sample", help="limit to one sample_name (testing)")
    ap.add_argument("--min-af", type=float, default=0.05, help="DRAGEN tumor AF floor")
    ap.add_argument("--min-span", type=int, default=8)
    ap.add_argument("--include-weak", action="store_true",
                    help="also include weak_evidence-filtered indels (depth-floor tail)")
    ap.add_argument("--no-normal-check", action="store_true")
    ap.add_argument("--snapshot-dir")
    ap.add_argument("--top", type=int, default=50)
    ap.add_argument("--out", default="wgs_offtarget_worklist_genomewide.csv")
    args = ap.parse_args()

    bundle = joblib.load(args.model)
    model, feat_names = bundle["model"], bundle["features"]
    check_sklearn_version(model, name=os.path.basename(args.model))
    crams = S.cram_index(args.cram_list)
    if args.sample:
        crams = {args.sample: crams[args.sample]}
    homidx = build_homology_index(args.homology_table) if args.homology_table else None

    rows, n_pass_indel = [], 0
    for s, tpath in crams.items():
        vcf = vcf_for_tumor(tpath)
        if not vcf:
            print(f"  {s}: no VCF", flush=True); continue
        tid = tumor_sample_id(vcf)
        npath = S.normal_cram_path(tpath)
        tbam = pysam.AlignmentFile(tpath, "rc", reference_filename=args.ref)
        ncache = {}
        v = pysam.VariantFile(vcf)
        sc = 0
        for rec in v.fetch():
            if not is_indel(rec):
                continue
            filt = list(rec.filter.keys())
            ispass = (not filt) or ("PASS" in filt)
            if not ispass and not (args.include_weak and "weak_evidence" in filt):
                continue
            af = rec.samples[tid].get("AF")
            af = (af[0] if isinstance(af, tuple) else af) or 0.0
            if af < args.min_af:
                continue
            chrom, pos = str(rec.chrom), int(rec.pos)
            # off-targets are genomic; skip vector/transgene + non-standard contigs
            # (the DRAGEN ref carries the CAR vector, which our hg38 REF can't decode)
            if chrom not in STD_CHROMS:
                continue
            n_pass_indel += 1
            try:
                rr = read_records(tbam, chrom, pos - 1, pos + 1)
                feats = None if rr is None else features_from_records(
                    rr[0], min_span=args.min_span, return_lowcov=True)
            except OSError:                      # CRAM/ref decode failure at this locus
                continue
            mm, tgt, dist = annotate_homology(homidx, chrom, pos)
            rec_d = {"sample": s, "chrom": chrom, "start": pos, "ref": rec.ref,
                     "alt": ",".join(rec.alts or []), "dragen_af": round(float(af), 3),
                     "min_mm": mm, "is_target": tgt, "hom_dist": dist, "truth": False}
            if feats is None or feats.get("lowcov"):
                vd, _ = S.verdict(0.0, feats)
                rec_d.update(score=np.nan, verdict=vd, spanning=(feats or {}).get("spanning", 0),
                             conc_ratio=np.nan, indel_frac=np.nan, modal_len=np.nan,
                             modal_pos=np.nan, ctrl_if=np.nan)
            else:
                psc = float(model.predict_proba(pd.DataFrame([{k: feats[k] for k in feat_names}]))[0, 1])
                opos = feats.get("modal_pos") or pos
                ctrl_if = None
                if not args.no_normal_check:
                    ctrl_if, _ = S.control_check(npath, chrom, opos, args.ref, bam_cache=ncache)
                vd, psc = S.verdict(psc, feats, ctrl_if)
                rec_d.update(score=psc, verdict=vd, spanning=feats["spanning"],
                             conc_ratio=round(feats["conc_ratio"], 3),
                             indel_frac=round(feats["indel_frac"], 3),
                             modal_len=feats.get("modal_len"), modal_pos=opos,
                             ctrl_if=(round(ctrl_if, 3) if ctrl_if is not None else np.nan))
            # built-in positive control
            if chrom == "chr12" and abs(pos - 32679408) <= 3 and "PLCB2" in s:
                rec_d["truth"] = True
            rows.append(rec_d); sc += 1
        v.close(); tbam.close()
        for b in ncache.values():
            b.close()
        print(f"  {s}: {sc} PASS somatic indels scored (AF>={args.min_af})", flush=True)

    res = S.add_recurrence(pd.DataFrame(rows))
    torder = {"A-review-first": 0, "B-review": 1, "C-recurrent-artifact": 2, "D-artifact": 3}
    res["_o"] = res["priority"].map(torder).fillna(3)
    res = res.sort_values(["_o", "score"], ascending=[True, False],
                          na_position="last").drop(columns="_o").reset_index(drop=True)
    res.insert(0, "rank", res.index + 1)
    res.to_csv(args.out, index=False)
    print(f"\ngenome-wide PASS somatic indels featurized: {n_pass_indel}")
    print(f"wrote {args.out} ({len(res)} rows)")

    cols = ["rank", "priority", "sample", "chrom", "start", "alt", "dragen_af",
            "indel_frac", "ctrl_if", "min_mm", "is_target", "conc_ratio", "score",
            "verdict", "n_guides_at_site", "truth"]
    cols = [c for c in cols if c in res.columns]
    print(f"\nTop {min(args.top, len(res))} (tiered):")
    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(res[cols].head(args.top).to_string(index=False))
    print("\ntiers: " + res["priority"].value_counts().to_string())
    # homology-blind survivors = the GUIDE-seq-analog finds
    le = res[res.verdict == "LIKELY EDIT"]
    blind = le[le.min_mm == 99]
    print(f"\nLIKELY EDIT: {len(le)}  | of those homology-blind (no table site within 25bp): {len(blind)}")
    if res["truth"].any():
        for _, k in res[res.truth].iterrows():
            print(f"POSITIVE CONTROL PLCB2: rank {int(k['rank'])}/{len(res)} score={k['score']} verdict='{k['verdict']}'")

    if args.snapshot_dir:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from pileup_snapshot import snapshot
        os.makedirs(args.snapshot_dir, exist_ok=True)
        surv = res[res.verdict == "LIKELY EDIT"]
        for _, k in surv.iterrows():
            tp = crams.get(k["sample"])
            if not tp:
                continue
            pos = int(k["modal_pos"]) if pd.notna(k.get("modal_pos")) else int(k["start"])
            npath = S.normal_cram_path(tp)
            fig, ax = plt.subplots(1, 2 if npath else 1, figsize=(15 if npath else 8, 6.5), squeeze=False)
            snapshot(ax[0][0], tp, str(k["chrom"]), pos, args.ref, window=45,
                     title=f"{k['sample']}_tumor  AF={k['dragen_af']} score={k['score']:.2f}")
            if npath:
                snapshot(ax[0][1], npath, str(k["chrom"]), pos, args.ref, window=45,
                         title=f"NORMAL  ctrl_if={k.get('ctrl_if')}")
            fig.suptitle(f"{k['sample']}  {k['chrom']}:{pos:,}  min_mm={k['min_mm']}  {k['verdict']}",
                         fontweight="bold")
            fig.tight_layout()
            fig.savefig(f"{args.snapshot_dir}/rank{int(k['rank']):03d}_{k['sample']}_{k['chrom']}_{pos}.png", dpi=120)
            plt.close(fig)
        print(f"\nrendered {len(surv)} snapshots -> {args.snapshot_dir}/")


if __name__ == "__main__":
    main()
