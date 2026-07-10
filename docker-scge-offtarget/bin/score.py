#!/usr/bin/env python3
"""
pileup/score.py — deploy the middle model: Stage-1 filter + learned Stage-2
shape ranker, on ANY candidate table (ECS .csv.gz or WGS .xlsx/.xlsm) + CRAMs.

Pipeline (matches issue_011 design):
  Stage 1 (no ML): keep off-target candidates with real, control-clean indel
    signal — is_target==0, min_mm<=MAX_MM, indel_fraction in (MIN_IFRAC,1], and
    control_indel_fraction<MAX_CONTROL. (homology + control live HERE, not in the
    model, so the model can't shortcut on the on/off axis.)
  Stage 2 (learned): walk the CRAM, gate on spanning coverage, compute count-free
    pileup shape, apply the depth-augmented shape ranker -> P(edit)-like score.

Verdicts:
  LIKELY EDIT           score >= HI, clonal
  POSSIBLE — review     MID <= score < HI
  ARTIFACT (shape)      score < MID  (scattered / low-MAPQ / not clonal)
  INSUFFICIENT COVERAGE < min-span spanning reads (cannot call — NOT artifact)
  NO CRAM               sample has no CRAM in the list

This RANKS candidate sites to review; it is not an off-target recall estimate
(that needs GUIDE-seq/CIRCLE-seq). On WGS, low-VAF off-targets can sit below the
depth detection floor and surface as INSUFFICIENT COVERAGE, by design.
"""
import os, sys, argparse
import numpy as np
import pandas as pd
import pysam
import joblib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from features import read_records, features_from_records, parse_target_info, MODEL_FEATURES, check_sklearn_version

REPO = "/storage2/fs1/dspencer/Active/clinseq/projects/scge"
CRAM_LIST = f"{REPO}/cram_list.txt"
REF = "/storage2/fs1/dspencer/Active/spencerlab/rotating_students/btyler/seqs/hg38_mgi_patch.fa"
MODEL = f"{os.path.dirname(os.path.abspath(__file__))}/middle_model.pkl"
HI, MID = 0.60, 0.25


def cram_index(path):
    """Map sample_name -> cram path. Accepts either a 1-column list of cram paths
    (key = basename without .cram) or a 2-column 'sample_name<TAB>path' map (e.g.
    pileup/wgs_cram_map.tsv, needed when table sample_name != cram basename)."""
    m = {}
    for line in open(path):
        line = line.rstrip("\n")
        if not line.strip():
            continue
        parts = line.split("\t") if "\t" in line else line.split()
        if len(parts) >= 2:
            m[parts[0]] = parts[1]
        else:
            m[os.path.basename(parts[0]).replace(".cram", "")] = parts[0]
    return m


def load_table(path, sample=None, sheet=None):
    """Load ECS csv(.gz) or WGS xlsx/xlsm into a common frame with a _truth col."""
    if path.endswith((".xlsx", ".xlsm")):
        xl = pd.ExcelFile(path)
        sh = sheet or ("gold_wgs" if "gold_wgs" in xl.sheet_names else xl.sheet_names[0])
        df = xl.parse(sh)
    else:
        df = pd.read_csv(path, low_memory=False)
    if "sample_name" not in df.columns:
        df["sample_name"] = sample or os.path.basename(path).split(".")[0]
    # unify a truth column across cohorts (present only for validation)
    truth = None
    for c in ("manual_review", "pass"):
        if c in df.columns:
            truth = c; break
    df["_truth"] = (df[truth].astype(str).str.strip().isin(["1", "1.0", "1?", "yes", "True", "true", "pass"])
                    if truth else False)
    df["_if"] = pd.to_numeric(df["indel_fraction"], errors="coerce")
    df["_cif"] = pd.to_numeric(df.get("control_indel_fraction"), errors="coerce").fillna(0.0)
    df["_tgt"] = df["is_target"].astype(str).str.strip().isin(["1", "1.0", "True", "true"])
    hom = df["target_info"].apply(parse_target_info)
    df["_min_mm"] = [h["min_mm"] if h["min_mm"] is not None else 99 for h in hom]
    return df


def stage1(df, max_mm, min_ifrac, max_control, include_ontarget, no_homology_gate=False):
    """Stage-1 candidate filter. The homology gate (min_mm<=max_mm) is what ties
    discovery to a predicted candidate list; --no-homology-gate drops it so a
    somatic, control-clean indel qualifies regardless of guide homology (WGS doing
    GUIDE-seq's nomination job). min_mm survives as a ranking annotation, not a gate."""
    m = (df["_if"] > min_ifrac) & (df["_if"] <= 1.0) & (df["_cif"] < max_control)
    if not no_homology_gate:
        m &= (df["_min_mm"] <= max_mm)
    if not include_ontarget:
        m &= ~df["_tgt"]
    return df[m].copy()


def add_recurrence(res):
    """Cross-sample recurrence: a real off-target is GUIDE-specific; a site recurring
    across many unrelated guides is germline/mapping, whatever its pileup shape. The
    per-locus model cannot see this — it's a cohort-level signal added here."""
    res = res.copy()
    res["site"] = res["chrom"].astype(str) + ":" + res["start"].astype(str)
    res["guide"] = res["sample"].str.replace(r"_[12]$", "", regex=True)
    res["n_samples_at_site"] = res.groupby("site")["sample"].transform("nunique")
    res["n_guides_at_site"] = res.groupby("site")["guide"].transform("nunique")

    def tier(r):
        if r["n_guides_at_site"] >= 3:
            return "C-recurrent-artifact"
        if r["verdict"] == "LIKELY EDIT":
            return "A-review-first"
        if str(r["verdict"]).startswith("POSSIBLE"):
            return "B-review"
        return "D-artifact"
    res["priority"] = res.apply(tier, axis=1)
    return res


# ── normal-subtraction verification (automates the IGV manual-review decision) ──
# A real somatic edit is present in the EDITED sample and absent in the matched
# NORMAL. The table's precomputed control column is measured at the ANNOTATED
# position and misses indels offset by a few bp; here we recompute in the normal
# CRAM at the OBSERVED modal indel position — exactly the check a human does in IGV.
NORMAL_MAX_IF = 0.03      # control indel_frac above this => germline/artifact, not edit
MIN_EDIT_RATIO = 3.0      # edited must exceed control by this ratio
# NOTE: modal_len is NOT gated. A 1bp indel is the MOST common Cas9 outcome (the
# confirmed PLCB2 off-target is 1bp) — length is reported as an annotation only.
# The matched-normal recompute is the real edit-vs-artifact arbiter here.


def normal_cram_path(tumor_path):
    """Matched normal for a dragen '..._tumor.cram' is '....cram' in the same dir."""
    if tumor_path.endswith("_tumor.cram"):
        cand = tumor_path[:-len("_tumor.cram")] + ".cram"
        return cand if os.path.exists(cand) else None
    return None


def control_check(normal_path, chrom, pos, ref, min_span=4, bam_cache=None):
    """Recompute indel_frac in the matched normal at the OBSERVED indel position.
    Returns (control_if, control_spanning) or (None, 0) if no normal / no coverage.
    Pass bam_cache (dict) to reuse open normal handles across candidates."""
    if normal_path is None:
        return None, 0
    if bam_cache is not None and normal_path in bam_cache:
        bam = bam_cache[normal_path]
    else:
        bam = pysam.AlignmentFile(normal_path, "rc", reference_filename=ref)
        if bam_cache is not None:
            bam_cache[normal_path] = bam
    rr = read_records(bam, str(chrom), int(pos) - 1, int(pos) + 1)
    if rr is None:
        return None, 0
    f = features_from_records(rr[0], min_span=min_span, return_lowcov=True)
    if f is None or f.get("lowcov"):
        return None, f.get("spanning", 0) if f else 0
    return float(f["indel_frac"]), int(f["spanning"])


def verdict(score, feats, ctrl_if=None):
    if feats is None or feats.get("lowcov"):
        return "INSUFFICIENT COVERAGE", np.nan
    if feats["modal_mapq"] < 20 and feats["indel_frac"] > 0:
        return "ARTIFACT (low-MAPQ repeat)", score
    # normal-subtraction: present in matched control => germline/artifact
    if ctrl_if is not None:
        if ctrl_if >= NORMAL_MAX_IF or feats["indel_frac"] < MIN_EDIT_RATIO * ctrl_if:
            return "GERMLINE/ARTIFACT (in normal)", score
    if score >= HI and feats["conc_ratio"] >= 0.5:
        return "LIKELY EDIT", score
    if score >= MID:
        return "POSSIBLE — review", score
    return "ARTIFACT (shape)", score


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--sample")
    ap.add_argument("--sheet")
    ap.add_argument("--cram-list", default=CRAM_LIST)
    ap.add_argument("--ref", default=REF)
    ap.add_argument("--model", default=MODEL)
    ap.add_argument("--max-mismatch", type=int, default=2)
    ap.add_argument("--min-ifrac", type=float, default=0.05)
    ap.add_argument("--max-control", type=float, default=0.02)
    ap.add_argument("--min-span", type=int, default=8)
    ap.add_argument("--include-ontarget", action="store_true")
    ap.add_argument("--no-homology-gate", action="store_true",
                    help="drop the min_mm<=max_mismatch gate (homology-free discovery); "
                         "min_mm kept as a ranking annotation only")
    ap.add_argument("--no-normal-check", action="store_true",
                    help="skip the matched-normal indel recompute (edited-vs-control)")
    ap.add_argument("--snapshot-dir",
                    help="render an IGV-style pileup PNG for each surviving LIKELY-EDIT "
                         "candidate (edited vs matched normal) into this dir")
    ap.add_argument("--top", type=int, default=40)
    ap.add_argument("--out")
    ap.add_argument("--validate", action="store_true",
                    help="report where known edits (_truth) land")
    args = ap.parse_args()

    bundle = joblib.load(args.model)
    model, feat_names = bundle["model"], bundle["features"]
    check_sklearn_version(model, name=os.path.basename(args.model))
    crams = cram_index(args.cram_list)
    df = load_table(args.table, args.sample, args.sheet)
    cand = stage1(df, args.max_mismatch, args.min_ifrac, args.max_control,
                  args.include_ontarget, args.no_homology_gate)
    gate = "NO homology gate" if args.no_homology_gate else f"min_mm<={args.max_mismatch}"
    print(f"Stage-1: {len(df)} rows -> {len(cand)} candidates "
          f"({gate}, if>{args.min_ifrac}, control<{args.max_control}"
          f"{', off-target only' if not args.include_ontarget else ''})", flush=True)

    bam_cache, normal_cache, rows = {}, {}, []
    for _, r in cand.iterrows():
        s = r["sample_name"]
        rec = {"sample": s, "chrom": r["chrom"], "start": int(r["start"]),
               "min_mm": int(r["_min_mm"]), "tbl_if": float(r["_if"]),
               "truth": bool(r["_truth"]), "is_target": int(bool(r["_tgt"]))}
        if s not in crams:
            rec.update(score=np.nan, verdict="NO CRAM", spanning=np.nan,
                       conc_ratio=np.nan, indel_frac=np.nan); rows.append(rec); continue
        if s not in bam_cache:
            bam_cache[s] = pysam.AlignmentFile(crams[s], "rc", reference_filename=args.ref)
        rr = read_records(bam_cache[s], str(r["chrom"]), int(r["start"]), int(r["end"]))
        feats = None if rr is None else features_from_records(
            rr[0], min_span=args.min_span, return_lowcov=True)
        if feats is None or feats.get("lowcov"):
            v, sc = verdict(0.0, feats)
            rec.update(score=np.nan, verdict=v, spanning=(feats or {}).get("spanning", 0),
                       conc_ratio=np.nan, indel_frac=np.nan, modal_len=np.nan,
                       modal_pos=np.nan, ctrl_if=np.nan)
        else:
            sc = float(model.predict_proba(pd.DataFrame([{k: feats[k] for k in feat_names}]))[0, 1])
            # normal-subtraction at the OBSERVED indel position (automates IGV review)
            opos = feats.get("modal_pos") or int(r["start"])
            ctrl_if, ctrl_span = (None, 0)
            if not args.no_normal_check:
                ctrl_if, ctrl_span = control_check(normal_cram_path(crams[s]),
                                                   r["chrom"], opos, args.ref,
                                                   bam_cache=normal_cache)
            v, sc = verdict(sc, feats, ctrl_if)
            rec.update(score=sc, verdict=v, spanning=feats["spanning"],
                       conc_ratio=round(feats["conc_ratio"], 3),
                       indel_frac=round(feats["indel_frac"], 3),
                       modal_len=feats.get("modal_len"), modal_pos=opos,
                       ctrl_if=(round(ctrl_if, 3) if ctrl_if is not None else np.nan),
                       ctrl_span=ctrl_span)
        rows.append(rec)

    res = add_recurrence(pd.DataFrame(rows))
    # tier first (A>B>C>D), then score within tier
    torder = {"A-review-first": 0, "B-review": 1, "C-recurrent-artifact": 2, "D-artifact": 3}
    res["_o"] = res["priority"].map(torder).fillna(3)
    res = res.sort_values(["_o", "score", "min_mm"], ascending=[True, False, True],
                          na_position="last").drop(columns="_o").reset_index(drop=True)
    res.insert(0, "rank", res.index + 1)

    cols = ["rank", "priority", "sample", "chrom", "start", "min_mm", "indel_frac",
            "ctrl_if", "modal_len", "conc_ratio", "score", "verdict", "n_guides_at_site"]
    cols = [c for c in cols if c in res.columns]
    if res["truth"].any():
        cols.append("truth")
    print(f"\nTop {min(args.top, len(res))} candidates (tiered):")
    with pd.option_context("display.max_rows", None, "display.width", 190):
        print(res[cols].head(args.top).to_string(index=False))
    print("\ntiers: " + res["priority"].value_counts().reindex(
        ["A-review-first", "B-review", "C-recurrent-artifact", "D-artifact"]).dropna().to_string())

    if args.out:
        res.to_csv(args.out, index=False)
        print(f"\nwrote {args.out} ({len(res)} scored)")

    # auto-render IGV-style edited-vs-normal snapshots for surviving LIKELY EDITs
    if args.snapshot_dir:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from pileup_snapshot import snapshot
        os.makedirs(args.snapshot_dir, exist_ok=True)
        surv = res[res["verdict"] == "LIKELY EDIT"]
        print(f"\nrendering {len(surv)} LIKELY-EDIT snapshots -> {args.snapshot_dir}/")
        for _, k in surv.iterrows():
            s = k["sample"]; tpath = crams.get(s)
            if not tpath:
                continue
            pos = int(k["modal_pos"]) if pd.notna(k.get("modal_pos")) else int(k["start"])
            npath = normal_cram_path(tpath)
            fig, axes = plt.subplots(1, 2 if npath else 1, figsize=(15 if npath else 8, 6.5),
                                     squeeze=False)
            snapshot(axes[0][0], tpath, str(k["chrom"]), pos, args.ref, window=45,
                     title=f"{s}_tumor (EDITED)  score={k['score']:.2f}")
            if npath:
                snapshot(axes[0][1], npath, str(k["chrom"]), pos, args.ref, window=45,
                         title=f"{s} (NORMAL)  ctrl_if={k.get('ctrl_if')}")
            fig.suptitle(f"{s}  {k['chrom']}:{pos:,}  min_mm={k['min_mm']}  {k['verdict']}",
                         fontweight="bold")
            fig.tight_layout()
            fig.savefig(f"{args.snapshot_dir}/rank{int(k['rank']):03d}_{s}_{k['chrom']}_{pos}.png", dpi=120)
            plt.close(fig)
        print(f"  done ({len(surv)} PNGs)")

    if args.validate and res["truth"].any():
        kn = res[res["truth"]]
        print(f"\n== VALIDATION: {len(kn)} known edit(s) among candidates ==")
        for _, k in kn.iterrows():
            print(f"  rank {int(k['rank'])}/{len(res)}  {k['sample']} {k['chrom']}:{k['start']}"
                  f"  is_target={k['is_target']}  score={k['score']}  verdict='{k['verdict']}'")

    print("\nNOTE: ranks sites to REVIEW; not an off-target recall estimate. Low-VAF "
          "off-targets below the depth floor surface as INSUFFICIENT COVERAGE.")


if __name__ == "__main__":
    main()
