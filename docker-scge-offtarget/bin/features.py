#!/usr/bin/env python3
"""
pileup/features.py — count-free, coverage-invariant pileup featurizer.

This is the featurizer for the read-based CRISPR-edit-vs-artifact model. It is a
cleaned successor to results/diag/pileup_probe.py:locus_features with two design
rules imposed for cross-assay (ECS deep-amplicon -> WGS ~30x) portability:

  1. FEATURES ARE COUNT-FREE. Every model feature is a fraction, a ratio, a bp
     spread, or a MAPQ — nothing that scales with sequencing depth. Raw read
     counts (`spanning`) are returned ONLY as a coverage-adequacy GATE, never as
     a model input. (`n_distinct_pos` from the old probe is dropped: it grows
     with depth and does not transfer.)

  2. DEPTH IS SIMULABLE. read_records() walks the CRAM once and returns per-read
     records; features_from_records() can then recompute features at full depth
     OR after downsampling to a target depth, WITHOUT re-fetching. This is the
     ECS->WGS bridge: an ECS locus with 500 reads can be evaluated as if it had
     30, so the model sees the noisy low-depth concordance it will face on WGS.

Concordance is POSITION-centric: a CRISPR cut puts indels at the SAME reference
position even when their lengths differ read-to-read; repeat slippage scatters
them. So the clonal signal is keyed on WHERE the indel sits, not its size.
"""
import re
import random
from collections import Counter
import numpy as np

WIN = 25          # bp around the locus to search for indel events
FETCH_PAD = 200   # bp padding for read fetch (reads can start before the locus)
MIN_SPAN = 8      # default min spanning reads to score a locus (coverage gate)


# ── guide-homology parsing (assay-independent; from Cas-OFFinder/CRISPRme etc.) ─
def parse_target_info(ti):
    """Extract guide-homology annotation from the `target_info` string.

    Returns dict(min_mm, n_tools, pam). min_mm = fewest protospacer mismatches
    across the tools that predicted this site (None if unannotated); n_tools =
    number of predicting tools (SRC=...); pam = the PAM if present.
    """
    s = "" if ti is None else str(ti)
    mm = None
    m = re.search(r"MISMATCHES=([0-9,]+)", s)
    if m:
        vals = [int(x) for x in m.group(1).split(",") if x != ""]
        mm = min(vals) if vals else None
    n_tools = 0
    m = re.search(r"SRC=([^;]+)", s)
    if m:
        n_tools = len([t for t in re.split(r"[|,]", m.group(1)) if t.strip()])
    pam = None
    m = re.search(r"PAM=([ACGTN]+)", s)
    if m:
        pam = m.group(1)
    return {"min_mm": mm, "n_tools": n_tools, "pam": pam}


# ── one CRAM walk -> per-read records ─────────────────────────────────────────
def read_records(bam, chrom, start, end, win=WIN, fetch_pad=FETCH_PAD):
    """Walk the reads over a candidate locus ONCE; return (records, mid) or None.

    records: list of dicts, one per usable read:
        {'spans': bool,                # covers the locus midpoint
         'indel': (ref_pos, length) or None,   # largest I/D in-window
         'mapq': int,
         'softclip': [ref_pos, ...]}   # in-window soft-clip ref edges
    Downstream featurizers subsample this list to simulate lower depth.
    """
    lo, hi = start - win, end + win
    fetch_lo = max(0, start - fetch_pad)
    mid = (start + end) // 2
    try:
        reads = bam.fetch(chrom, fetch_lo, end + fetch_pad)
    except (ValueError, KeyError):
        return None
    records = []
    for r in reads:
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
            continue
        spans = r.reference_start <= mid <= (r.reference_end or r.reference_start)
        refpos = r.reference_start
        best = None            # (abs_ref_pos, length) of largest in-window indel
        sc = []
        for op, ln in (r.cigartuples or []):
            if op in (0, 7, 8):       # M/=/X consume ref+query
                refpos += ln
            elif op == 2:             # D consumes ref
                if lo <= refpos <= hi and (best is None or ln > best[1]):
                    best = (refpos, ln)
                refpos += ln
            elif op == 1:             # I consumes query only
                if lo <= refpos <= hi and (best is None or ln > best[1]):
                    best = (refpos, ln)
            elif op == 4:             # S soft-clip; record ref edge
                if lo <= refpos <= hi:
                    sc.append(refpos)
        records.append({"spans": spans, "indel": best,
                        "mapq": r.mapping_quality, "softclip": sc})
    return records, mid


def _subsample(records, downsample_to, rng):
    """Thin records to ~downsample_to spanning reads, modelling lower depth.

    Keeps each read independently with prob downsample_to/n_spanning, so both
    spanning and indel reads are thinned together and fractions are preserved in
    expectation while acquiring realistic low-depth sampling noise.
    """
    n_span = sum(1 for r in records if r["spans"])
    if downsample_to is None or n_span <= downsample_to:
        return records
    p = downsample_to / n_span
    return [r for r in records if rng.random() < p]


# ── records -> count-free feature dict ────────────────────────────────────────
def features_from_records(records, min_span=MIN_SPAN, downsample_to=None,
                          rng=None, return_lowcov=False):
    """Compute count-free features from read records (optionally downsampled).

    Returns a dict, or None (or {'spanning','lowcov':True} if return_lowcov) when
    spanning coverage < min_span. `spanning` is included for GATING only — do not
    feed it to a model.
    """
    if downsample_to is not None:
        records = _subsample(records, downsample_to, rng or random.Random(0))

    spanning = sum(1 for r in records if r["spans"])
    if spanning < min_span:
        return {"spanning": spanning, "lowcov": True} if return_lowcov else None

    indel = [r for r in records if r["indel"] is not None]
    n_indel = len(indel)
    if n_indel:
        pos = np.array([r["indel"][0] for r in indel])
        length = np.array([r["indel"][1] for r in indel])
        mapq = np.array([r["mapq"] for r in indel])
        center = Counter(pos.tolist()).most_common(1)[0][0]   # modal cut (robust)
        near = np.abs(pos - center) <= 2
        conc_ratio = float(near.sum()) / n_indel              # count-free clonality
        pos_conc = float(near.sum()) / spanning               # frac of all reads
        pos_mad = float(np.median(np.abs(pos - center)))
        modal_len = int(Counter(length[near].tolist()).most_common(1)[0][0])
        modal_mapq = float(np.median(mapq[near]))
        modal_pos = int(center)                               # OBSERVED indel ref pos
    else:
        conc_ratio = pos_conc = pos_mad = modal_mapq = 0.0
        modal_len = 0
        modal_pos = None

    sc_all = [p for r in records for p in r["softclip"]]
    sc_modal = (Counter(sc_all).most_common(1)[0][1] / spanning) if sc_all else 0.0

    return dict(
        spanning=spanning,               # GATE ONLY (not a model feature)
        lowcov=False,
        indel_frac=n_indel / spanning,   # fraction ✓
        conc_ratio=conc_ratio,           # of indel reads, frac at modal cut ✓ count-free
        pos_conc=pos_conc,               # of spanning reads, frac at modal cut ✓
        pos_mad=pos_mad,                 # bp spread of indel positions ✓
        modal_len=modal_len,             # bp size of modal indel ✓
        modal_mapq=modal_mapq,           # repeat => low MAPQ ✓
        softclip_frac=sc_modal,          # fraction ✓
        modal_pos=modal_pos,             # OBSERVED indel ref pos (for control recompute/snapshot)
    )


# ── convenience one-shot wrapper (fetch + compute) ────────────────────────────
def locus_features(bam, chrom, start, end, min_span=MIN_SPAN,
                   downsample_to=None, rng=None, return_lowcov=False):
    """Fetch reads at a locus and return count-free features (or None/lowcov)."""
    rr = read_records(bam, chrom, start, end)
    if rr is None:
        return None
    records, _ = rr
    return features_from_records(records, min_span=min_span,
                                 downsample_to=downsample_to, rng=rng,
                                 return_lowcov=return_lowcov)


# count-free features safe to feed a model (spanning/lowcov are gates, not inputs)
MODEL_FEATURES = ["indel_frac", "conc_ratio", "pos_conc", "pos_mad",
                  "modal_len", "modal_mapq", "softclip_frac"]
HOMOLOGY_FEATURES = ["min_mm", "n_tools"]


def check_sklearn_version(model, name="model", raise_on_backward=True):
    """Loudly flag a scikit-learn version skew between a pickled model and the
    runtime, so a container/library bump can't silently mis-score.

    No-op for non-sklearn models (e.g. xgboost) or when the pickle did not record
    a version. Behaviour on a (major, minor) mismatch:
      * runtime OLDER than the pickle  -> raise RuntimeError (the backward-
        incompatible load is the case that silently corrupts predictions);
      * any other skew (e.g. the intended 1.6.1 model under a 1.8.0 runtime)
        -> RuntimeWarning on stderr and proceed.
    """
    import warnings
    pickled = getattr(model, "_sklearn_version", None)
    if pickled is None:  # unwrap a Pipeline and inspect its final estimator
        steps = getattr(model, "steps", None)
        if steps:
            try:
                pickled = getattr(steps[-1][1], "_sklearn_version", None)
            except Exception:
                pickled = None
    if pickled is None:
        return  # not an sklearn estimator, or version not recorded in the pickle
    try:
        import sklearn
        runtime = sklearn.__version__
        r = tuple(int(x) for x in runtime.split(".")[:2])
        p = tuple(int(x) for x in pickled.split(".")[:2])
    except Exception:
        return
    if r == p:
        return
    msg = (f"[sklearn version guard] {name}: model pickled under scikit-learn "
           f"{pickled} but the runtime is {runtime}.")
    if raise_on_backward and r < p:
        raise RuntimeError(
            msg + " The runtime is OLDER than the model — loading it can silently "
            f"produce wrong scores. Rebuild the image with scikit-learn>={pickled} "
            "or re-pickle the model.")
    warnings.warn(
        msg + " Proceeding with a forward-compatible load; sanity-check outputs if "
        "this skew is unexpected.", RuntimeWarning)
