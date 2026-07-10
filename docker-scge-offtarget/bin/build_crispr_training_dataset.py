#!/usr/bin/env python3
"""
build_crispr_training_dataset.py

Extract per-read features and site-level labels for CRISPR training data.

For each (edited_cram, control_cram, target_file) pair in the mastersheet:
  1. Parse target positions from VCF.gz or CSV target files
  2. Compute per-site indel fractions in edited and control BAMs
  3. Assign site label: 1 if edited_frac >= min_indel_frac AND ctrl_frac <= max_ctrl_frac
                        0 if at low-indel positions in edited BAM (background)
  4. Extract the same 16 features used by extract_variant_reads_ML.py
  5. Sample background negatives from random genomic windows in control BAM

Output: parquet with columns: 14 features + label, sample_id, dataset,
        individual_id, chrom, pos

Usage:
  python build_crispr_training_dataset.py \\
    --mastersheet assets/all_samples_mastersheet.csv \\
    --output-file training_data/crispr_training_v2.parquet \\
    --window 150 \\
    --min-indel-fraction 0.05 \\
    --max-control-fraction 0.01 \\
    --n-background-sites 5000 \\
    --threads 8
"""

import argparse
import importlib.util
import os
import sys
import subprocess
import random
import re
import concurrent.futures
import time
import traceback

# Add repo root to sys.path so locally-installed packages are found
_repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

# Ensure required packages are available; install to repo root if missing.
# scipy is checked via scipy.stats (the submodule that loads the broken C extension)
# so a broken conda scipy is caught here rather than at first use.
for _import_check, _pip_name in [("pyarrow", "pyarrow"),
                                   ("biotite", "biotite"),
                                   ("pyranges", "pyranges"),
                                   ("joblib", "joblib"),
                                   ("scipy.stats", "scipy")]:
    try:
        __import__(_import_check)
    except ImportError:
        # Evict any partially-loaded version from the module cache
        _base = _import_check.split(".")[0]
        for _k in list(sys.modules.keys()):
            if _k == _base or _k.startswith(_base + "."):
                del sys.modules[_k]
        subprocess.run([sys.executable, "-m", "pip", "install", "-q",
                        "--target", _repo_root, _pip_name], check=True)

import numpy as np
import pandas as pd
import pysam

# ---------------------------------------------------------------------------
# Import feature functions from extract_variant_reads_ML.py
# ---------------------------------------------------------------------------

def _load_ml_module():
    """Load extract_variant_reads_ML.py from the same bin/ directory."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    ml_path = os.path.join(script_dir, "extract_variant_reads_ML.py")
    if not os.path.exists(ml_path):
        raise FileNotFoundError(f"Cannot find extract_variant_reads_ML.py at {ml_path}")
    spec = importlib.util.spec_from_file_location("extract_variant_reads_ML", ml_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

_ML = _load_ml_module()

parse_cigar_once             = _ML.parse_cigar_once
count_mismatches_fast        = _ML.count_mismatches_fast
calculate_control_fractions  = _ML.calculate_control_fractions
calculate_exclusivity_features = _ML.calculate_exclusivity_features

# ---------------------------------------------------------------------------
# Import shared feature helpers from crispr_ml_features.py
# ---------------------------------------------------------------------------

def _load_features_module():
    """Load crispr_ml_features.py from the same bin/ directory."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    feat_path = os.path.join(script_dir, "crispr_ml_features.py")
    if not os.path.exists(feat_path):
        raise FileNotFoundError(
            f"Cannot find crispr_ml_features.py at {feat_path}"
        )
    spec = importlib.util.spec_from_file_location("crispr_ml_features", feat_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

_FEAT = _load_features_module()

_softclip_pam_distance        = _FEAT._softclip_pam_distance
_precompute_site_ref_features = _FEAT._precompute_site_ref_features

FEATURE_NAMES = [
    'read_pair_gap', 'read_softclip',
    'read_del_vs_control', 'read_ins_vs_control', 'read_mismatch_vs_control',
    'deletion_exclusive_to_edited',
    'is_at_any_target_site', 'distance_to_closest_pam',
    'total_indel_size', 'indel_size_category', 'indel_complexity_score',
    'softclip_dist_to_PAM', 'is_in_homopolymer', 'dist_to_microsatellite',
]

# Reference FASTA (same as used by the pipeline)
DEFAULT_FASTA = ("/storage2/fs1/dspencer/Active/clinseq/projects/scge/data/refdata/"
                 "singh_v4.3.6/hg38_PLVM_CD19_CARv4_cd34.fa")

# Autosomes + chrX for background sampling
BACKGROUND_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

# Minimum gap (bp) between a background window and any target site
MIN_GAP_FROM_TARGET = 500_000


# ---------------------------------------------------------------------------
# Target file parsing
# ---------------------------------------------------------------------------

def load_targets(target_file: str) -> pd.DataFrame:
    """
    Parse target sites from a VCF.gz or CSV file.

    Returns a DataFrame with at minimum: Chromosome, Start (0-based), Ontarget.
    """
    if target_file.endswith(".vcf.gz") or target_file.endswith(".vcf"):
        return _load_vcf_targets(target_file)
    elif target_file.endswith(".csv"):
        return _load_csv_targets(target_file)
    else:
        raise ValueError(f"Unsupported target file format: {target_file}")


def _load_vcf_targets(vcf_path: str) -> pd.DataFrame:
    rows = []
    vcf_in = pysam.VariantFile(vcf_path)
    for rec in vcf_in:
        info_dict = dict(rec.info)
        rows.append({
            'Chromosome': rec.chrom,
            'Start': rec.pos - 1,   # 0-based
            'Pos': rec.pos,
            'Ontarget': 1 if info_dict.get('TARGET', False) else 0,
        })
    return pd.DataFrame(rows)


def _load_csv_targets(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    # Rename Start column if necessary; CSV uses 1-based 'Start'
    if 'Start' not in df.columns and 'Chromosome' in df.columns:
        raise ValueError(f"CSV target file missing 'Start' column: {csv_path}")
    df = df.rename(columns={'Start': 'Start_1based'})
    df['Start'] = df['Start_1based'] - 1   # convert to 0-based
    df['Pos'] = df['Start_1based']
    if 'On_target' in df.columns:
        df['Ontarget'] = df['On_target']
    elif 'Ontarget' not in df.columns:
        df['Ontarget'] = 0
    return df[['Chromosome', 'Start', 'Pos', 'Ontarget']]


# ---------------------------------------------------------------------------
# Site-level indel fraction calculation
# ---------------------------------------------------------------------------

def compute_site_indel_fraction(bam: pysam.AlignmentFile,
                                 chrom: str, pos: int, window: int = 150,
                                 min_mapq: int = 1):
    """
    Count indel-bearing reads and total reads in [pos-window, pos+window].
    Returns (indel_fraction, total_reads).
    """
    start = max(0, pos - window)
    end = pos + window
    total = indel_count = 0
    for read in bam.fetch(chrom, start, end):
        if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue
        total += 1
        if read.cigartuples and any(op in (1, 2) for op, _ in read.cigartuples):
            indel_count += 1
    frac = indel_count / total if total > 0 else 0.0
    return frac, total


def compute_site_labels(edited_bam: pysam.AlignmentFile,
                         control_bam: pysam.AlignmentFile,
                         targets_df: pd.DataFrame,
                         window: int,
                         min_indel_frac: float,
                         max_ctrl_frac: float) -> dict:
    """
    Returns dict mapping (chrom, pos_0based) -> (label, is_on_target)
      label=1  : edited_frac >= min_indel_frac AND ctrl_frac <= max_ctrl_frac
      label=0  : edited_frac <= 0.01  (clear background in edited sample)
      label=-1 : ambiguous, skip
    """
    labels = {}
    for _, row in targets_df.iterrows():
        chrom = row['Chromosome']
        pos = int(row['Start'])
        is_on_target = int(row.get('Ontarget', 0))

        edited_frac, _ = compute_site_indel_fraction(edited_bam, chrom, pos, window)
        ctrl_frac, _   = compute_site_indel_fraction(control_bam, chrom, pos, window)

        if edited_frac >= min_indel_frac and ctrl_frac <= max_ctrl_frac:
            label = 1
        elif edited_frac <= 0.01:
            label = 0
        else:
            label = -1   # ambiguous; skip

        labels[(chrom, pos)] = (label, is_on_target)
    return labels


# ---------------------------------------------------------------------------
# Per-read feature extraction
# ---------------------------------------------------------------------------

def _fetch_ctrl_pool(control_bam, chrom, pos, window=100):
    """
    Fetch control reads in [pos-window, pos+window] once and return lightweight
    tuples (start, end, has_ins, has_del).  Call once per site; pass the result
    to _exclusivity_from_pool for each read to avoid per-read CRAM fetches.
    """
    pool = []
    start = max(0, pos - window)
    for cr in control_bam.fetch(chrom, start, pos + window):
        if cr.is_unmapped or cr.is_duplicate or not cr.cigartuples:
            continue
        pool.append((
            cr.reference_start,
            cr.reference_end,
            any(op == 1 for op, _ in cr.cigartuples),
            any(op == 2 for op, _ in cr.cigartuples),
        ))
    return pool


def _exclusivity_from_pool(read, ctrl_pool):
    """Compute exclusivity features against a pre-fetched control-read pool."""
    result = {'deletion_exclusive_to_edited': 0, 'control_has_same_variant': 0}
    if not read.cigartuples:
        return result
    has_ins = any(op == 1 for op, _ in read.cigartuples)
    has_del = any(op == 2 for op, _ in read.cigartuples)
    r_start = read.reference_start
    r_end   = read.reference_end
    ctrl_has_ins = ctrl_has_del = False
    for ctrl_start, ctrl_end, c_has_ins, c_has_del in ctrl_pool:
        if r_end >= ctrl_start and ctrl_end >= r_start:
            if has_ins and c_has_ins:
                ctrl_has_ins = True
            if has_del and c_has_del:
                ctrl_has_del = True
    result['deletion_exclusive_to_edited'] = 1 if (has_del and not ctrl_has_del) else 0
    result['control_has_same_variant']     = 1 if (ctrl_has_ins or ctrl_has_del) else 0
    return result


def extract_read_features(read, control_bam, chrom, pos, pam_positions, is_on_target,
                           _ctrl_fracs=None, _ctrl_pool=None, _ref_site_feats=None):
    """
    Extract the 14 CRISPR read features for a single pysam.AlignedSegment.
    Returns a dict or None if read is unusable.

    _ctrl_fracs and _ctrl_pool are optional pre-computed site-level values.
    When provided (from extract_features_for_site), no CRAM fetch is needed per
    read — the control data is reused across all reads at the same site.
    When absent (background sampling, where every pos is unique), falls back to
    live CRAM calls as before.

    _ref_site_feats is the dict returned by _precompute_site_ref_features()
    (keys: is_in_homopolymer, dist_to_microsatellite). None triggers fallback
    sentinel values (0 and 999 respectively).
    """
    if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary:
        return None
    if read.cigartuples is None:
        return None

    cigar_data = parse_cigar_once(read)
    mismatch   = count_mismatches_fast(read)

    read_start = read.reference_start
    read_end   = read.reference_end or read_start

    # is_at_any_target_site: within 25 bp of any PAM position
    is_at_any = 0
    if pam_positions:
        for p in pam_positions:
            if read_start <= p + 25 and read_end >= p - 25:
                is_at_any = 1
                break

    if _ctrl_fracs is None:
        ctrl_fracs = calculate_control_fractions(control_bam, chrom, pos, window=50)
    else:
        ctrl_fracs = _ctrl_fracs

    if _ctrl_pool is None:
        excl = calculate_exclusivity_features(read, control_bam, chrom, pos, window=100)
    else:
        excl = _exclusivity_from_pool(read, _ctrl_pool)

    has_del = 1 if cigar_data['deletions'] > 0 else 0
    has_ins = 1 if cigar_data['insertions'] > 0 else 0
    has_mm  = 1 if mismatch > 0 else 0

    total_indel = cigar_data['insertions'] + cigar_data['deletions']

    if total_indel == 0:
        indel_cat = 0
    elif total_indel <= 3:
        indel_cat = 1
    elif total_indel <= 10:
        indel_cat = 2
    else:
        indel_cat = 3

    indel_ops = [op for op, _ in read.cigartuples if op in (1, 2)]
    complexity = min(len(indel_ops) + (total_indel / 10.0), 10.0)

    dist_to_pam = min(abs(p - pos) for p in pam_positions) if pam_positions else -1

    ref_feats = _ref_site_feats or {}

    return {
        'read_pair_gap': abs(read.template_length) if (
            hasattr(read, 'template_length') and read.is_paired and read.is_proper_pair
        ) else -1,
        'read_softclip':  cigar_data['softclips'],
        'read_del_vs_control':      has_del - ctrl_fracs['fraction_control_reads_del'],
        'read_ins_vs_control':      has_ins - ctrl_fracs['fraction_control_reads_ins'],
        'read_mismatch_vs_control': has_mm  - ctrl_fracs['fraction_control_reads_mismatch'],
        'deletion_exclusive_to_edited': excl['deletion_exclusive_to_edited'],
        'is_on_target_site':      is_on_target,
        'is_at_any_target_site':  is_at_any,
        'distance_to_closest_pam': dist_to_pam,
        'total_indel_size':        total_indel,
        'indel_size_category':     indel_cat,
        'indel_complexity_score':  complexity,
        'softclip_dist_to_PAM':    _softclip_pam_distance(read, pam_positions),
        'is_in_homopolymer':       ref_feats.get('is_in_homopolymer', 0),
        'dist_to_microsatellite':  ref_feats.get('dist_to_microsatellite', 999),
    }


def extract_features_for_site(edited_bam, control_bam, chrom, pos,
                                pam_positions, label, is_on_target, window,
                                fasta_file=None):
    """Extract features for all reads at a target site.

    Control-BAM data is fetched once per site (not once per read) using
    pre-computed ctrl_fracs and ctrl_pool, then reused for every read in the
    window.  For sites with hundreds of reads this eliminates ~N_reads CRAM
    round-trips and replaces them with 2 fetches regardless of read depth.

    fasta_file: optional open pysam.FastaFile handle.  When provided,
    reference-sequence site features (is_in_homopolymer, dist_to_microsatellite)
    are pre-computed once and passed to every extract_read_features call.
    When None, fallback sentinel values are used (0 and 999).
    """
    # Pre-compute site-level control features once
    ctrl_fracs = calculate_control_fractions(control_bam, chrom, pos, window=50)
    ctrl_pool  = _fetch_ctrl_pool(control_bam, chrom, pos, window=100)

    # Pre-compute site-level reference features once (if FASTA available)
    ref_site_feats = None
    if fasta_file is not None:
        ref_site_feats = _precompute_site_ref_features(fasta_file, chrom, pos)

    rows = []
    start = max(0, pos - window)
    end   = pos + window
    for read in edited_bam.fetch(chrom, start, end):
        feat = extract_read_features(
            read, control_bam, chrom, pos, pam_positions, is_on_target,
            _ctrl_fracs=ctrl_fracs, _ctrl_pool=ctrl_pool,
            _ref_site_feats=ref_site_feats,
        )
        if feat is not None:
            feat['label'] = label
            feat['chrom'] = chrom
            feat['pos']   = pos
            rows.append(feat)
    return rows


# ---------------------------------------------------------------------------
# Background negative sampling
# ---------------------------------------------------------------------------

def _chrom_lengths(bam: pysam.AlignmentFile) -> dict:
    return {sq['SN']: sq['LN'] for sq in bam.header['SQ'] if sq['SN'] in BACKGROUND_CHROMS}


def sample_background_negatives(iter_bam: pysam.AlignmentFile,
                                  control_bam: pysam.AlignmentFile,
                                  target_positions: list,
                                  n: int,
                                  window: int,
                                  rng: random.Random) -> list:
    """
    Sample n reads from random genomic windows far from any target.
    All resulting rows get label=0.

    iter_bam   -- dedicated handle used ONLY for the outer fetch loop
    control_bam -- handle passed to extract_read_features for inner feature calls
    Using separate handles avoids nested fetches on the same CRAM file descriptor,
    which corrupts the shared cram_fd state and causes CRC32 failures.
    """
    chrom_lens = _chrom_lengths(iter_bam)
    if not chrom_lens:
        return []

    # Build fast exclusion set: {chrom: sorted list of positions}
    excl = {}
    for (chrom, pos) in target_positions:
        excl.setdefault(chrom, []).append(pos)

    rows = []
    attempts = 0
    max_attempts = n * 20

    while len(rows) < n and attempts < max_attempts:
        attempts += 1
        chrom = rng.choice(list(chrom_lens.keys()))
        chrom_len = chrom_lens[chrom]
        if chrom_len < 2 * window + 1:
            continue
        center = rng.randint(window, chrom_len - window - 1)

        # Check distance from any target on this chromosome
        if chrom in excl:
            if any(abs(p - center) < MIN_GAP_FROM_TARGET for p in excl[chrom]):
                continue

        for read in iter_bam.fetch(chrom, max(0, center - window), center + window):
            if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary:
                continue
            if read.cigartuples is None:
                continue
            feat = extract_read_features(
                read, control_bam, chrom, center,
                pam_positions=None, is_on_target=0
            )
            if feat is not None:
                feat['label'] = 0
                feat['chrom'] = chrom
                feat['pos']   = center
                rows.append(feat)
                if len(rows) >= n:
                    break

    return rows[:n]


# ---------------------------------------------------------------------------
# Global target map (for cross-sample negatives)
# ---------------------------------------------------------------------------

def load_all_targets(mastersheet_df: pd.DataFrame, offtarget_csv: str = None) -> dict:
    """
    Build global target map from all samples in the mastersheet.

    Returns:
        dict mapping (chrom, pos_0based) -> {
            'sample_ids': set of sample IDs that own this site,
            'pam_pos':    representative PAM position (1-based Pos from target file),
            'is_offtarget': bool,
        }

    Off-target hotspot CSV (optional) columns: sample_id, chrom, pos[, pam].
    Off-target entries are merged with is_offtarget=True.
    """
    global_targets = {}

    for _, row in mastersheet_df.iterrows():
        sample_id   = row['id']
        target_file = row['target_file']
        if not os.path.exists(target_file):
            print(f"  WARNING: target file missing for {sample_id}: {target_file}", flush=True)
            continue
        try:
            targets_df = load_targets(target_file)
        except Exception as e:
            print(f"  WARNING: could not load targets for {sample_id}: {e}", flush=True)
            continue

        for _, trow in targets_df.iterrows():
            chrom   = trow['Chromosome']
            pos     = int(trow['Start'])        # 0-based
            pam_pos = int(trow.get('Pos', pos + 1))  # 1-based
            key = (chrom, pos)
            if key not in global_targets:
                global_targets[key] = {
                    'sample_ids': set(),
                    'pam_pos': pam_pos,
                    'is_offtarget': False,
                }
            global_targets[key]['sample_ids'].add(sample_id)

    # Merge validated off-target hotspots if provided
    if offtarget_csv and os.path.exists(offtarget_csv):
        ot_df = pd.read_csv(offtarget_csv)
        for _, row in ot_df.iterrows():
            chrom     = str(row['chrom'])
            pos       = int(row['pos'])
            sample_id = str(row['sample_id'])
            # Use explicit pam column if present, otherwise fall back to pos
            pam_pos   = int(row['pam']) if 'pam' in ot_df.columns else pos
            key = (chrom, pos)
            if key not in global_targets:
                global_targets[key] = {
                    'sample_ids': set(),
                    'pam_pos': pam_pos,
                    'is_offtarget': True,
                }
            global_targets[key]['sample_ids'].add(sample_id)
            global_targets[key]['is_offtarget'] = True

    return global_targets


def sample_cross_sample_negatives(
    edited_bam: pysam.AlignmentFile,
    ctrl_bam: pysam.AlignmentFile,
    cross_sample_sites: list,
    max_sites: int,
    min_indel_frac: float = 0.01,
    window: int = 150,
    rng=None,
) -> list:
    """
    Extract reads from other samples' target positions that are confirmed
    unedited in this sample's edited BAM.  These are hard negatives: the model
    sees is_at_any_target_site=1 with no real edit signal.

    cross_sample_sites: list of (chrom, pos_0based, pam_pos) from other samples.
    Returns rows with label=0, is_on_target_site=0.
    """
    sites = list(cross_sample_sites)
    if rng is not None:
        rng.shuffle(sites)

    rows = []
    sites_used = 0

    for chrom, pos, pam_pos in sites:
        if sites_used >= max_sites:
            break

        # Confirm unedited in this sample's edited BAM
        indel_frac, total_reads = compute_site_indel_fraction(edited_bam, chrom, pos, window)
        if indel_frac >= min_indel_frac:
            continue   # possibly edited here; skip
        if total_reads < 5:
            continue   # insufficient coverage; skip

        # Pre-compute control features once for this site
        ctrl_fracs = calculate_control_fractions(ctrl_bam, chrom, pos, window=50)
        ctrl_pool  = _fetch_ctrl_pool(ctrl_bam, chrom, pos, window=100)

        start = max(0, pos - window)
        end   = pos + window
        site_rows = []
        for read in edited_bam.fetch(chrom, start, end):
            feat = extract_read_features(
                read, ctrl_bam, chrom, pos,
                pam_positions=[pam_pos],
                is_on_target=0,
                _ctrl_fracs=ctrl_fracs,
                _ctrl_pool=ctrl_pool,
            )
            if feat is not None:
                feat['label'] = 0
                feat['chrom'] = chrom
                feat['pos']   = pos
                site_rows.append(feat)

        if site_rows:
            rows.extend(site_rows)
            sites_used += 1

    return rows


# ---------------------------------------------------------------------------
# Per-sample processing
# ---------------------------------------------------------------------------

def process_sample(row: pd.Series, args, global_targets: dict = None):
    """
    Process one mastersheet row.

    Writes results directly to a per-sample parquet chunk file inside the
    worker process and returns the chunk file path (or None if nothing was
    extracted).  This avoids shipping millions of dicts back to the main
    process over IPC, which was the cause of OOM failures.

    global_targets: dict from load_all_targets(), used for cross-sample
    negatives and off-target hotspot positives.  None disables both.
    """
    sample_id   = row['id']
    dataset     = row['dataset']
    indiv_id    = row['individual_id']
    edited_path = row['edited_cram']
    control_path= row['control_cram']
    target_path = row['target_file']

    print(f"  Processing {sample_id} ...", flush=True)

    # Validate paths
    for p in (edited_path, control_path, target_path):
        if not os.path.exists(p):
            print(f"    WARNING: missing file {p} — skipping {sample_id}", flush=True)
            return None

    try:
        edited_bam   = pysam.AlignmentFile(edited_path,  "rc", reference_filename=args.fasta)
        control_bam  = pysam.AlignmentFile(control_path, "rc", reference_filename=args.fasta)
        # Separate handle used only for outer iteration in sample_background_negatives.
        # Avoids nested fetch() calls on the same CRAM file descriptor (shared cram_fd),
        # which causes CRC32 failures when inner feature calls seek to overlapping regions.
        control_iter = pysam.AlignmentFile(control_path, "rc", reference_filename=args.fasta)
        fasta_file   = pysam.FastaFile(args.fasta)
        targets_df   = load_targets(target_path)
    except Exception as e:
        print(f"    ERROR opening files for {sample_id}: {e}", flush=True)
        return None

    if len(targets_df) == 0:
        print(f"    WARNING: no targets in {target_path} — skipping", flush=True)
        return None

    # --- Compute site labels ---
    site_labels = compute_site_labels(
        edited_bam, control_bam, targets_df,
        window=args.window,
        min_indel_frac=args.min_indel_fraction,
        max_ctrl_frac=args.max_control_fraction,
    )

    pos_sites = sum(1 for (l, _) in site_labels.values() if l == 1)
    neg_sites = sum(1 for (l, _) in site_labels.values() if l == 0)
    print(f"    Sites: {pos_sites} positive, {neg_sites} negative, "
          f"{len(site_labels)-pos_sites-neg_sites} ambiguous", flush=True)

    if pos_sites == 0:
        print(f"    WARNING: no CRISPR-positive sites found for {sample_id}", flush=True)

    # PAM positions per chrom for distance calculation
    pam_by_chrom = {}
    for _, trow in targets_df.iterrows():
        pam_by_chrom.setdefault(trow['Chromosome'], []).append(trow['Pos'])

    # --- Extract features for labeled target sites ---
    all_rows = []
    for (chrom, pos), (label, is_on_target) in site_labels.items():
        if label == -1:
            continue
        pam_positions = pam_by_chrom.get(chrom, [])
        site_rows = extract_features_for_site(
            edited_bam, control_bam, chrom, pos,
            pam_positions, label, is_on_target, args.window,
            fasta_file=fasta_file,
        )
        all_rows.extend(site_rows)

    # --- Off-target hotspot positives (owning sample only) ---
    if global_targets is not None and getattr(args, 'offtarget_csv', None):
        ot_rows = []
        for (chrom, pos), entry in global_targets.items():
            if not entry.get('is_offtarget', False):
                continue
            if sample_id not in entry['sample_ids']:
                continue
            if (chrom, pos) in site_labels:
                continue   # already handled as a labeled target site

            edited_frac, _ = compute_site_indel_fraction(edited_bam, chrom, pos, args.window)
            ctrl_frac, _   = compute_site_indel_fraction(control_bam, chrom, pos, args.window)

            if edited_frac >= args.min_indel_fraction and ctrl_frac <= args.max_control_fraction:
                ot_label = 1
            elif edited_frac <= 0.01:
                ot_label = 0
            else:
                continue   # ambiguous

            pam_pos = entry['pam_pos']
            site_rows_ot = extract_features_for_site(
                edited_bam, control_bam, chrom, pos,
                pam_positions=[pam_pos], label=ot_label,
                is_on_target=0, window=args.window,
                fasta_file=fasta_file,
            )
            ot_rows.extend(site_rows_ot)

        all_rows.extend(ot_rows)
        if ot_rows:
            n_ot_pos = sum(1 for r in ot_rows if r['label'] == 1)
            n_ot_neg = len(ot_rows) - n_ot_pos
            print(f"    Off-target hotspots: {n_ot_pos} pos, {n_ot_neg} neg rows", flush=True)

    # --- Background negatives ---
    target_positions = list(site_labels.keys())
    rng = random.Random(hash(sample_id) & 0xFFFFFFFF)
    bg_rows = sample_background_negatives(
        control_iter, control_bam, target_positions,
        n=args.n_background_sites,
        window=args.window,
        rng=rng,
    )
    all_rows.extend(bg_rows)

    # --- Cross-sample negatives ---
    if getattr(args, 'cross_sample_negatives', False) and global_targets is not None:
        cross_sites = [
            (chrom, pos, entry['pam_pos'])
            for (chrom, pos), entry in global_targets.items()
            if sample_id not in entry['sample_ids']
        ]
        cross_rows = sample_cross_sample_negatives(
            edited_bam, control_bam,
            cross_sites,
            max_sites=args.cross_sample_sites,
            min_indel_frac=0.01,
            window=args.window,
            rng=rng,
        )
        all_rows.extend(cross_rows)
        print(f"    Cross-sample negatives: {len(cross_rows)} rows from "
              f"up to {args.cross_sample_sites} sites", flush=True)

    edited_bam.close()
    control_bam.close()
    control_iter.close()
    fasta_file.close()

    if not all_rows:
        return None

    # Tag rows with sample metadata
    for r in all_rows:
        r['sample_id']    = sample_id
        r['dataset']      = dataset
        r['individual_id'] = indiv_id

    n_pos = sum(1 for r in all_rows if r['label'] == 1)
    n_neg = sum(1 for r in all_rows if r['label'] == 0)
    print(f"    Extracted {len(all_rows)} rows ({n_pos} pos, {n_neg} neg)", flush=True)

    # Write chunk parquet directly from within the worker process.
    # Returning only a file path over IPC (instead of millions of dicts)
    # keeps peak memory inside each worker and eliminates large pickling.
    # is_on_target_site is kept as metadata (for post-hoc analysis) but is NOT
    # in FEATURE_NAMES and therefore won't be used during model training.
    meta_cols = ['sample_id', 'dataset', 'individual_id', 'chrom', 'pos', 'is_on_target_site']
    col_order = meta_cols + FEATURE_NAMES + ['label']
    out_dir = os.path.dirname(os.path.abspath(args.output_file))
    chunk_path = os.path.join(out_dir, '{}.tmp.parquet'.format(sample_id))
    pd.DataFrame(all_rows).reindex(columns=col_order, fill_value=0).to_parquet(
        chunk_path, index=False)
    print(f"    Chunk written: {chunk_path}", flush=True)
    return chunk_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Build CRISPR training dataset from CRAM files")
    p.add_argument('--mastersheet',      required=True,
                   help='Consolidated mastersheet CSV (assets/all_samples_mastersheet.csv)')
    p.add_argument('--output-file',      required=True,
                   help='Output parquet file (training_data/crispr_training_v2.parquet)')
    p.add_argument('--fasta',            default=DEFAULT_FASTA,
                   help='Reference FASTA file')
    p.add_argument('--window',           type=int,   default=150,
                   help='Window (bp) around each target for feature extraction (default: 150)')
    p.add_argument('--min-indel-fraction', type=float, default=0.05,
                   help='Min edited indel fraction for positive label (default: 0.05)')
    p.add_argument('--max-control-fraction', type=float, default=0.01,
                   help='Max control indel fraction for positive label (default: 0.01)')
    p.add_argument('--n-background-sites', type=int, default=5000,
                   help='Background reads to sample per sample from control BAM (default: 5000)')
    p.add_argument('--threads',          type=int,   default=1,
                   help='Parallel sample processing threads (default: 1)')
    p.add_argument('--samples',          nargs='+',  default=None,
                   help='Only process these sample IDs (optional subset)')
    p.add_argument('--cross-sample-negatives', action='store_true',
                   help='Enable cross-sample hard negative generation: extract reads from '
                        "other samples' target sites that are unedited in the current sample")
    p.add_argument('--cross-sample-sites', type=int, default=500,
                   help='Max target sites per sample to use as cross-sample negatives '
                        '(default: 500)')
    p.add_argument('--offtarget-csv',    default=None,
                   help='CSV of validated off-target hotspots with columns: '
                        'sample_id, chrom, pos[, pam].  Owning sample treats them as '
                        'additional labeled sites; all other samples use them as '
                        'cross-sample negatives.')
    return p.parse_args()




def process_sample_with_retry(row, args, global_targets=None, max_retries=2):
    last_exc = None
    for attempt in range(max_retries + 1):
        try:
            return process_sample(row, args, global_targets)
        except OSError as e:
            msg = str(e).lower()
            if attempt < max_retries and ('truncated' in msg or 'crc' in msg):
                wait = 5 * (attempt + 1)
                print(f"    I/O error on {row['id']} (attempt {attempt+1}), "
                      f"retrying in {wait}s: {e}", flush=True)
                time.sleep(wait)
                last_exc = e
            else:
                raise
    raise last_exc


def main():
    args = parse_args()

    print(f"Loading mastersheet: {args.mastersheet}", flush=True)
    mastersheet = pd.read_csv(args.mastersheet)
    print(f"  {len(mastersheet)} samples", flush=True)

    if args.samples:
        mastersheet = mastersheet[mastersheet['id'].isin(args.samples)]
        print(f"  Filtered to {len(mastersheet)} samples: {args.samples}", flush=True)

    if not os.path.exists(args.fasta):
        print(f"ERROR: Reference FASTA not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)

    out_dir = os.path.dirname(os.path.abspath(args.output_file))
    os.makedirs(out_dir, exist_ok=True)

    # Build global target map before workers if cross-sample negatives are requested.
    # This loads every target file once up-front and is then passed (pickled) to each
    # worker as a read-only lookup table.
    global_targets = None
    if args.cross_sample_negatives or args.offtarget_csv:
        print("Building global target map ...", flush=True)
        global_targets = load_all_targets(mastersheet, args.offtarget_csv)
        print(f"  {len(global_targets)} unique target sites across all samples", flush=True)

    # Workers write their own chunk parquet files and return the path.
    # The main process only collects file paths — no large row data in RAM.
    chunk_files = []
    failed_samples = []

    if args.threads > 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as ex:
            futs = {ex.submit(process_sample_with_retry, row, args, global_targets): row['id']
                    for _, row in mastersheet.iterrows()}
            for fut in concurrent.futures.as_completed(futs):
                sid = futs[fut]
                try:
                    path = fut.result()
                    if path:
                        chunk_files.append(path)
                except Exception as e:
                    print(f"  ERROR processing {sid}: {e}", flush=True)
                    traceback.print_exc()
                    failed_samples.append(sid)
    else:
        for _, row in mastersheet.iterrows():
            try:
                path = process_sample_with_retry(row, args, global_targets)
                if path:
                    chunk_files.append(path)
            except Exception as e:
                print(f"  ERROR processing {row['id']}: {e}", flush=True)
                traceback.print_exc()
                failed_samples.append(row['id'])

    if failed_samples:
        print(f"\nWARNING: {len(failed_samples)} sample(s) failed:", flush=True)
        for s in failed_samples:
            print(f"  - {s}", flush=True)
        failed_file = args.output_file.replace('.parquet', '.failed_samples.txt')
        with open(failed_file, 'w') as fh:
            fh.write('\n'.join(failed_samples) + '\n')
        print(f"Failed IDs written to: {failed_file}", flush=True)

    if not chunk_files:
        print("ERROR: No rows extracted. Check paths and target files.", file=sys.stderr)
        sys.exit(1)

    print(f"\nConcatenating {len(chunk_files)} chunk files ...", flush=True)
    df = pd.concat([pd.read_parquet(f) for f in chunk_files], ignore_index=True)

    # Clean up temp files
    for f in chunk_files:
        try:
            os.remove(f)
        except OSError:
            pass

    n_pos = (df['label'] == 1).sum()
    n_neg = (df['label'] == 0).sum()
    print(f"\nTotal rows: {len(df):,} ({n_pos:,} positive, {n_neg:,} negative)", flush=True)
    print(f"Class ratio (pos/neg): {n_pos/max(n_neg,1):.3f}", flush=True)
    print(f"Datasets: {df['dataset'].value_counts().to_dict()}", flush=True)
    print(f"Individuals: {df['individual_id'].nunique()} unique", flush=True)

    print(f"\nSaving to {args.output_file} ...", flush=True)
    df.to_parquet(args.output_file, index=False)
    print("Done.", flush=True)

    if failed_samples:
        sys.exit(2)  # partial success — caller can detect partial run


if __name__ == '__main__':
    main()
