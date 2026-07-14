#!/usr/bin/env python3
"""
build_from_curated_tsv.py

Build a CRISPR training parquet from David's pre-curated TP/TN TSV files.

Each TSV row: Label<TAB>VCF_Chrom<TAB>VCF_Pos<TAB><full SAM record, tab-delimited>

Reads are already quality-filtered (MAPQ>=20, <=4 SNP mismatches, no dup/secondary/supp,
guaranteed to contain indels, large soft-clips, or supplementary alignments).

The output parquet has the same schema as training_data/crispr_training_v2.parquet
produced by build_crispr_training_dataset.py:
  sample_id | dataset | individual_id | chrom | pos | <14 features> | label

Usage:
  python build_from_curated_tsv.py \\
    --tp-tsv /storage2/.../AAVS1_training_tp.tsv \\
    --tn-tsv /storage2/.../AAVS1_training_tn.tsv \\
    --control-cram /storage2/.../GM24385-unedited.cram \\
    --fasta /storage2/.../hg38_PLVM_CD19_CARv4_cd34.fa \\
    --output-file training_data/aavs1_tptn_only.parquet
"""

import argparse
import importlib.util
import os
import sys

# Add repo root to sys.path so locally-installed packages (pyarrow, joblib, etc.) are found
_repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

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
        raise FileNotFoundError("Cannot find extract_variant_reads_ML.py at {}".format(ml_path))
    spec = importlib.util.spec_from_file_location("extract_variant_reads_ML", ml_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_ML = _load_ml_module()

parse_cigar_once      = _ML.parse_cigar_once
count_mismatches_fast = _ML.count_mismatches_fast

# ---------------------------------------------------------------------------
# Import shared feature helpers from crispr_ml_features.py
# ---------------------------------------------------------------------------

def _load_features_module():
    """Load crispr_ml_features.py from the same bin/ directory."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    feat_path = os.path.join(script_dir, "crispr_ml_features.py")
    if not os.path.exists(feat_path):
        raise FileNotFoundError(
            "Cannot find crispr_ml_features.py at {}".format(feat_path)
        )
    spec = importlib.util.spec_from_file_location("crispr_ml_features", feat_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_FEAT = _load_features_module()

_softclip_pam_distance        = _FEAT._softclip_pam_distance
_precompute_site_ref_features = _FEAT._precompute_site_ref_features

# ---------------------------------------------------------------------------
# Position-level caches — populated once per unique (chrom, pos)
# ---------------------------------------------------------------------------
# calculate_control_fractions depends only on (chrom, pos), not the read.
# calculate_exclusivity_features does a CRAM fetch per read, but the fetched
# control-read pool is identical for all reads at the same (chrom, pos).
# We cache the pool as lightweight tuples and re-run the exclusivity logic
# inline, eliminating ~N_reads CRAM fetches and replacing them with ~N_positions.
#
# Impact on AAVS1 data:
#   TP: 77k reads at 2 positions  → 2 CRAM fetches (vs 77k)
#   TN: 72k reads at 6037 positions → 6037 CRAM fetches (vs 72k)

_ctrl_fracs_cache = {}  # (chrom, pos) -> fractions dict
_ctrl_reads_cache = {}  # (chrom, pos) -> list of (start, end, has_ins, has_del)
_site_ref_cache   = {}  # (chrom, pos) -> dict from _precompute_site_ref_features()


def _cached_ctrl_fracs(control_bam, chrom, pos, window=50):
    key = (chrom, pos)
    if key not in _ctrl_fracs_cache:
        _ctrl_fracs_cache[key] = _ML.calculate_control_fractions(
            control_bam, chrom, pos, window=window)
    return _ctrl_fracs_cache[key]


def _get_ctrl_reads(control_bam, chrom, pos, window=100):
    """Return cached list of (start, end, has_ins, has_del) for control reads at pos."""
    key = (chrom, pos)
    if key not in _ctrl_reads_cache:
        pool = []
        start = max(0, pos - window)
        for cr in control_bam.fetch(chrom, start, pos + window):
            if cr.is_unmapped or cr.is_duplicate or not cr.cigartuples:
                continue
            pool.append((
                cr.reference_start,
                cr.reference_end,
                any(op == 1 for op, _ in cr.cigartuples),  # has_ins
                any(op == 2 for op, _ in cr.cigartuples),  # has_del
            ))
        _ctrl_reads_cache[key] = pool
    return _ctrl_reads_cache[key]


def _cached_exclusivity(read, control_bam, chrom, pos, window=100):
    """Exclusivity features using the cached control-read pool (no per-read CRAM fetch)."""
    result = {'deletion_exclusive_to_edited': 0}
    if read.is_unmapped or read.is_duplicate or not read.cigartuples:
        return result

    has_del = any(op == 2 for op, _ in read.cigartuples)
    r_start = read.reference_start
    r_end   = read.reference_end

    ctrl_has_del = False
    for ctrl_start, ctrl_end, c_has_ins, c_has_del in _get_ctrl_reads(
            control_bam, chrom, pos, window):
        # overlapping reads
        if r_end >= ctrl_start and ctrl_end >= r_start:
            if has_del and c_has_del:
                ctrl_has_del = True

    result['deletion_exclusive_to_edited'] = 1 if (has_del and not ctrl_has_del) else 0
    return result


def _cached_site_ref_features(fasta_file, chrom, pos):
    """Return cached ref-sequence features for (chrom, pos).

    pos is the VCF/1-based coordinate from the TSV; converted to 0-based for
    _precompute_site_ref_features.
    """
    key = (chrom, pos)
    if key not in _site_ref_cache:
        _site_ref_cache[key] = _precompute_site_ref_features(
            fasta_file, chrom, pos - 1)  # VCF 1-based → pysam 0-based
    return _site_ref_cache[key]


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


# ---------------------------------------------------------------------------
# Per-read feature extraction (mirrors build_crispr_training_dataset.py)
# ---------------------------------------------------------------------------

def extract_read_features(read, control_bam, chrom, pos, pam_positions,
                           fasta_file=None):
    """
    Extract the 14 CRISPR read features for a single pysam.AlignedSegment.
    Returns a dict or None if the read is unusable.
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

    ctrl_fracs = _cached_ctrl_fracs(control_bam, chrom, pos, window=50)
    excl       = _cached_exclusivity(read, control_bam, chrom, pos, window=100)

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

    if fasta_file is not None:
        ref_feats = _cached_site_ref_features(fasta_file, chrom, pos)
    else:
        ref_feats = {}

    return {
        'read_pair_gap': abs(read.template_length) if (
            hasattr(read, 'template_length') and read.is_paired and read.is_proper_pair
        ) else -1,
        'read_softclip':  cigar_data['softclips'],
        'read_del_vs_control':      has_del - ctrl_fracs['fraction_control_reads_del'],
        'read_ins_vs_control':      has_ins - ctrl_fracs['fraction_control_reads_ins'],
        'read_mismatch_vs_control': has_mm  - ctrl_fracs['fraction_control_reads_mismatch'],
        'deletion_exclusive_to_edited': excl['deletion_exclusive_to_edited'],
        'is_at_any_target_site':   is_at_any,
        'distance_to_closest_pam': dist_to_pam,
        'total_indel_size':        total_indel,
        'indel_size_category':     indel_cat,
        'indel_complexity_score':  complexity,
        'softclip_dist_to_PAM':    _softclip_pam_distance(read, pam_positions),
        'is_in_homopolymer':       ref_feats.get('is_in_homopolymer', 0),
        'dist_to_microsatellite':  ref_feats.get('dist_to_microsatellite', 999),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Build CRISPR training parquet from pre-curated TP/TN TSV files"
    )
    p.add_argument('--tp-tsv',        required=True,
                   help='Path to curated TP TSV (Label|VCF_Chrom|VCF_Pos|SAM fields...)')
    p.add_argument('--tn-tsv',        required=True,
                   help='Path to curated TN TSV (Label|VCF_Chrom|VCF_Pos|SAM fields...)')
    p.add_argument('--control-cram',  required=True,
                   help='Unedited control CRAM for feature extraction')
    p.add_argument('--fasta',         default=DEFAULT_FASTA,
                   help='Reference FASTA (default: pipeline reference)')
    p.add_argument('--output-file',   required=True,
                   help='Output parquet path (e.g. training_data/aavs1_tptn_only.parquet)')
    p.add_argument('--dataset',       default='AAVS1_curated',
                   help='Dataset label for metadata column (default: AAVS1_curated)')
    p.add_argument('--individual-id', default='GM24385',
                   help='Individual ID for CV grouping (default: GM24385)')
    p.add_argument('--sample-id',     default='AAVS1_curated',
                   help='Sample ID for metadata column (default: AAVS1_curated)')
    return p.parse_args()


def main():
    args = parse_args()

    # Validate inputs
    for label, path in [('--tp-tsv', args.tp_tsv),
                         ('--tn-tsv', args.tn_tsv),
                         ('--control-cram', args.control_cram),
                         ('--fasta', args.fasta)]:
        if not os.path.exists(path):
            print("ERROR: {} not found: {}".format(label, path), file=sys.stderr)
            sys.exit(1)

    os.makedirs(os.path.dirname(os.path.abspath(args.output_file)), exist_ok=True)

    print("Opening control CRAM: {}".format(args.control_cram), flush=True)
    control_bam = pysam.AlignmentFile(
        args.control_cram, "rc", reference_filename=args.fasta
    )
    header = control_bam.header

    print("Opening reference FASTA: {}".format(args.fasta), flush=True)
    fasta = pysam.FastaFile(args.fasta)

    # TSV sources: (path, label_int)
    sources = [
        (args.tp_tsv, 1),
        (args.tn_tsv, 0),
    ]

    all_rows = []
    skipped  = 0

    for tsv_path, label_int in sources:
        source_name = "TP" if label_int == 1 else "TN"
        print("Processing {} TSV: {}".format(source_name, tsv_path), flush=True)

        processed = 0
        kept = 0

        with open(tsv_path) as fh:
            for line in fh:
                line = line.rstrip('\n')
                if not line or line.startswith('Label'):
                    continue

                parts = line.split('\t')
                if len(parts) < 4:
                    skipped += 1
                    continue

                chrom   = parts[1]
                try:
                    pos = int(parts[2])
                except ValueError:
                    skipped += 1
                    continue

                sam_line = '\t'.join(parts[3:])

                try:
                    read = pysam.AlignedSegment.fromstring(sam_line, header)
                except Exception as e:
                    skipped += 1
                    if skipped <= 5:
                        print("  WARNING: could not parse SAM record ({}): {}".format(
                            e, sam_line[:80]), flush=True)
                    continue

                processed += 1

                # VCF_Pos is the PAM position; use directly
                feat = extract_read_features(
                    read, control_bam, chrom, pos,
                    pam_positions=[pos],
                    fasta_file=fasta,
                )
                if feat is None:
                    skipped += 1
                    continue

                feat['label']         = label_int
                feat['chrom']         = chrom
                feat['pos']           = pos
                feat['sample_id']     = args.sample_id
                feat['dataset']       = parts[0]        # "TP" or "TN" from TSV Label column
                feat['individual_id'] = args.individual_id

                all_rows.append(feat)
                kept += 1

                if processed % 10000 == 0:
                    print("  {} reads processed, {} kept so far ...".format(
                        processed, kept), flush=True)

        print("  {} reads processed, {} kept, {} skipped".format(
            processed, kept, skipped), flush=True)

    fasta.close()
    control_bam.close()

    if not all_rows:
        print("ERROR: No rows extracted. Check TSV paths and SAM record format.",
              file=sys.stderr)
        sys.exit(1)

    # Build DataFrame with the canonical column order
    meta_cols = ['sample_id', 'dataset', 'individual_id', 'chrom', 'pos']
    col_order = meta_cols + FEATURE_NAMES + ['label']

    df = pd.DataFrame(all_rows).reindex(columns=col_order, fill_value=0)

    n_pos = int((df['label'] == 1).sum())
    n_neg = int((df['label'] == 0).sum())
    print("\nTotal rows: {:,} ({:,} positive / TP, {:,} negative / TN)".format(
        len(df), n_pos, n_neg), flush=True)
    print("Class ratio (pos/neg): {:.3f}".format(n_pos / max(n_neg, 1)), flush=True)
    print("Columns: {}".format(df.columns.tolist()), flush=True)
    print("Dataset values: {}".format(df['dataset'].value_counts().to_dict()), flush=True)
    print("Individual IDs: {}".format(df['individual_id'].unique().tolist()), flush=True)

    print("\nSaving to {} ...".format(args.output_file), flush=True)
    df.to_parquet(args.output_file, index=False)
    print("Done.", flush=True)


if __name__ == '__main__':
    main()
