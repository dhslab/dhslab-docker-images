#!/usr/bin/env python3
"""
build_training_from_tsv.py

Build AAVS1 CRISPR ML training data from pre-curated TP/TN TSV files.

David pre-curated:
  TP TSV: ~77K reads at 2 on-target AAVS1 positions (chr19:55115732,
          chr19:55115752) from edited CRAMs — confirmed CRISPR edits.
  TN TSV: ~72K reads at 6,038 positions genome-wide from the unedited
          control CRAM — reads that look like edits (indels/SA tags) at
          predicted off-target sites but are confirmed non-edits.

This script:
  1. Parses TP and TN TSVs
  2. Builds a temporary sorted+indexed BAM from the embedded SAM records
  3. Extracts the same 16 features as build_crispr_training_dataset.py
  4. Optionally adds cross-site hard negatives from site14 CRAMs
  5. Outputs a parquet with the same schema as the CART training data

Usage:
  python bin/build_training_from_tsv.py \\
    --tp-tsv AAVS1_training_tp.tsv \\
    --tn-tsv AAVS1_training_tn.tsv \\
    --control-cram GM24385-unedited.cram \\
    --target-vcf AAVS1_site5.targets.vcf \\
    --fasta hg38.fa \\
    --output-file training_data/aavs1_training.parquet
"""

import argparse
import importlib.util
import os
import sys
import random
import traceback

import pandas as pd
import pysam


# ---------------------------------------------------------------------------
# Import functions from build_crispr_training_dataset.py
# ---------------------------------------------------------------------------

def _load_build_module():
    """Load build_crispr_training_dataset.py from the same bin/ directory."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    build_path = os.path.join(script_dir, "build_crispr_training_dataset.py")
    if not os.path.exists(build_path):
        raise FileNotFoundError(
            f"Cannot find build_crispr_training_dataset.py at {build_path}"
        )
    spec = importlib.util.spec_from_file_location("build_crispr_training_dataset", build_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_BUILD = _load_build_module()

extract_read_features         = _BUILD.extract_read_features
extract_features_for_site     = _BUILD.extract_features_for_site
_fetch_ctrl_pool              = _BUILD._fetch_ctrl_pool
sample_cross_sample_negatives = _BUILD.sample_cross_sample_negatives
compute_site_indel_fraction   = _BUILD.compute_site_indel_fraction
load_targets                  = _BUILD.load_targets
FEATURE_NAMES                 = _BUILD.FEATURE_NAMES


# ---------------------------------------------------------------------------
# TSV parsing
# ---------------------------------------------------------------------------

def parse_tsv(path):
    """
    Parse a TP or TN TSV file.

    Returns a list of (label_int, chrom, vcf_pos_1based, sam_fields) tuples.
      label_int:      1 for TP, 0 for TN
      chrom:          VCF chromosome string (e.g. 'chr19')
      vcf_pos_1based: int, 1-based VCF position of the target site
      sam_fields:     list[str] — columns[3:] onward (QNAME through optional tags)

    TSV columns (tab-separated):
      Label | VCF_Chrom | VCF_Pos | QNAME | FLAG | RNAME | POS | MAPQ |
      CIGAR | RNEXT | PNEXT | TLEN | SEQ | QUAL | [optional tags...]
    """
    label_map = {'TP': 1, 'TN': 0, '1': 1, '0': 0}
    records = []

    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            cols = line.split('\t')
            # Minimum: 3 TSV meta cols + 11 mandatory SAM fields = 14
            if len(cols) < 14:
                continue

            label_str = cols[0].strip()
            if label_str not in label_map:
                continue  # header row or unrecognised label

            try:
                vcf_pos = int(cols[2])
            except ValueError:
                continue  # skip header-like rows where VCF_Pos is non-numeric

            records.append((label_map[label_str], cols[1], vcf_pos, cols[3:]))

    return records


# ---------------------------------------------------------------------------
# Build temporary sorted+indexed BAM
# ---------------------------------------------------------------------------

def build_temp_bam(all_records, control_cram_path, fasta_path, out_bam_path):
    """
    Build a coordinate-sorted, indexed BAM from SAM records embedded in TSVs.

    Steps:
      1. Open control CRAM to copy its SAM header (full hg38 @SQ lines).
      2. Write a temp SAM: header + one line per record (cols[3:] joined by tab).
      3. pysam.sort → out_bam_path.
      4. pysam.index → out_bam_path.bai.
      5. Remove temp SAM.

    Parameters
    ----------
    all_records:        list of (label_int, chrom, vcf_pos_1based, sam_fields)
    control_cram_path:  path to unedited control CRAM (header source)
    fasta_path:         reference FASTA for CRAM decoding
    out_bam_path:       destination sorted BAM path
    """
    ctrl = pysam.AlignmentFile(control_cram_path, "rc", reference_filename=fasta_path)
    header_text = str(ctrl.header)
    ctrl.close()

    tmp_sam = out_bam_path + ".tmp.sam"
    with open(tmp_sam, 'w') as fh:
        fh.write(header_text)
        for _, _, _, sam_fields in all_records:
            fh.write('\t'.join(sam_fields) + '\n')

    pysam.sort("-o", out_bam_path, tmp_sam)
    pysam.index(out_bam_path)

    try:
        os.remove(tmp_sam)
    except OSError:
        pass

    print(f"  Temp BAM written and indexed: {out_bam_path}", flush=True)


# ---------------------------------------------------------------------------
# Main extraction logic
# ---------------------------------------------------------------------------

def extract_all_features(args):
    """
    Parse TSVs, build temp BAM, extract per-read features, return row dicts.

    Returns a list of feature-dict rows ready for pd.DataFrame conversion.
    """
    print(f"Parsing TP TSV: {args.tp_tsv}", flush=True)
    tp_records = parse_tsv(args.tp_tsv)
    print(f"  {len(tp_records):,} TP records", flush=True)

    print(f"Parsing TN TSV: {args.tn_tsv}", flush=True)
    tn_records = parse_tsv(args.tn_tsv)
    print(f"  {len(tn_records):,} TN records", flush=True)

    # Build temp BAM in output directory
    out_dir = os.path.dirname(os.path.abspath(args.output_file))
    os.makedirs(out_dir, exist_ok=True)
    tmp_bam_path = os.path.join(out_dir, f"{args.sample_id}.tmp.bam")

    print("Building temporary sorted BAM from TSV records ...", flush=True)
    build_temp_bam(tp_records + tn_records, args.control_cram, args.fasta, tmp_bam_path)

    # Load target VCF for PAM positions and on-target flags
    print(f"Loading target VCF: {args.target_vcf}", flush=True)
    targets_df = load_targets(args.target_vcf)
    print(f"  {len(targets_df):,} target sites", flush=True)

    # on_target_set: (chrom, pos_0based) for sites with TARGET flag
    on_target_set = set()
    for _, row in targets_df[targets_df['Ontarget'] == 1].iterrows():
        on_target_set.add((row['Chromosome'], int(row['Start'])))
    print(f"  {len(on_target_set)} on-target sites (TARGET flag)", flush=True)

    # pam_by_chrom: chrom -> list of 1-based PAM positions (for distance feature)
    pam_by_chrom = {}
    for _, row in targets_df.iterrows():
        pam_by_chrom.setdefault(row['Chromosome'], []).append(int(row['Pos']))

    # Collect unique sites preserving insertion order (dict keys)
    tp_sites = {}
    for _, chrom, vcf_pos_1based, _ in tp_records:
        tp_sites[(chrom, vcf_pos_1based - 1)] = True

    tn_sites = {}
    for _, chrom, vcf_pos_1based, _ in tn_records:
        tn_sites[(chrom, vcf_pos_1based - 1)] = True

    print(f"  Unique TP sites: {len(tp_sites)}", flush=True)
    print(f"  Unique TN sites: {len(tn_sites)}", flush=True)

    # Open BAMs for feature extraction
    temp_bam = pysam.AlignmentFile(tmp_bam_path, "rb")
    ctrl_bam = pysam.AlignmentFile(args.control_cram, "rc", reference_filename=args.fasta)
    fasta_file = pysam.FastaFile(args.fasta)

    all_rows = []

    # --- TP sites (label=1) ---
    print(f"\nExtracting features for {len(tp_sites)} TP sites ...", flush=True)
    n_tp_sites = len(tp_sites)
    for i, (chrom, pos_0based) in enumerate(tp_sites, 1):
        is_on_target = 1 if (chrom, pos_0based) in on_target_set else 0
        pam_positions = pam_by_chrom.get(chrom, [])
        try:
            site_rows = extract_features_for_site(
                temp_bam, ctrl_bam, chrom, pos_0based,
                pam_positions, label=1, is_on_target=is_on_target,
                window=args.window,
                fasta_file=fasta_file,
            )
            all_rows.extend(site_rows)
        except Exception as e:
            print(f"  WARNING: TP site ({chrom},{pos_0based}) failed: {e}", flush=True)
        if i % 10 == 0 or i == n_tp_sites:
            print(f"  TP {i}/{n_tp_sites} sites done — {len(all_rows):,} rows", flush=True)

    tp_row_count = len(all_rows)
    print(f"  TP rows total: {tp_row_count:,}", flush=True)

    # --- TN sites (label=0, is_on_target=0) ---
    print(f"\nExtracting features for {len(tn_sites)} TN sites ...", flush=True)
    n_tn_sites = len(tn_sites)
    tn_start = len(all_rows)
    for i, (chrom, pos_0based) in enumerate(tn_sites, 1):
        pam_positions = pam_by_chrom.get(chrom, [])
        try:
            site_rows = extract_features_for_site(
                temp_bam, ctrl_bam, chrom, pos_0based,
                pam_positions, label=0, is_on_target=0,
                window=args.window,
                fasta_file=fasta_file,
            )
            all_rows.extend(site_rows)
        except Exception as e:
            print(f"  WARNING: TN site ({chrom},{pos_0based}) failed: {e}", flush=True)
        if i % 500 == 0 or i == n_tn_sites:
            tn_rows_so_far = len(all_rows) - tn_start
            print(f"  TN {i}/{n_tn_sites} sites done — {tn_rows_so_far:,} TN rows",
                  flush=True)

    tn_row_count = len(all_rows) - tp_row_count
    print(f"  TN rows total: {tn_row_count:,}", flush=True)

    temp_bam.close()
    ctrl_bam.close()
    fasta_file.close()

    # Clean up temp BAM (index file is .bam.bai)
    for suffix in ('', '.bai'):
        try:
            os.remove(tmp_bam_path + suffix)
        except OSError:
            pass

    # Tag main-dataset rows with sample metadata
    for r in all_rows:
        r['sample_id']     = args.sample_id
        r['dataset']       = args.dataset
        r['individual_id'] = args.individual_id

    # --- Cross-site hard negatives (optional) ---
    # site5 TP positions should be unedited in site14 CRAMs → hard negatives:
    # the model sees is_at_any_target_site=1 with no real edit signal.
    if args.cross_site_crams:
        # Use all predicted target sites (from target VCF) as cross-site candidates.
        # This gives ~5k–7k genome-wide positions; only the handful that overlap the
        # other guide's on-target cut will be filtered by indel_frac, leaving genuine
        # hard negatives.  Using the actual VCF Pos column avoids the p+1 approximation.
        cross_sites = [
            (row['Chromosome'], int(row['Start']), int(row['Pos']))
            for _, row in targets_df.iterrows()
        ]

        print(f"\nGenerating cross-site hard negatives from "
              f"{len(args.cross_site_crams)} CRAM pair(s) ...", flush=True)

        for pair_str in args.cross_site_crams:
            pair_str = pair_str.strip()
            parts = pair_str.split(':', 1)
            if len(parts) != 2:
                print(f"  WARNING: skipping malformed --cross-site-crams entry: "
                      f"'{pair_str}' (expected 'edited.cram:ctrl.cram')", flush=True)
                continue

            edited_cram_path = parts[0].strip()
            cross_ctrl_path  = parts[1].strip()
            sample_stem = os.path.splitext(os.path.basename(edited_cram_path))[0]
            print(f"  Pair: {sample_stem}", flush=True)

            try:
                cross_edited_bam = pysam.AlignmentFile(
                    edited_cram_path, "rc", reference_filename=args.fasta
                )
                cross_ctrl_bam = pysam.AlignmentFile(
                    cross_ctrl_path, "rc", reference_filename=args.fasta
                )
            except Exception as e:
                print(f"  WARNING: could not open CRAMs for {sample_stem}: {e}", flush=True)
                continue

            try:
                pair_rows = sample_cross_sample_negatives(
                    edited_bam=cross_edited_bam,
                    ctrl_bam=cross_ctrl_bam,
                    cross_sample_sites=cross_sites,
                    max_sites=args.cross_site_sites,
                    min_indel_frac=0.01,
                    window=args.window,
                    rng=random.Random(42),
                )
            except Exception as e:
                print(f"  WARNING: cross-site extraction failed for {sample_stem}: {e}",
                      flush=True)
                traceback.print_exc()
                pair_rows = []
            finally:
                cross_edited_bam.close()
                cross_ctrl_bam.close()

            for r in pair_rows:
                r['sample_id']     = sample_stem
                r['dataset']       = args.dataset
                r['individual_id'] = args.individual_id

            all_rows.extend(pair_rows)
            print(f"    {len(pair_rows):,} cross-site negative rows from {sample_stem}",
                  flush=True)

    return all_rows


# ---------------------------------------------------------------------------
# argparse
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Build AAVS1 CRISPR training data from pre-curated TP/TN TSV files"
    )
    p.add_argument('--tp-tsv',        required=True,
                   help='TSV of TP (confirmed-edit) SAM records with site labels')
    p.add_argument('--tn-tsv',        required=True,
                   help='TSV of TN (confirmed non-edit) SAM records with site labels')
    p.add_argument('--control-cram',  required=True,
                   help='Unedited control CRAM — supplies SAM header and control features')
    p.add_argument('--target-vcf',    required=True,
                   help='Target VCF with PAM positions and TARGET info flag '
                        '(e.g. AAVS1_site5.targets.vcf)')
    p.add_argument('--fasta',         required=True,
                   help='Reference FASTA for CRAM/BAM decoding')
    p.add_argument('--output-file',   required=True,
                   help='Output parquet path (e.g. training_data/aavs1_training.parquet)')
    p.add_argument('--sample-id',     default='AAVS1_site5_curated',
                   help='sample_id tag for main TSV-derived rows '
                        '(default: AAVS1_site5_curated)')
    p.add_argument('--dataset',       default='AAVS1_curated',
                   help='dataset tag for all output rows (default: AAVS1_curated)')
    p.add_argument('--individual-id', default='GM24385',
                   help='individual_id tag (default: GM24385)')
    p.add_argument('--window',        type=int, default=150,
                   help='Window (bp) around each target site for feature extraction '
                        '(default: 150)')
    p.add_argument('--cross-site-crams', nargs='*', default=None,
                   metavar='EDITED_CRAM:CTRL_CRAM',
                   help='Space-separated "edited.cram:ctrl.cram" pairs for cross-site '
                        'hard negative generation (site5 TP positions confirmed unedited '
                        'in site14 CRAMs)')
    p.add_argument('--cross-site-sites', type=int, default=500,
                   help='Max TP sites to use as cross-site negatives per CRAM pair '
                        '(default: 500)')
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    # Validate required input paths
    for flag, path in [
        ('--tp-tsv',       args.tp_tsv),
        ('--tn-tsv',       args.tn_tsv),
        ('--control-cram', args.control_cram),
        ('--target-vcf',   args.target_vcf),
        ('--fasta',        args.fasta),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {flag} not found: {path}", file=sys.stderr)
            sys.exit(1)

    print(f"TP TSV:       {args.tp_tsv}")
    print(f"TN TSV:       {args.tn_tsv}")
    print(f"Control CRAM: {args.control_cram}")
    print(f"Target VCF:   {args.target_vcf}")
    print(f"Output:       {args.output_file}")
    print(f"sample_id:    {args.sample_id}")
    print(f"dataset:      {args.dataset}")
    print(f"individual_id:{args.individual_id}")
    print(f"window:       {args.window}", flush=True)

    all_rows = extract_all_features(args)

    if not all_rows:
        print("ERROR: no rows extracted. Check input files and paths.", file=sys.stderr)
        sys.exit(1)

    # Build DataFrame with same schema as build_crispr_training_dataset.py
    # is_on_target_site is metadata (not in FEATURE_NAMES); kept for post-hoc analysis.
    meta_cols = ['sample_id', 'dataset', 'individual_id', 'chrom', 'pos', 'is_on_target_site']
    col_order = meta_cols + FEATURE_NAMES + ['label']

    df = pd.DataFrame(all_rows).reindex(columns=col_order, fill_value=0)

    n_pos = int((df['label'] == 1).sum())
    n_neg = int((df['label'] == 0).sum())
    n_ontarget_pos = int(
        ((df['label'] == 1) & (df['is_on_target_site'] == 1)).sum()
    )

    print(f"\nTotal rows:        {len(df):,}")
    print(f"  label=1 (TP):    {n_pos:,}")
    print(f"  label=0 (TN):    {n_neg:,}")
    print(f"  Class ratio (pos/neg): {n_pos / max(n_neg, 1):.4f}")
    print(f"  On-target positives (is_on_target_site=1, label=1): {n_ontarget_pos:,}")
    print(f"  Unique sample_ids: {df['sample_id'].nunique()}", flush=True)

    out_dir = os.path.dirname(os.path.abspath(args.output_file))
    os.makedirs(out_dir, exist_ok=True)

    print(f"\nSaving to {args.output_file} ...", flush=True)
    df.to_parquet(args.output_file, index=False, engine='pyarrow')
    print("Done.", flush=True)


if __name__ == '__main__':
    main()
