#!/usr/bin/env python3
"""
crispr_ml_features.py

Shared feature computation utilities for CRISPR ML training and inference.

Imported via importlib by both build_crispr_training_dataset.py and
extract_variant_reads_ML.py to ensure training and inference use identical
feature logic.

Sentinel values for missing distances: 999 (not -1).  This puts the
missing-value token far above any realistic distance, allowing tree models
to isolate it with a single high-threshold split rather than confusing it
with short distances near zero.
"""


def _softclip_pam_distance(read, pam_positions):
    """
    Min distance (bp) from any softclip reference junction in this read to
    the nearest PAM position.

    A softclip junction is the reference coordinate immediately adjacent to a
    soft-clipped block: the first (or last) aligned reference base bordering
    the clip.  For a leading clip the junction equals read.reference_start;
    for a trailing clip it equals the running reference cursor after all
    preceding reference-consuming ops.

    Parameters
    ----------
    read : pysam.AlignedSegment
    pam_positions : list[int]
        PAM positions in whatever coordinate convention the caller uses
        (typically 1-based VCF Pos, consistent with distance_to_closest_pam).

    Returns
    -------
    int : min distance in bp, or 999 if no softclips or pam_positions is empty.
    """
    if not pam_positions or read.cigartuples is None:
        return 999

    junctions = []
    ref_pos = read.reference_start  # 0-based

    for op, length in read.cigartuples:
        if op == 4:  # S — soft clip; record the current reference cursor
            junctions.append(ref_pos)
        # Advance reference cursor for ops that consume the reference
        if op in (0, 2, 3, 7, 8):  # M, D, N, =, X
            ref_pos += length

    if not junctions:
        return 999

    return min(abs(j - p) for j in junctions for p in pam_positions)


def _find_tandem_repeats(seq, min_unit=1, max_unit=6, min_copies=4):
    """
    Locate tandem repeat regions in *seq*.

    Scans left-to-right for runs of unit strings with length in
    [min_unit, max_unit] repeated at least min_copies consecutive times.

    Parameters
    ----------
    seq : str
        Reference sequence (uppercase recommended).
    min_unit : int
        Minimum repeat-unit length (bp).
    max_unit : int
        Maximum repeat-unit length (bp).
    min_copies : int
        Minimum number of consecutive copies to qualify as a repeat.

    Returns
    -------
    list of (start, end) tuples — 0-based, half-open, in seq coordinates.
    """
    repeats = []
    n = len(seq)
    for unit_len in range(min_unit, max_unit + 1):
        i = 0
        while i <= n - unit_len * min_copies:
            unit = seq[i:i + unit_len]
            j = i + unit_len
            copies = 1
            while j + unit_len <= n and seq[j:j + unit_len] == unit:
                copies += 1
                j += unit_len
            if copies >= min_copies:
                repeats.append((i, j))
                i = j   # skip past this repeat to avoid double-counting
            else:
                i += 1
    return repeats


def _precompute_site_ref_features(fasta_file, chrom, pos_0based, window=500):
    """
    Compute per-site reference-sequence features.

    Called once per target site; results are reused for all reads at that site.

    Parameters
    ----------
    fasta_file : pysam.FastaFile
        Open FASTA handle; fetch() uses 0-based half-open coordinates.
    chrom : str
    pos_0based : int
        0-based genomic position (pysam convention).
    window : int
        Number of bp to scan on each side when looking for microsatellites.

    Returns
    -------
    dict with keys:
        'is_in_homopolymer'     : int, 1 if pos is inside a run of ≥5 identical bases
        'dist_to_microsatellite': int, distance (bp) to nearest tandem repeat;
                                  0 = pos is inside a repeat; 999 = none found
    """
    start = max(0, pos_0based - window)
    end   = pos_0based + window + 1
    try:
        ref_seq = fasta_file.fetch(chrom, start, end).upper()
    except Exception:
        return {'is_in_homopolymer': 0, 'dist_to_microsatellite': 999}

    site_idx = pos_0based - start  # index within ref_seq for the target position

    # ── Homopolymer detection ─────────────────────────────────────────────────
    is_in_homopolymer = 0
    if 0 <= site_idx < len(ref_seq):
        base = ref_seq[site_idx]
        if base in 'ACGT':
            run_lo = site_idx
            while run_lo > 0 and ref_seq[run_lo - 1] == base:
                run_lo -= 1
            run_hi = site_idx
            while run_hi + 1 < len(ref_seq) and ref_seq[run_hi + 1] == base:
                run_hi += 1
            if (run_hi - run_lo + 1) >= 5:
                is_in_homopolymer = 1

    # ── Microsatellite / tandem-repeat distance ───────────────────────────────
    repeats    = _find_tandem_repeats(ref_seq)
    dist_to_ms = 999
    for rep_start, rep_end in repeats:
        if rep_start <= site_idx < rep_end:
            dist_to_ms = 0
            break
        d = min(abs(site_idx - rep_start), abs(site_idx - (rep_end - 1)))
        if d < dist_to_ms:
            dist_to_ms = d

    return {
        'is_in_homopolymer':      is_in_homopolymer,
        'dist_to_microsatellite': dist_to_ms,
    }
