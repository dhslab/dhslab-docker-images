#!/usr/bin/env python3
import sys
import argparse
import pysam
import re
import edlib

def reverse_complement(seq):
    """Returns the reverse complement of a DNA sequence."""
    if not seq: return ""
    complement = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(complement)[::-1]

def get_bnd_parts(alt_str):
    """
    Parses VCF BND strings (e.g., ]chr1:123]T) into components.
    Returns: (pre_bases, bracket_char, remote_chrom, remote_pos, post_bases)
    """
    # Regex captures: 1=Pre, 2=Bracket, 3=RemoteLoc, 4=Bracket, 5=Post
    pattern = re.compile(r"([ACGTNacgtn]*)(\[|\])([^:]+:\d+)(\[|\])([ACGTNacgtn]*)")
    match = pattern.fullmatch(alt_str)
    if not match:
        return None
    return match.groups()

def generate_contig(record, ref_file, flank_length):
    """
    Constructs the theoretical alternate sequence (Contig) using the Reference genome.
    Logic:
      1. Identify SV type (Linear vs BND).
      2. Fetch Local Flank (Left or Right depending on break orientation).
      3. Fetch Remote Flank (and RC if necessary based on bracket direction).
      4. Stitch them together.
    """
    chrom = record.chrom
    pos = record.pos # 1-based
    ref_base = record.ref
    alt_base = record.alts[0]
    svtype = record.info.get('SVTYPE', 'Unknown')
    strand = record.info.get('STRAND', '')
    
    # Pysam 0-based conversion
    start_idx = pos - 1 
    end_idx = start_idx + len(ref_base)

    try:
        # --- Handle Breakends (BND) ---
        if svtype == 'BND' or '[' in alt_base or ']' in alt_base:
            parts = get_bnd_parts(alt_base)
            if not parts:
                return "Error: Unparseable BND string"

            pre_bases, bracket, remote_loc, _, post_bases = parts
            r_chrom, r_pos = remote_loc.split(':')
            r_pos = int(r_pos)
            
            # 1. Get Local Sequence
            if pre_bases: 
                # Local is upstream (Left) -> Break
                local_seq = ref_file.fetch(chrom, max(0, start_idx - flank_length), end_idx - 1)

                if bracket == '[':
                    r_start_0 = r_pos - 1
                    remote_seq = ref_file.fetch(r_chrom, r_start_0, r_start_0 + flank_length)
                else: # bracket == ']'
                    # Remote points Reverse (Left/Upstream of r_pos) -> RC needed
                    r_end_0 = r_pos
                    raw_remote = ref_file.fetch(r_chrom, max(0, r_end_0 - flank_length), r_end_0)
                    remote_seq = reverse_complement(raw_remote)
                    
                return local_seq + pre_bases + remote_seq if strand=="+" else reverse_complement(local_seq + pre_bases + remote_seq)

            else: # joined before t (ref is downstream)
                if bracket == '[':
                    local_seq = ref_file.fetch(chrom, start_idx, end_idx + flank_length)
                    # Remote points Forward (Right/Downstream of r_pos)
                    r_start_0 = r_pos - 1
                    raw_remote = ref_file.fetch(r_chrom, r_start_0 + 1, r_start_0 + flank_length)
                    remote_seq = reverse_complement(raw_remote)

                else: # bracket == ']'
                    local_seq = ref_file.fetch(chrom, start_idx + 1, end_idx + flank_length)
                    # Remote points Reverse (Left/Upstream of r_pos) -> RC needed
                    r_end_0 = r_pos
                    remote_seq = ref_file.fetch(r_chrom, max(0, r_end_0 - flank_length), r_end_0)

                return remote_seq + post_bases + local_seq if strand=="+" else reverse_complement(remote_seq + post_bases + local_seq)

        # --- Handle Linear SVs (DEL, INS, etc) ---
        else:
            if len(alt_base) > len(ref_base) or record.info.get('SVLEN', 0) > 0: # DUP            
                left_flank = ref_file.fetch(chrom, max(0, start_idx - flank_length), start_idx)
                right_flank = ref_file.fetch(chrom, end_idx, end_idx + flank_length)
                # For linear SVs, we sandwich the ALT string between the flanks
            
            elif len(alt_base) < len(ref_base) or record.info.get('SVLEN', 0) < 0: #DEL
                left_flank = ref_file.fetch(chrom, max(0, start_idx - flank_length), start_idx)
                right_flank = ref_file.fetch(chrom, end_idx, end_idx + flank_length)
                # For linear SVs, we sandwich the ALT string between the flanks

            seq = f"{left_flank}{alt_base}{right_flank}" if strand=="+" else reverse_complement(f"{left_flank}{alt_base}{right_flank}")
            
            return seq 

    except KeyError as e:
        return f"Error: Chromosome/Region not found in reference ({e})"
    except Exception as e:
        return f"Error: {str(e)}"

def print_alignment(query_seq, target_seq):
    """
    Aligns Query (Generated) to Target (VCF SEQ) using Edlib and prints visualization.
    Mode: NW (Global) - forces end-to-end alignment to show flank discrepancies.
    """
    # Use Needleman-Wunsch (global) to force end-to-end comparison
    # This highlights if the flanks in the VCF SEQ differ from the Reference flanks
    result = edlib.align(query_seq, target_seq, mode="NW", task="path")
    
    if result['editDistance'] == 0:
        return True, 0

    # Get visual alignment
    nice = edlib.getNiceAlignment(result, query_seq, target_seq)
    
    print("\n    Alignment (Top: Generated, Bottom: VCF SEQ):")
    print(f"    Q: {nice['query_aligned']}")
    print(f"       {nice['matched_aligned']}")
    print(f"    T: {nice['target_aligned']}")
    
    return False, result['editDistance']

def main():
    parser = argparse.ArgumentParser(description="Reconstruct SV contigs and align to VCF SEQ using Edlib.")
    parser.add_argument("vcf")
    parser.add_argument("-r", "--reference", required=True, help="Path to indexed reference FASTA (.fa)")
    parser.add_argument("-f", "--flank", type=int, default=200, help="Flank length for reconstruction (default: 50)")
    args = parser.parse_args()

    # Open Reference
    try:
        ref_file = pysam.FastaFile(args.reference)
    except ValueError:
        sys.stderr.write(f"CRITICAL ERROR: Could not open/index reference file: {args.reference}\n")
        sys.exit(1)

    # Open VCF from Stdin
    try:
        vcf_in = pysam.VariantFile(args.vcf)
    except Exception as e:
        sys.stderr.write(f"CRITICAL ERROR: Could not read VCF from stdin. {e}\n")
        sys.exit(1)

    # Header
    print(f"{'CHROM':<10} {'POS':<10} {'SVTYPE':<10} {'EDITS':<6} {'STATUS':<10}")
    print("-" * 80)

    for record in vcf_in:
        vcf_seq = record.info.get('SEQ', None)
        if record.info.get("STRAND","+") == "-":
            vcf_seq = reverse_complement(vcf_seq) 

        svtype = record.info.get('SVTYPE', 'UNK')
        
        # Generate the theoretical sequence from Reference
        generated_seq = generate_contig(record, ref_file, args.flank)

        # Handle Errors in generation
        if "Error" in generated_seq:
             print(f"{record.chrom:<10} {str(record.pos):<10} {svtype:<10} {'-':<6} {'ERROR':<10}")
             print(f"  > {generated_seq}\n")
             continue

        if not vcf_seq:
            print(f"{record.chrom:<10} {str(record.pos):<10} {svtype:<10} {'-':<6} {'NO_SEQ':<10}")
            continue

        # Perform Alignment
        # We print the row first to keep logs organized, then print alignment details if mismatched
        
        # Run alignment silently first to get score
        result = edlib.align(vcf_seq, generated_seq, mode="HW", task="path")
        edit_dist = result['editDistance']
        status = "MATCH" if edit_dist == 0 else "MISMATCH"

        if status == "MISMATCH":
            print(f"{str(record.id)}\t{record.chrom:<10} {str(record.pos):<10} {record.alts[0]} {svtype:<10} {str(edit_dist):<6} {status:<10}")

            # Call helper to print the nice visualization
            align = edlib.getNiceAlignment(result, vcf_seq, generated_seq, gapSymbol='-')
            print(align['query_aligned'])
            print(align['matched_aligned'])
            print(align['target_aligned'])
            print(f"SEQ: {vcf_seq}")
            print(f"contig: {generated_seq}")
            print("") # spacer

if __name__ == "__main__":
    main()