import argparse
import pysam
import numpy as np
import gzip
import random
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Simulate reads for a region defined by a center point and flanking size.")
    parser.add_argument("-r", "--reference", required=True, help="Indexed Reference FASTA file")
    parser.add_argument("-c", "--chrom", required=True, help="Chromosome")
    
    # NEW ARGUMENTS: Center + Flank
    parser.add_argument("-p", "--position", type=int, required=True, help="Center genomic coordinate (1-based)")
    parser.add_argument("-f", "--flank", type=int, required=True, help="Flanking size on either side of the position")
    
    parser.add_argument("--prob-del", type=float, default=0.1, help="Probability of a read having a Deletion")
    parser.add_argument("--prob-dup", type=float, default=0.6, help="Probability of a read having a Duplication")
    
    parser.add_argument("--mean-insert", type=float, default=350, help="Mean insert size (default: 500)")
    parser.add_argument("--std-insert", type=float, default=100, help="Std dev of insert size (default: 50)")
    parser.add_argument("--read-count", type=int, default=10000, help="Total number of read pairs to simulate")
    
    # SV Size Control
    parser.add_argument("--mean-sv-size", type=float, default=100, help="Mean (scale) for SV size (exponential)")
    parser.add_argument("--min-sv-size", type=int, default=1, help="Hard minimum SV size")
    parser.add_argument("--max-sv-size", type=int, default=5000, help="Hard maximum SV size")
    
    # SV Position Control
    parser.add_argument("--mean-offset", type=float, default=3.0, help="Mean distance of SV center from the target Position (exponential)")
    
    parser.add_argument("--out-prefix", default="simulated_flanked", help="Output file prefix")
    return parser.parse_args()

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def get_exponential_sv_size(mean, min_val, max_val):
    while True:
        val = int(np.random.exponential(scale=mean))
        if min_val <= val <= max_val:
            return val

def generate_mutated_sequence(ref_seq, sv_type, args, region_start):
    region_len = len(ref_seq)
    buffer = 500
    
    # 1. Determine SV Size
    effective_max = min(args.max_sv_size, region_len - (2 * buffer))
    if effective_max < args.min_sv_size:
        return ref_seq, 0, 0, 0
        
    sv_size = get_exponential_sv_size(args.mean_sv_size, args.min_sv_size, effective_max)
    
    # 2. Determine Position
    # Center logic: We want the SV to be near the MIDDLE of the fetched seq.
    midpoint = region_len // 2
    offset_dist = int(np.random.exponential(scale=args.mean_offset))
    direction = random.choice([-1, 1])
    
    sv_center = midpoint + (direction * offset_dist)
    
    rel_start = sv_center - (sv_size // 2)
    rel_end = rel_start + sv_size
    
    # 3. Boundary Checks (Clamping)
    min_valid_start = buffer
    max_valid_start = region_len - buffer - sv_size
    
    if rel_start < min_valid_start: rel_start = min_valid_start
    if rel_start > max_valid_start: rel_start = max_valid_start
    rel_end = rel_start + sv_size
    
    # Absolute coords
    abs_start = region_start + rel_start
    abs_end = region_start + rel_end
    
    # 4. Construct Sequence
    if sv_type == 'DEL':
        new_seq = ref_seq[:rel_start] + ref_seq[rel_end:]
    elif sv_type == 'DUP':
        dup_block = ref_seq[rel_start:rel_end]
        new_seq = ref_seq[:rel_end] + dup_block + ref_seq[rel_end:]
    else:
        new_seq = ref_seq
        abs_start = 0
        abs_end = 0
        sv_size = 0

    return new_seq, abs_start, abs_end, sv_size

def main():
    args = parse_args()
    
    # --- Calculate Start/End from Position/Flank ---
    region_start = max(1, args.position - args.flank)
    region_end = args.position + args.flank
    
    print(f"Target: {args.chrom}:{args.position} +/- {args.flank}bp")
    print(f"Fetching region: {args.chrom}:{region_start}-{region_end}")

    prob_wt = 1.0 - (args.prob_del + args.prob_dup)
    if prob_wt < 0: prob_wt = 0.0

    fasta = pysam.FastaFile(args.reference)
    try:
        ref_seq = fasta.fetch(args.chrom, region_start - 1, region_end)
    except ValueError:
        sys.exit(f"Error fetching region")
        
    r1_fn = f"{args.out_prefix}_R1.fastq.gz"
    r2_fn = f"{args.out_prefix}_R2.fastq.gz"
    read_len = 151
    
    print(f"Simulating {args.read_count} reads...")
    print(f"Probabilities -> WT: {prob_wt:.2f}, DEL: {args.prob_del:.2f}, DUP: {args.prob_dup:.2f}")

    with gzip.open(r1_fn, 'wt') as f1, gzip.open(r2_fn, 'wt') as f2:
        for i in range(args.read_count):
            
            choice = random.choices(['WT', 'DEL', 'DUP'], weights=[prob_wt, args.prob_del, args.prob_dup])[0]
            
            if choice == 'WT':
                current_seq = ref_seq
                sv_info = "WT:0-0"
            else:
                current_seq, sv_s, sv_e, sv_len = generate_mutated_sequence(ref_seq, choice, args, region_start)
                sv_info = f"{choice}:{sv_s}-{sv_e}"

            valid_frag = False
            attempts = 0
            seq_len = len(current_seq)
            
            while not valid_frag and attempts < 10:
                insert_size = int(np.random.normal(args.mean_insert, args.std_insert))
                if insert_size < read_len: insert_size = read_len
                
                if seq_len > insert_size:
                    frag_start = random.randint(0, seq_len - insert_size)
                    frag_end = frag_start + insert_size
                    fragment = current_seq[frag_start:frag_end]
                    valid_frag = True
                attempts += 1
            
            if not valid_frag: continue
            
            seq_r1 = fragment[:read_len]
            seq_r2 = reverse_complement(fragment[-read_len:])
            
            if len(seq_r1) < read_len: seq_r1 = seq_r1.ljust(read_len, 'N')
            if len(seq_r2) < read_len: seq_r2 = seq_r2.ljust(read_len, 'N')

            qual = 'I' * read_len
            header_base = f"@{i}:{args.chrom}:{region_start}-{region_end}_{sv_info}"
            
            f1.write(f"{header_base}/1\n{seq_r1}\n+\n{qual}\n")
            f2.write(f"{header_base}/2\n{seq_r2}\n+\n{qual}\n")

    print(f"Done. Wrote to {r1_fn} and {r2_fn}")

if __name__ == "__main__":
    main()