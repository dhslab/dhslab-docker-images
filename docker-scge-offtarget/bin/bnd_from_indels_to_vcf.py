#!/usr/bin/env python3
import sys
import csv
import argparse
from collections import defaultdict

def breakend_alleles(ref_base, mate_chr, mate_pos, strands):
    # strands: '++' or '+-' per extract_variant_reads.py
    # Use simple VCF BND allele conventions
    # For simplicity, use ref_base = 'N' if unknown
    if not ref_base or ref_base == '.':
        ref_base = 'N'
    mate = f"{mate_chr}:{mate_pos}"
    if strands == '++':
        # right clip -> A]chr:pos]
        alt = f"{ref_base}]{mate}]"
    else:
        # '+-' use A[chr:pos[
        alt = f"{ref_base}[{mate}["
    return ref_base, alt

def main(args):
    vcf_path = args.outfile
    with open(args.indels_path) as fh:
        header = fh.readline().strip().split('\\t')
        hdr = {k:i for i,k in enumerate(header)}
        # No BNDs: emit empty VCF header and exit
        records = []
        for line in fh:
            cols = line.rstrip('\\n').split('\\t')
            if len(cols) < len(header):
                continue
            bnd_count = cols[hdr.get('bnd_count', -1)] if 'bnd_count' in hdr else '0'
            if bnd_count in ('0', 0, '', None):
                continue
            bnd_info = cols[hdr.get('bnd_info', -1)] if 'bnd_info' in hdr else '.'
            if not bnd_info or bnd_info == '.':
                continue
            for key in bnd_info.split(';'):
                # chrom:pos:chrom2:pos2:strands:counts:control_alt_counts:Distance
                parts = key.split('|')
                if len(parts) < 5:
                    continue
                chrom, pos, chrom2, pos2, strands = parts[:5]
                try:
                    pos = int(pos)
                    pos2 = int(pos2)
                except Exception:
                    continue
                records.append((chrom, pos, chrom2, pos2, strands))

    # build paired mates deterministically by index
    lines = []
    for idx, (c1, p1, c2, p2, st) in enumerate(records, start=1):
        id1 = f"BND{idx}_1"
        id2 = f"BND{idx}_2"
        r1, a1 = breakend_alleles('N', c2, p2, st)
        # reverse mate strands for reciprocal allele rendering
        st_mate = '++' if st == '+-' else '+-'
        r2, a2 = breakend_alleles('N', c1, p1, st_mate)
        info1 = f"SVTYPE=BND;MATEID={id2};EVENT=EVT{idx}"
        info2 = f"SVTYPE=BND;MATEID={id1};EVENT=EVT{idx}"
        lines.append((c1, p1, id1, r1, a1, '.', 'PASS', info1))
        lines.append((c2, p2, id2, r2, a2, '.', 'PASS', info2))

    with open(vcf_path, 'w') as out:
        out.write("##fileformat=VCFv4.2\\n")
        out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\\n')
        out.write('##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakend">\\n')
        out.write('##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of associated event">\\n')
        out.write("#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n")
        for c, p, i, r, a, q, f, info in lines:
            out.write(f"{c}\\t{p}\\t{i}\\t{r}\\t{a}\\t{q}\\t{f}\\t{info}\\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert indels.txt to BND VCF.')
    parser.add_argument('--meta_id', required=True, help='Sample ID')
    parser.add_argument('--indels_path', required=True, help='Path to indels.txt file')
    # add output file arg
    parser.add_argument('--outfile', required=True, help='Path to output VCF file')
    args = parser.parse_args()
    main(args)
