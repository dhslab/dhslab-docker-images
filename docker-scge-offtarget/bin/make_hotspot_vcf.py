#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import pyranges as pr
import pysam
from pathlib import Path
from typing import List

# Define the script version
__version__ = "1.0.0"

# --- Helper Functions ---

def check_file(path: str) -> Path:
    """Checks if a file path exists and is a file."""
    p = Path(path)
    if not p.is_file():
        raise argparse.ArgumentTypeError(f"File not found: {path}")
    return p

def add_sequence_column(row: pd.Series, fasta_handle: pysam.FastaFile) -> str:
    """
    Fetches a genomic sequence for a given row using an open pysam.FastaFile handle.
    """
    try:
        # pysam.fetch is 0-based, and BED Start is 0-based.
        # It fetches the interval [start, end), which matches the BED standard.
        sequence = fasta_handle.fetch(row["Chromosome"], row["Start"], row["End"])
        return sequence
    except (ValueError, KeyError) as e:
        print(f"Warning: Could not fetch sequence for {row['Chromosome']}:{row['Start']}-{row['End']}. Reason: {e}", file=sys.stderr)
        return "N" * (row['End'] - row['Start'])


def dataframe_to_vcf(df: pd.DataFrame, fasta_handle: pysam.FastaFile, outfile_path: str) -> None:
    """
    Converts a DataFrame with genomic positions into a VCF file using pysam.
    Records are sorted based on the contig order in the provided fasta_handle.
    
    Args:
        df: DataFrame containing 'Chromosome', 'Position', and 'Sequence' columns.
        fasta_handle: Open pysam.FastaFile object.
        outfile_path: Path to write the output VCF.
    """
    # 1. Enforce Sorting Order based on Reference FASTA
    # We create a categorical type for the 'Chromosome' column using the 
    # ordered list of references from the pysam FastaFile.
    df_sorted = df.copy()
    df_sorted['Chromosome'] = pd.Categorical(
        df_sorted['Chromosome'], 
        categories=fasta_handle.references, 
        ordered=True
    )
    
    # Sort by Chromosome (ref order) then Position. 
    # DropNA ensures we don't crash on chroms not in the reference.
    df_sorted = df_sorted.sort_values(by=['Chromosome', 'Position']).dropna(subset=['Chromosome'])

    # 2. Create the VCF Header
    header = pysam.VariantHeader()
    
    # Add contigs to header with lengths from the FASTA index
    for contig in fasta_handle.references:
        length = fasta_handle.get_reference_length(contig)
        header.contigs.add(contig, length=length)
    
    # 3. Write directly to the specified output file
    # pysam.VariantFile handles opening the file for writing.
    with pysam.VariantFile(outfile_path, 'w', header=header) as vcf_out:
        for _, row in df_sorted.iterrows():
            # Create a new record object
            rec = vcf_out.new_record()
            
            # Fill fields
            rec.chrom = str(row['Chromosome'])
            rec.pos = int(row['Position'])
            rec.id = '.'
            rec.ref = row['Sequence']
            rec.alts = ('N')
            rec.filter.add('PASS')
            
            # Write record
            vcf_out.write(rec)


# --- Main Application Logic ---

def main():
    """
    Main function to generate a hotspot VCF from a BED file of genomic regions.
    """
    parser = argparse.ArgumentParser(
        description="Prepare a hotspot VCF from a BED file of genomic regions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--bed", type=check_file, required=True, help="BED file containing hotspot regions.")
    parser.add_argument("--targets", type=check_file, required=True, help="Editing target file.")
    parser.add_argument("--window",type=int,default=100,help="Window around editing targets to include in VCF file")
    parser.add_argument("--fasta", type=check_file, required=True, help="Path to the indexed reference genome FASTA file.")
    parser.add_argument("--outfile", required=True, help="Output VCF file.")
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args()

    try:
        # --- 1. Open FASTA file once ---
        print(f"Opening reference FASTA: {args.fasta}", file=sys.stderr)
        fasta_handle = pysam.FastaFile(str(args.fasta))

        # --- 2. Read and Merge Genomic Regions ---
        print(f"Reading and processing BED file: {args.bed}", file=sys.stderr)
        bed_df = pd.read_csv(args.bed, sep='\t', usecols=[0, 1, 2], names=['Chromosome', 'Start', 'End'])
        bed_df['Start'] = bed_df['Start'] + 1

        # Read in editing targets in VCF format
        print(f"Reading editing targets: {args.targets}", file=sys.stderr)
        vcf_in = pysam.VariantFile(args.targets)
        vcf_data = [(record.chrom, record.pos) for record in vcf_in]
        targets_df = pd.DataFrame(vcf_data, columns=['Chromosome', 'Start'])
        targets_df['Start'] = targets_df['Start'] - args.window
        targets_df['End'] = targets_df['Start'] + (args.window)

        bed_df = pd.concat([bed_df,targets_df],ignore_index=True)

        merged_df = pr.PyRanges(bed_df).merge().sort().df

        # --- 3. Fetch Sequences for Merged Regions ---
        print("Fetching sequences for merged regions...", file=sys.stderr)
        merged_df['sequences'] = merged_df.apply(
            lambda row: add_sequence_column(row, fasta_handle), 
            axis=1
        )

        # --- 4. Expand Regions into a VCF-like DataFrame (Optimized) ---
        print("Expanding genomic regions into VCF format...", file=sys.stderr)
        all_chroms, all_positions, all_sequences = [], [], []

        for _, row in merged_df.iterrows():
            num_bases = len(row['sequences'])
            all_chroms.extend([row['Chromosome']] * num_bases)
            all_positions.extend(range(row['Start'], row['End']))
            all_sequences.extend(list(row['sequences']))

        vcf_df = pd.DataFrame({
            'Chromosome': all_chroms,
            'Position': all_positions,
            'Sequence': all_sequences
        })

        # --- 5. Convert to VCF and Write to File ---
        dataframe_to_vcf(vcf_df, fasta_handle, args.outfile)

    except Exception as e:
        sys.exit(f"An error occurred: {e}")
    finally:
        # Ensure the FASTA file handle is closed
        if 'fasta_handle' in locals() and fasta_handle:
            fasta_handle.close()

if __name__ == "__main__":
    main()

