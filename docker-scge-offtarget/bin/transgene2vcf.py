#!/usr/bin/env python3

import sys
import argparse
import os
import pandas as pd

__version__ = "1.2.0"

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert a transgene TSV file into a VCF file."
    )
    parser.add_argument(
        "transgene_path", 
        help="Path to the input transgene TSV file."
    )
    parser.add_argument(
        "output_file", 
        help="Path for the output VCF file."
    )
    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show the version number and exit."
    )
    return parser.parse_args()

def write_vcf(path, records=None):
    """
    Writes VCF header and records to path. 
    If records is None or empty, only the header is written.
    """
    if records is None:
        records = []
        
    try:
        with open(path, 'w') as out:
            # Write VCF Header
            out.write("##fileformat=VCFv4.2\n")
            out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

            # Write Records
            for idx, (chrom, pos, ref, alt) in enumerate(records, start=1):
                variant_id = f"TRANSGENE{idx}"
                qual = "."
                filter_val = "PASS"
                info = "SVTYPE=TRANSGENE"
                out.write(f"{chrom}\t{pos}\t{variant_id}\t{ref}\t{alt}\t{qual}\t{filter_val}\t{info}\n")
                
    except IOError as e:
        sys.exit(f"Error writing to output file '{path}': {e}")

def main():
    args = parse_args()

    transgene_path = args.transgene_path
    vcf_path = args.output_file

    # Check if input file exists
    if not os.path.exists(transgene_path):
        sys.exit(f"Error: Input file '{transgene_path}' not found.")

    try:
        # Attempt to read the file using pandas
        # This automatically handles empty files by raising EmptyDataError
        df = pd.read_csv(transgene_path, sep='\t')
    
    except pd.errors.EmptyDataError:
        # File is empty (or only contains whitespace)
        write_vcf(vcf_path)
        return
    except Exception as e:
        sys.exit(f"Error reading input file: {e}")

    # Check if DataFrame is logically empty (has header but no rows)
    if df.empty:
        write_vcf(vcf_path)
        return

    # Clean up column names (strip whitespace)
    df.columns = df.columns.str.strip()

    # Handle common header variations (e.g. #Chromosome -> Chromosome)
    if '#Chromosome' in df.columns and 'Chromosome' not in df.columns:
        df.rename(columns={'#Chromosome': 'Chromosome'}, inplace=True)

    # Validate required columns exist
    required_cols = ['Chromosome', 'Start']
    missing_cols = [col for col in required_cols if col not in df.columns]

    if missing_cols:
        found_headers = ", ".join(f"'{h}'" for h in df.columns)
        msg = f"Error: Required column(s) {missing_cols} missing from input file header.\nFound headers: {found_headers}"
        if len(df.columns) <= 1:
            msg += "\nHint: The header contains very few columns. Ensure the file is Tab-Separated (TSV), not Comma-Separated (CSV)."
        sys.exit(msg)

    # Extract records
    records = []
    try:
        for _, row in df.iterrows():
            chrom = row['Chromosome']
            pos = row['Start']
            ref = 'N'
            alt = '<INS>'
            records.append((chrom, pos, ref, alt))
    except Exception as e:
        print(f"Warning: Error processing a row: {e}", file=sys.stderr)

    # Write populated VCF
    write_vcf(vcf_path, records)

if __name__ == '__main__':
    main()