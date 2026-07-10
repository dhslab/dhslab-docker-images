#!/usr/bin/env python3

import argparse
import os
import sys
import pandas as pd
import re

__version__ = "1.0.0"

def process_paired_format(df, input_dir):
    """
    Processes a DataFrame where Read 1 and Read 2 are in separate columns.
    Returns a new DataFrame in the DRAGEN fastq_list.csv format.
    """
    new_rows = []
    
    # Clean up column names by stripping whitespace
    df.columns = df.columns.str.strip()
    
    # Define required columns for this format
    required_cols = ['FASTQ Path - Read 1', 'FASTQ Path - Read 2', 'Flowcell Lane', 'Library Name', 'Index Sequence', 'Flowcell ID']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column for paired format: '{col}'")
            
    for _, row in df.iterrows():
        read1_filename = row['FASTQ Path - Read 1']
        read2_filename = row['FASTQ Path - Read 2']
        
        read1_path = os.path.join(input_dir, read1_filename)
        read2_path = os.path.join(input_dir, read2_filename)
        
        # --- Validation: Check if files exist ---
        if not os.path.exists(read1_path):
            print(f"Warning: Read 1 file not found, skipping row: {read1_path}", file=sys.stderr)
            continue
        if not os.path.exists(read2_path):
            print(f"Warning: Read 2 file not found, skipping row: {read2_path}", file=sys.stderr)
            continue

        # --- Data Extraction ---
        lane = row['Flowcell Lane']
        rgsm = row['Library Name']
        flowcell = row['Flowcell ID']
        
        index_seq = str(row['Index Sequence'])
        indices = [i for i in index_seq.split('-') if i]
        index1 = indices[0] if len(indices) > 0 else ''
        index2 = indices[1] if len(indices) > 1 else ''

        # --- Construct new fields ---
        rgid_parts = [flowcell, index1, index2, str(lane)]
        rglb_parts = [rgsm, index1, index2]
        
        rgid = '.'.join(filter(None, rgid_parts))
        rglb = '.'.join(filter(None, rglb_parts))

        new_rows.append({
            'RGID': rgid,
            'RGSM': rgsm,
            'RGLB': rglb,
            'Lane': lane,
            'Read1File': read1_path,
            'Read2File': read2_path
        })
        
    return pd.DataFrame(new_rows)

def process_single_format(df, input_dir):
    """
    Processes a DataFrame where each read is on its own row.
    Groups R1 and R2 reads and returns a DataFrame in the DRAGEN format.
    """
    # Clean up column names
    df.columns = df.columns.str.strip()

    required_cols = ['FASTQ', 'Flowcell Lane', 'Library Name', 'Index Sequence', 'Flowcell ID']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column for single format: '{col}'")

    # --- Data Extraction and Grouping Key Creation ---
    # Extract read number and create a key to group R1/R2 pairs
    df['read_num'] = df['FASTQ'].apply(lambda x: 'R1' if '_R1_' in x else ('R2' if '_R2_' in x else 'Unknown'))
    df = df[df['read_num'] != 'Unknown'] # Filter out non-R1/R2 files
    
    df['grouping_key'] = df['Library Name'].astype(str) + '_' + df['Flowcell Lane'].astype(str)
    
    # --- Pivot the data to get R1 and R2 in the same row ---
    pivoted = df.pivot_table(
        index=['grouping_key', 'Library Name', 'Flowcell Lane', 'Index Sequence', 'Flowcell ID'],
        columns='read_num',
        values='FASTQ',
        aggfunc='first'
    ).reset_index()

    pivoted = pivoted.rename(columns={'R1': 'ReadFile1_fname', 'R2': 'ReadFile2_fname'})

    # Now process the pivoted data like the paired format
    new_rows = []
    for _, row in pivoted.iterrows():
        read1_filename = row['ReadFile1_fname']
        read2_filename = row['ReadFile2_fname']

        if pd.isna(read1_filename) or pd.isna(read2_filename):
            print(f"Warning: Missing R1 or R2 for group {row['grouping_key']}, skipping.", file=sys.stderr)
            continue

        read1_path = os.path.join(input_dir, read1_filename)
        read2_path = os.path.join(input_dir, read2_filename)

        if not os.path.exists(read1_path) or not os.path.exists(read2_path):
            print(f"Warning: One or both FASTQ files not found for group {row['grouping_key']}, skipping.", file=sys.stderr)
            continue

        # --- Data Extraction ---
        lane = row['Flowcell Lane']
        rgsm = row['Library Name']
        flowcell = row['Flowcell ID']

        index_seq = str(row['Index Sequence'])
        indices = [i for i in index_seq.split('-') if i]
        index1 = indices[0] if len(indices) > 0 else ''
        index2 = indices[1] if len(indices) > 1 else ''
        
        # --- Construct new fields ---
        rgid_parts = [flowcell, index1, index2, str(lane)]
        rglb_parts = [rgsm, index1, index2]

        rgid = '.'.join(filter(None, rgid_parts))
        rglb = '.'.join(filter(None, rglb_parts))

        new_rows.append({
            'RGID': rgid,
            'RGSM': rgsm,
            'RGLB': rglb,
            'Lane': lane,
            'Read1File': read1_path,
            'Read2File': read2_path
        })
        
    return pd.DataFrame(new_rows)


def main():
    parser = argparse.ArgumentParser(
        description="Converts a custom sequencing stats CSV into a standard Illumina fastq_list.csv format.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input_csv", help="Path to the input CSV file from the sequencer.")
    parser.add_argument("output_csv", help="Path for the output fastq_list.csv file.")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    
    args = parser.parse_args()

    try:
        # Determine the directory of the input file to prepend to FASTQ paths
        input_dir = os.path.dirname(os.path.abspath(args.input_csv))
        
        print(f"Reading input file: {args.input_csv}", file=sys.stderr)
        df = pd.read_csv(args.input_csv)
        
        # --- Auto-detect file format ---
        header = {h.strip() for h in df.columns}
        if 'FASTQ Path - Read 1' in header and 'FASTQ Path - Read 2' in header:
            print("Detected paired-column format.", file=sys.stderr)
            final_df = process_paired_format(df, input_dir)
        elif 'FASTQ' in header:
            print("Detected single-column (long) format.", file=sys.stderr)
            final_df = process_single_format(df, input_dir)
        else:
            raise ValueError("Unrecognized input format. Could not find 'FASTQ' or 'FASTQ Path - Read 1' columns.")

        # --- Write Final Output ---
        output_columns = ['RGID', 'RGSM', 'RGLB', 'Lane', 'Read1File', 'Read2File']
        
        if not final_df.empty:
            final_df = final_df[output_columns] # Ensure correct column order
            print(f"Writing {len(final_df)} records to {args.output_csv}", file=sys.stderr)
            final_df.to_csv(args.output_csv, index=False)
            print("Conversion complete.", file=sys.stderr)
        else:
            print("Warning: No valid records were processed to write to the output file.", file=sys.stderr)

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()

