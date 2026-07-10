#!/usr/bin/env python3

import json, sys, os, csv, argparse, gzip, re
from collections import Counter
import pandas as pd

__version__ = '1.0.0'

def revcomp(dna):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}  # DNA complement pairs
    reverse_complement = "".join(complement.get(base, base) for base in reversed(dna))
    return reverse_complement

def checkfile(file_path):
    """Check if a file exists at the given path."""
    if not os.path.exists(file_path):
        raise argparse.ArgumentTypeError(f"The file {file_path} does not exist.")
    return file_path

def process_fastq_list(row):
    row = row.to_dict()

    df = pd.read_csv(row['fastq_list'])
    df = df[df['RGSM']==row['id']]
    row['read1'] = df['Read1File'].tolist()
    row['read2'] = df['Read2File'].tolist()
    
    del row['fastq_list']

    # Dummy implementation, replace with actual file processing
    return pd.Series(row)

def process_reads(row):

    row = row.to_dict()

    read1path = row['read1']
    read2path = row['read2']
    
    with gzip.open(read1path, 'rt') as gz_file:
        readName = gz_file.readline().strip()

    read_header_parts = readName.split(':')

    # if indexes are present in the read header, get the first 1000
    # and find the most common one (to account for mismatches/errors in index read)
    index1 = 'UNKNOWN'
    index2 = 'UNKNOWN'
    indexes = readName.split(' ')
    if len(indexes) > 1:
        indexlist = []
        with gzip.open(read1path, 'rt') as file:
            for i, line in enumerate(file):
                if i % 4 == 0:  # Read names are on every 4th line starting from 0
                    read_name = line.strip()
                    indexes = read_name.split(' ')
                    if len(indexes) > 1:
                        index = indexes[1].split(':')[-1]
                        indexlist = indexlist + [index]
        
                if i == 3999:  # Stop after reading the first 100 entries (0 to 399)
                    break

        counter = Counter(indexlist)
        most_common_index, _ = counter.most_common(1)[0]
        index1, index2 = most_common_index.split('+')

        outdict = row
        outdict['i7index'] = row.get('i7index') or index1
        outdict['i5index'] = row.get('i5index') or revcomp(index2)
        outdict['flowcell'] = row.get('flowcell') or read_header_parts[2]
        outdict['lane'] = row.get('lane') or read_header_parts[3]
        
        return(pd.Series(outdict))

def reformat_mgi_samplesheet(file,filepath):

    # get path of input file if no path was passed (read files are basenames only)
    if filepath is None:
        filepath = os.path.dirname(file)

    df = pd.read_csv(file, sep=None, engine='python')

    if 'FASTQ' in df.columns: # 1 read file per line format

        df['FASTQ'] = df['FASTQ'].apply(lambda x: os.path.join(filepath,x))

        df['Read'] = df['FASTQ'].apply(lambda x: 'R1' if '_R1_' in x else 'R2')
        df_read1 = df[df['Read']=='R1'].copy().rename(columns={'FASTQ':'read1'})
        df_read2 = df[df['Read']=='R2'].copy().rename(columns={'FASTQ':'read2'})
        
        df = df_read1.merge(df_read2,on=['Library Name','Flowcell ID','Index Sequence','Flowcell Lane'])

    elif 'FASTQ Path - Read 1' in df.columns: # paired reads per line format

        df['read1'] = df['FASTQ Path - Read 1'].apply(lambda x: os.path.join(filepath,x))
        df['read2'] = df['FASTQ Path - Read 2'].apply(lambda x: os.path.join(filepath,x))

        df = df[['read1','read2','Library Name','Flowcell ID','Index Sequence','Flowcell Lane']]

    df[['i7index', 'i5index']] = df['Index Sequence'].str.split('-', expand=True)
    df[['uid', 'sample','assay']] = df['Library Name'].str.split('-', expand=True)[[0,1,2]]
    df[['id','flowcell', 'lane']] = df[['Library Name','Flowcell ID', 'Flowcell Lane']]
    df = df[['id','uid','sample_type','sample_id','assay','i7index', 'i5index','flowcell', 'lane','read1','read2']]

    return df


def main():

    parser = argparse.ArgumentParser(description="Validate samplesheet for dragen multiworkflow")
    parser.add_argument('-m', '--mgi', action='store_true', help="MGI samplesheet style")
    parser.add_argument('-d', '--dir', type=str, default=None,help="File directory location for read1 and read2")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s: ' + __version__)
    parser.add_argument('samplesheet',type=str, help="Samplesheet to validate")

    args = parser.parse_args()

    df = pd.DataFrame()

    df_reads = pd.DataFrame()
    df_fastqlist = pd.DataFrame()
    df_bam = pd.DataFrame()
    df_cram = pd.DataFrame()
    df_dragen = pd.DataFrame()

    if args.mgi:
        df_reads = reformat_mgi_samplesheet(args.samplesheet,args.dir)

    else:
        df = pd.read_csv(args.samplesheet, sep=None, engine='python')

        # Check 'id' field
        if df['id'].isnull().any():
            raise ValueError("'id' field is required and cannot contain null values.")

        # Process 'fastq_list' by getting the reads
        if 'fastq_list' in df.columns and df['fastq_list'].notna().any():
            df_fastqlist = df[df['fastq_list'].notna()].copy().drop(columns=['read1','read2','cram','bam'],errors='ignore')
            df_fastqlist = df_fastqlist.apply(lambda row: process_fastq_list(row),axis=1)
            df_fastqlist['reads'] = df_fastqlist.apply(lambda row: list(zip(row['read1'], row['read2'])), axis=1)
            df_fastqlist = df_fastqlist.explode('reads')
            df_fastqlist['read1'], df_fastqlist['read2'] = zip(*df_fastqlist['reads'])
            df_fastqlist.drop(columns=['reads'], inplace=True)

            df = pd.concat([df[df['fastq_list'].isna()].copy(),df_fastqlist],axis=0,ignore_index=True)

        if 'read1' in df.columns:
            df_reads = df[df['read1'].notna()].copy().drop(columns=['fastq_list','cram','bam'],errors='ignore')
            df_reads = df_reads.apply(lambda row: process_reads(row),axis=1)

        # Process 'cram'
        if 'cram' in df.columns and df['cram'].notna().any():
            df_cram = df[df['cram'].notna()].copy().drop(columns=['read1','read2','fastq_list','bam'],errors='ignore')

        # Process 'bam'
        if 'bam' in df.columns and df['bam'].notna().any():
            df_bam = df[df['bam'].notna()].copy().drop(columns=['read1','read2','fastq_list','cram'],errors='ignore')

        if 'dragen_path' in df.columns and df['dragen_path'].notna().any():
            df_dragen = df[df['dragen_path'].notna()].copy().drop(columns=['read1','read2','fastq_list','cram','bam'],errors='ignore')

    df = pd.concat([df_reads,df_cram,df_bam,df_dragen],axis=0,ignore_index=True)

    # Check file existence for files that should be staged or mounted
    for col in ['cram', 'bam', 'read1', 'read2']:
        if col in df.columns:
            if df[col].notna().any() and not df[df[col].notna()][col].apply(os.path.exists).all():
                raise FileNotFoundError(f"Some files in the '{col}' column do not exist.")

    # 'dragen_path' may be an absolute host path not mounted inside the container.
    # Do not hard-fail here; emit a warning and proceed. Downstream processes will stage required files.
    if 'dragen_path' in df.columns and df['dragen_path'].notna().any():
        exists_mask = df[df['dragen_path'].notna()]['dragen_path'].apply(os.path.exists)
        if not exists_mask.all():
            missing_count = (~exists_mask).sum()
            print(
                f"Warning: {missing_count} 'dragen_path' entries are not accessible in the current environment. "
                "Proceeding, but ensure downstream steps can access DRAGEN outputs.",
                file=sys.stderr,
            )

    if df.shape[0] == 0:
        sys.exit(f"No rows in samplesheet {args.samplesheet} after checking.")

    # Write the DataFrame to a CSV file
    output_filename = 'samplesheet.valid.csv'
    df.to_csv(output_filename, index=False)

if __name__ == "__main__":
    main()

