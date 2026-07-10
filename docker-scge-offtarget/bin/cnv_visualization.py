#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pyranges as pr
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde
import argparse, gzip, random, sys

__version__ = '1.0.0'

def find_local_maxima(scores):
    # Convert scores to a numpy array if it's not already
    scores_array = np.array(scores)
    kde = gaussian_kde(scores)
    x_grid = np.linspace(0, 1, 1000)
    density = kde.evaluate(x_grid)

    # Finding peaks in the density distribution
    peaks_indices = find_peaks(density,distance=10)[0]
    peak_values = [ round(x,4) for x in x_grid[peaks_indices].tolist() ]

    return peak_values

def read_bedgraph(file_path, window_size=500000):

    def process_window(vafs, chromosome, start, end):
        if len(vafs)>20:
            maxima = find_local_maxima(vafs)
            if isinstance(maxima,list):
                return ([chromosome, start, end, ','.join([str(i) for i in maxima])])
            else:
                return ([chromosome, start, end, None])
        else:
            return ([chromosome, start, end, None])

    ranges = []

    open_func = gzip.open if file_path.endswith('.gz') else open
    with open_func(file_path, 'rt') as file:
        window_data = []
        prev_chr = None
        window_start = 0
        end = 0

        for line in file:
            parts = line.strip().split()
            if len(parts) == 4:
                chromosome, start, end, score = parts
                start, end, score = int(start), int(end), float(score)

                if prev_chr is None:
                    prev_chr = chromosome

                if chromosome != prev_chr or end > window_start + window_size:
                    row = process_window(window_data, prev_chr, window_start, window_start + window_size)
                    if row:
                        ranges.append(row)

                    window_data = []
                    if chromosome != prev_chr:
                        prev_chr = chromosome
                        window_start = 0
                    else:
                        window_start = window_start + window_size

                window_data.append(score)

        # Process the last window after EOF
        row = process_window(window_data, prev_chr, window_start, window_start + window_size)
        if row:
            ranges.append(row)

        return pd.DataFrame(ranges,columns=['Chromosome','Start','End','VAF'])


def main():
    parser = argparse.ArgumentParser(description="Perform KMeans clustering on BEDGraph data.")
    parser.add_argument("-o", "--outfile", type=str, help="Path to the BEDGraph file (can be gzipped)")
    parser.add_argument("vaf_file", type=str, help="Path to the BEDGraph VAF file (can be gzipped)")
    parser.add_argument("cov_file", type=str, help="Path to the normalized coverage file (can be gzipped)")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s: ' + __version__)

    args = parser.parse_args()

    vafs = read_bedgraph(args.vaf_file)
    vafs = pr.PyRanges(vafs)
    gr = pd.read_csv(args.cov_file, sep='\t', skiprows=3, names=['Chromosome','Start','End','Name','NormalizedCoverage','Pairs'])
    gr = pr.PyRanges(gr)
    df = vafs.join(gr).df
    df['Chromosome'] = df['Chromosome'].astype(str)

    df = df.groupby(['Chromosome', 'Start', 'End']).agg({'VAF':'first','NormalizedCoverage': 'mean'}).reset_index() 
    df = pr.PyRanges(df).sort().df.reset_index()

    if args.outfile:
        df.to_csv(args.outfile,sep="\t",index=False)
    else:
        df.to_csv(sys.stdout,sep="\t",index=False)
        

if __name__ == "__main__":
    main()

    