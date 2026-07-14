#!/usr/bin/env python3

import csv
import argparse

parser = argparse.ArgumentParser(description='Transform transgene data.')
parser.add_argument('--input', type=str, help='Input file path')
parser.add_argument('--output', type=str, help='Output file path')
args = parser.parse_args()

input_file = args.input
output_file = args.output

# Open the input file and process it
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    reader = csv.DictReader(infile, delimiter='\t')
    outfile.write("#chr start end [options]\n")  # Write the header for the output file

    for row in reader:
        # Extract chromosome, start, and end positions
        chromosome = row['Chromosome'].replace('chr', 'hs')  # Convert 'chr' to 'hs'
        start = row['Start']
        end = row['End']
        annotation = row['Var']  # Use the 'Var' column for annotation

        # Write the transformed row to the output file
        outfile.write(f"{chromosome}\t{start}\t{end}\tannotation={annotation}\n")
