#!/usr/bin/env python3

import sys
import os
import csv
import gzip

def extract_info_from_fastq(fastq_file):
    if fastq_file.endswith('.gz'):
        # If the file is gzipped, open it with gzip.open
        with gzip.open(fastq_file, 'rt') as f:
            first_line = f.readline().strip()
    else:
        # If the file is not gzipped, open it normally
        with open(fastq_file, 'r') as f:
            first_line = f.readline().strip()

    parts = first_line.split(':')
    flowcell_id = parts[2]
    lane_number = parts[3]
    index_sequence = parts[-1].split('+')
    index1 = index_sequence[0]
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    seq = index_sequence[1]
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    index = index1 + '-' + reverse_complement

    return flowcell_id, lane_number, index

def main():
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <sample_sheet_file>")
        sys.exit(1)

    sample_sheet_file = sys.argv[1]

    if not os.path.isfile(sample_sheet_file):
        print(f"Error: File '{sample_sheet_file}' does not exist.")
        sys.exit(1)

    normal_fastqlist = []
    tumor_fastqlist = []

    with open(sample_sheet_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            sample_id = row['id']
            fastq1 = row['read1']
            fastq2 = row['read2']

            if not (os.path.isfile(fastq1) and os.path.isfile(fastq2)):
                print(f"Error: FASTQ files not found for sample ID {sample_id}")
                continue

            flowcell_id_1, lane_number_1, index_sequence_1 = extract_info_from_fastq(fastq1)
            flowcell_id_2, lane_number_2, index_sequence_2 = extract_info_from_fastq(fastq2)

            if flowcell_id_1 != flowcell_id_2 or lane_number_1 != lane_number_2:
                print(f"Error: FASTQ files for sample ID {sample_id} have different flowcell ID or lane number.")
                continue

            rgid = flowcell_id_1 + '.' + index_sequence_1 + '.' + lane_number_1
            rglb = sample_id + '.' + index_sequence_1
            rgsm = sample_id

            fastqlist_entry = [rgid, rglb, rgsm, lane_number_1, os.path.abspath(fastq1), os.path.abspath(fastq2)]

            if row['type'] == 'normal':
                normal_fastqlist.append(fastqlist_entry)
            elif row['type'] == 'tumor':
                tumor_fastqlist.append(fastqlist_entry)

    if normal_fastqlist:
        with open('normal_fastqlist.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["RGID", "RGLB", "RGSM", "Lane", "Read1File", "Read2File"])
            writer.writerows(normal_fastqlist)

    if tumor_fastqlist:
        with open('tumor_fastqlist.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["RGID", "RGLB", "RGSM", "Lane", "Read1File", "Read2File"])
            writer.writerows(tumor_fastqlist)

if __name__ == "__main__":
    main()
