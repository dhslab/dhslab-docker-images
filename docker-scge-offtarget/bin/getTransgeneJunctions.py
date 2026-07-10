#!/usr/bin/env python3

import argparse
import re
import sys
import pysam
import pandas as pd
import pyranges as pr

# reverse complement function
def revcomp(seq):
    tab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh') # maketrans <- maps the reverse complement
    return seq.translate(tab)[::-1] # translate(x)[::-1] <- works backward through the string, effectively reversing the string

import re

def cigar_to_aligned_positions(cigar_string):
    """Converts a CIGAR string to a set of aligned base pair positions, adjusting for leading soft clippings."""
    cigar_operations = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)

    aligned_positions = set()
    position = 0

    for length, operation in cigar_operations:
        length = int(length)
        if operation in ['S','H']:  # Adjust start position for clipping
            position += length
        elif operation in ['M', '=', 'X']:  # Match or mismatch
            for i in range(position, position + length):
                aligned_positions.add(i)
            position += length
        elif operation in ['D', 'N']:  # Deletion or skipped region: move position in reference
            position += length
        # Insertions (I), hard clippings (H), and padding (P) do not affect the position

    return aligned_positions

def count_cigar_overlap(cigar1, cigar2):
    """Counts the number of overlapping base pairs between two alignments."""
    positions1 = cigar_to_aligned_positions(cigar1)
    positions2 = cigar_to_aligned_positions(cigar2)

    overlap = positions1.intersection(positions2)
    return len(overlap)

# function to get indels for one amplicon from bam file
def get_chimeras(bam,contig,exclude=None,minSoftClip=20,minMq=1,maxMismatches=1):

    df = pd.DataFrame(columns=['Chromosome','Start','End','Var','Strand','Info'])

    if not contig in bam.references:
        return df
    
    # format of output: 
    # chr pos1 pos2 strand readname
    # iterate through once and get split reads and first end of discordant reads
    for read in bam.fetch(contig,multiple_iterators=True):
                
        # only consider primary alignments to the transgene sequence
        if read.is_duplicate is True or \
            read.is_secondary is True or \
            read.is_supplementary is True or \
            read.mapping_quality == 0:
            continue

        if exclude is not None and read.reference_start >= int(exclude.split(',')[0]) and read.reference_start <= int(exclude.split(',')[1]):
            continue
        
        # get left and right softclip lengths
        cigar = read.cigartuples
        leftSoftClip = 0
        rightSoftClip = 0
        if cigar[0][0] == 4:
            leftSoftClip = cigar[0][1]

        if cigar[-1][0] == 4:
            rightSoftClip = cigar[-1][1]
        # check for a chimeric read
        if read.is_proper_pair is False and read.mate_is_mapped is True and read.reference_name != read.next_reference_name:
            mate_pos = read.next_reference_start
            mate_strand = '-'
            # if mate is forward, then add 151 because the junction is at the end of the read
            if read.mate_is_forward:
                mate_pos = mate_pos + 151
                mate_strand = '+'

            mate = bam.mate(read)
            info = ['ID='+read.query_name, 'Type=PR', 'Read1Seq=' + read.query_sequence, 'Read2Seq=' + mate.query_sequence]

            df = pd.concat([df,pd.DataFrame([{'Chromosome':read.next_reference_name,'Start':mate_pos-1,'End':mate_pos,'Var':'INS','Strand':mate_strand,'Info':';'.join(info)}])]).reset_index(drop=True)

        # if this is a forward read and the left end is clipped then check for a supplementary alignment
        elif read.is_proper_pair and read.is_forward and leftSoftClip >= minSoftClip:
            if read.has_tag('SA'):
                for sa in read.get_tag('SA').rstrip(';').split(';'):
                    sChr, sPos, sStrand, sCigar, sMq, sNm = sa.split(',')
                    sPos = int(sPos)
                    if sChr != read.reference_name and int(sMq)>=minMq and int(sNm) <= maxMismatches:
                        readAligned = cigar_to_aligned_positions(read.cigarstring)
                        saAligned = cigar_to_aligned_positions(sCigar)
                        if len(readAligned.intersection(saAligned)) / len(readAligned) < 0.2:
                            if sStrand == '+':
                                sPos = sPos + len(saAligned)
                            
                            mateseq = '.'
                            if read.reference_name != read.next_reference_name:
                                mate = bam.mate(read)
                                mateseq = mate.query_sequence

                            info = ['ID='+read.query_name, 'Type=PR', 'Read1Seq=' + read.query_sequence, 'Read2Seq=' + mateseq]
                            df = pd.concat([df,pd.DataFrame([{'Chromosome':sChr,'Start':sPos-1,'End':sPos,'Var':'INS','Strand':sStrand,'Info':';'.join(info)}])]).reset_index(drop=True)
                        
        elif read.is_proper_pair and read.is_reverse and rightSoftClip >= minSoftClip:
            if read.has_tag('SA'):
                for sa in read.get_tag('SA').rstrip(';').split(';'):
                    sChr, sPos, sStrand, sCigar, sMq, sNm = sa.split(',')
                    if sChr != read.reference_name and int(sMq)>=minMq and int(sNm) <= maxMismatches:
                        readAligned = cigar_to_aligned_positions(read.cigarstring)
                        saAligned = cigar_to_aligned_positions(sCigar)
                        if len(readAligned.intersection(saAligned)) / len(readAligned) < 0.2:
                            if sStrand == '+':
                                sPos = int(sPos) + len(saAligned)

                            mateseq = '.'
                            if read.reference_name != read.next_reference_name:
                                mate = bam.mate(read)
                                mateseq = mate.query_sequence
                            
                            info = ['ID='+read.query_name, 'Type=PR', 'Read1Seq=' + read.query_sequence, 'Read2Seq=' + mateseq]
                            df = pd.concat([df,pd.DataFrame([{'Chromosome':sChr,'Start':int(sPos)-1,'End':int(sPos),'Var':'INS','Strand':sStrand,'Info':';'.join(info)}])]).reset_index(drop=True)
        else:
            continue

    return(pr.PyRanges(df).sort().df)

parser = argparse.ArgumentParser(description='Find split and discordant reads that partially map to a transgene sequence')
parser.add_argument('-n', '--name',type=str,help='Name of transgene contig in reference FASTA')
parser.add_argument('-r','--reference',type=str,default=None,help='Reference FASTA file')
parser.add_argument('-x','--exclude',type=str,default=None,help='Coordinates to exclude from transgene contig')
parser.add_argument('-o','--outfile',type=str,default=None,help='Output to file [stdout]')
parser.add_argument('expbamfile',type=str,help='BAM file')

args = parser.parse_args()

contig = args.name

# open bam file(s)
expsamfile = pysam.AlignmentFile(args.expbamfile,"rc",reference_filename=args.reference)

chimericReads = get_chimeras(expsamfile,contig,exclude=args.exclude)

expsamfile.close()

if args.outfile:
    chimericReads.to_csv(args.outfile,sep="\t",index=False)
else:
    chimericReads.to_csv(sys.stdout,sep="\t",index=False)

