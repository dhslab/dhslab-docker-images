#!/usr/bin/env python3

from __future__ import division
import os
import importlib.util
import pysam
import biotite.sequence.align as align
import biotite.sequence as seq
import argparse
import re
import string
import csv, sys
import scipy.stats as stats
import pandas as pd, pyranges as pr
import joblib
import numpy as np

# ---------------------------------------------------------------------------
# Import shared feature helpers from crispr_ml_features.py
# ---------------------------------------------------------------------------

def _load_crispr_ml_features():
    """Load crispr_ml_features.py from the same bin/ directory."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    feat_path = os.path.join(script_dir, "crispr_ml_features.py")
    if not os.path.exists(feat_path):
        raise FileNotFoundError(
            f"Cannot find crispr_ml_features.py at {feat_path}"
        )
    spec = importlib.util.spec_from_file_location("crispr_ml_features", feat_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

_FEAT = _load_crispr_ml_features()

_softclip_pam_distance        = _FEAT._softclip_pam_distance
_precompute_site_ref_features = _FEAT._precompute_site_ref_features

# ============================================================================
# SECTION 1: ORIGINAL INDEL CALCULATION FUNCTIONS
# ============================================================================

def revcomp(seq):
    tab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh') # maketrans <- maps the reverse complement
    return seq.translate(tab)[::-1] # translate(x)[::-1] <- works backward through the string, effectively reversing the string

def make_cigar_tuples(cigar):
    cigar_dict = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
    # use regex to parse cigar string into tuples
    cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    cigar_tuples = [(cigar_dict[operation],int(length)) for length, operation in cigar_tuples]
    return cigar_tuples

def cigar_summary(cigar):
    cigar_dict = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
    # use regex to parse cigar string into tuples
    cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    cigar_sum = { 'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0 }
    for length, operation in cigar_tuples:
        cigar_sum[operation] += int(length)

    return cigar_sum

def make_info_string(data):
    """
    Converts a single dict or a list of dicts into a VCF INFO string.
    Aggregates values for shared keys across a list.
    Handles pysam types (tuples for lists, booleans for flags).
    """
    if not data:
        return "."

    # Normalize input: ensure we always iterate over a list of dicts
    if isinstance(data, dict):
        data_list = [data]
    elif isinstance(data, list):
        data_list = data
    else:
        return "."

    info_parts = []
    
    # 1. Identify all unique keys present across the input
    all_keys = sorted(set().union(*(d.keys() for d in data_list)))

    for k in all_keys:
        vals = []
        is_flag = False

        for d in data_list:
            if k in d:
                v = d[k]
                # Handle pysam Flags (boolean True implies presence)
                if isinstance(v, bool):
                    if v: 
                        is_flag = True
                # Handle pysam Tuples/Lists (e.g., AF=0.1,0.2)
                elif isinstance(v, (tuple, list)):
                    vals.extend([str(x) for x in v])
                # Handle scalars (int, str, float)
                else:
                    vals.append(str(v))

        # 2. Format string based on content
        # If it was a flag and we collected no value data (pure flag)
        if is_flag and not vals:
            info_parts.append(k)
        # If we have values, join them (e.g. KEY=val1,val2)
        elif vals:
            info_parts.append(f"{k}={','.join(vals)}")

    return ';'.join(info_parts)    

def indels_from_aligned_pairs(pairs,readseq):

    # remove leading soft clips
    while pairs[0][0] is None:
        pairs = pairs[1:]

    # remove training soft clips
    while pairs[-1][0] is None:
        pairs = pairs[:-1]

    pairs = [ (x[0],x[1],x[2].upper(),None) for x in pairs ]
    for i in range(len(pairs)):
        if pairs[i][0] is not None:
            pairs[i] = (pairs[i][0],readseq[pairs[i][0]],pairs[i][1],pairs[i][2])

    variant_start_index = j-1
    variant_end_index = j-1
    j = 0
    while j < len(pairs):
        if pairs[j][0] is not None and pairs[j][2] is not None:
            j+=1
            continue

        if pairs[j][0] is None or pairs[j][2] is None: # indel
            variant_start_index = j-1
            variant_end_index = len(pairs)
            break

        j+=1

    # construct indel tuple
    indel = (pairs[variant_start_index][2],''.join([ x[3] if x[3] else '' for x in pairs[variant_start_index:variant_end_index]]),''.join([ x[1] if x[1] else '' for x in pairs[variant_start_index:variant_end_index]]))
    
    return indel

# function to count the number of reads that support an indel or BND event
def add_normal_counts(df, reads, fasta, maxreads=0,handicap=5,window=25,flank=150):

    # make ref and alt sequences for each indel/BND
    df['refseq'] = ''
    df['altseq'] = ''
    for it, row in df.iterrows():
        df.at[it,'refseq'] = fasta.fetch(row['chrom'],row['pos']-1-flank,row['pos']+len(row['ref'])+flank)
        if row['type'] in ['INDEL','DEL','INS']:
            df.at[it,'altseq'] = fasta.fetch(row['chrom'],row['pos']-1-flank,row['pos']-1) + row['alt'] + fasta.fetch(row['chrom'],row['pos']+len(row['ref'])-1,row['pos']+len(row['ref'])+flank)

        elif row['type'] == 'BND':
            df.at[it,'altseq'] = fasta.fetch(row['chrom'],row['pos']-1-flank,row['pos']-1) + (fasta.fetch(row['chrom2'],row['pos2']-1,row['pos2']+flank) if row['strands']=='++' else revcomp(fasta.fetch(row['chrom2'],row['pos2']-1-flank,row['pos2'])))

    df['control_alt_counts'] = 0
    df['control_total_counts'] = 0

    total_reads = set()

    # iterate through reads in bam file for the first position
    for read in reads:
        # skip if not primary alignment or a duplicate or alignment doesnt overlap start,end
        if read.is_mapped is False or \
            read.is_duplicate is True or \
            read.is_secondary is True or \
            read.is_supplementary is True or \
            read.mapping_quality == 0:
            continue

        total_reads.add(read.query_name)

        # for speed: if the read has no indels, add it to the ref_count set
        if len(read.cigartuples) == 1 and read.cigartuples[0][0] == 0 and read.cigartuples[0][1] == len(read.query_sequence):
            continue
        
        for it, row in df.iterrows():
            if read.cigartuples[0][0]==0 and read.cigartuples[-1][0]==0:
                if (len(row['ref']) > len(row['alt']) and 'D' not in read.cigarstring) or (len(row['alt']) > len(row['ref']) and 'I' not in read.cigarstring):
                    continue

            if row['type'] == 'BND' and not read.has_tag('SA'):
                continue

            refseq = row['refseq']
            altseq = row['altseq']

            ref_align = align.align_optimal(seq.NucleotideSequence(row['refseq']),seq.NucleotideSequence(read.query_sequence),matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),gap_penalty=(-10,-1),local=True,terminal_penalty=False,max_number=1)[0]
        
            alt_align = align.align_optimal(seq.NucleotideSequence(row['altseq']),seq.NucleotideSequence(read.query_sequence),matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),gap_penalty=(-10,-1),local=True,terminal_penalty=False,max_number=1)[0]

            if alt_align.score > 0 and alt_align.score - handicap > ref_align.score:
                ref_align_cigar_sum = cigar_summary(align.write_alignment_to_cigar(ref_align))
                alt_align_cigar_sum = cigar_summary(align.write_alignment_to_cigar(alt_align))

                if ((len(row['ref'])-len(row['alt']) > 0 and ref_align_cigar_sum['D']==len(row['ref'])-len(row['alt'])) or 
                (len(row['alt'])-len(row['ref']) > 0 and ref_align_cigar_sum['I']==len(row['alt'])-len(row['ref']))):
                    df.at[it,'control_alt_counts'] += 1
                    # print to stderr: found control alt count
                    print(f"\tFound control alt count for {row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}:{read.query_name}", file=sys.stderr)

    df['control_total_counts'] = len(total_reads)
    return df.copy()

# function to get indels for one amplicon from bam file
def get_indels(bam,controlbam,chr,start,end,fasta,window=100,distance=25,pam_positions=pd.NA,minSecMapQual=10,maxNM=5,svDistanceThreshold=100000,saAlignmentTolerance=5,minSoftClipLength=5,deletionDistanceThreshold=1000,minreads=1,maxcontrol=0,strict=False,verbose=False):

    # window == flanking sequence to add to read search
    # distance == max distance from start/end to consider an indel, BND, or call a read wild-type.

    cigarVarDict = {'0':'M','1':'I','2':'D','3':'N','4':'S','5':'H','6':'P','7':'=','8':'X'}
    readaln = pd.DataFrame(columns=['read','chrom','pos','chrom2','pos2','ref','alt','strands','type'])
    
    # if verbose, print to stderr: getting reads that align within window of start and end
    if verbose:
        print(f"\tGetting reads that align within window of {chr}:{start}-{end}", file=sys.stderr)

    regionStart = max(start - svDistanceThreshold,1)
    regionEnd = min(end + svDistanceThreshold,fasta.get_reference_length(chr))
    regionSeq = fasta.fetch(chr,regionStart-1,regionEnd)

    # get reads that align within window of start and end
    for read in bam.fetch(chr, start-window, end+window, multiple_iterators = True):

        # skip if not primary alignment or a duplicate or alignment doesnt overlap start,end
        if read.is_mapped is False or \
            read.is_duplicate is True or \
            read.is_secondary is True or \
            read.is_supplementary is True or \
            read.mapping_quality == 0:
            continue

        cigar = read.cigartuples # get cigar info
        mate_cigar = make_cigar_tuples(read.get_tag('MC')) if read.has_tag('MC') else None

        # Quality filter: Get mismatches from NM tag and subtract deleted/skipped bases from cigar
        # and skip if too many mismatches
        if read.has_tag('NM') and int(read.get_tag('NM')) > maxNM:
            if (int(read.get_tag('NM')) - sum([x[1] for x in cigar if x[0] == 2 or x[0]==3])) > maxNM:
                continue

        read_strand = "+"
        if read.is_reverse is True:
            read_strand = "-"

        leftSoftClip = 0
        rightSoftClip = 0
        read_reference_start = read.reference_start + 1
        read_reference_end = read.reference_end
        read_next_reference_start = read.next_reference_start + 1
        read_next_reference_end = None if mate_cigar is None else read_next_reference_start + sum([x[1] for x in mate_cigar if x[0] == 0 or x[0] == 2 or x[0] == 3])
        
        read_query_start = 0
        read_query_end = len(read.query_sequence)

        # setup indelinfo dictionary
        indelinfo = {'read':read.query_name,'chrom':'','pos':None,'chrom2':None,'pos2':None,'strands': '','ref':'','alt':'','type':''}

        # Get left soft clip length
        if cigar[0][0] == 4:
            leftSoftClip = cigar[0][1]

        # Get right soft clip length
        if cigar[-1][0] == 4:
            rightSoftClip = cigar[-1][1]

        if leftSoftClip > rightSoftClip:
            sChr, sPos, sStrand, sCigar, sMq, sNm = read.get_tag('SA').split(';')[0].split(',') if read.has_tag('SA') else [None,None,None,None,None,None]

            if (sChr != None and int(sMq)>=minSecMapQual and
                int(sNm)<maxNM and sChr == read.reference_name and
                sStrand == read_strand and int(sPos) < read_reference_start and
                abs(int(sPos) - read_reference_start) < svDistanceThreshold):
                
                sPos = int(sPos)
                sCigarTuples = make_cigar_tuples(sCigar)

                if sCigarTuples[-1][0] >= 4 and sCigarTuples[0][0] == 0 and sCigarTuples[0][1] == cigar[0][1]:
                    cigar = sCigarTuples[:-1] + [(2,read_reference_start-sPos-1)] + cigar[1:]
                else:
                    cigar = make_cigar_tuples(align.write_alignment_to_cigar(align.align_optimal(seq.NucleotideSequence(regionSeq[sPos-regionStart:read_reference_end-regionStart+1]),
                                                                                seq.NucleotideSequence(read.query_sequence),
                                                                                matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                                                                                gap_penalty=(-10,-1),local=False,max_number=1)[0]))
                read_reference_start = sPos

            else:
                if leftSoftClip >= minSoftClipLength:
                    clippedSeq = read.query_sequence[:leftSoftClip]
                    refSeq = regionSeq[read_reference_start - regionStart - window + 1:read_reference_start - regionStart +  1]
                    leftClipPos = refSeq.rfind(clippedSeq)
                    if leftClipPos != -1:
                        delSeq = refSeq[leftClipPos + len(clippedSeq)-1:]
                        cigar = [(0,len(clippedSeq))] + [(2,len(delSeq))] + cigar[1:]
                        read_reference_start = read_reference_end - sum([x[1] for x in cigar if x[0] == 0 or x[0] == 2 or x[0] == 3]) + 1

        # Process right soft clip, if present
        if rightSoftClip > leftSoftClip:
            sChr, sPos, sStrand, sCigar, sMq, sNm = read.get_tag('SA').split(';')[0].split(',') if read.has_tag('SA') else [None,None,None,None,None,None]

            if (sChr != None):
                sPos = int(sPos)
                sCigarTuples = make_cigar_tuples(sCigar)

                if (int(sMq)>=minSecMapQual and sCigarTuples[0][0] >= 4 and sCigarTuples[-1][0] == 0 and
                    int(sNm)<maxNM and sChr == read.reference_name and 
                    sStrand == read_strand and #int(sPos) > read_reference_end and 
                    abs(int(sPos) - read_reference_end) < svDistanceThreshold):
                
                    if sCigarTuples[0][0] >= 4 and sCigarTuples[-1][0] == 0 and sCigarTuples[-1][1] == cigar[-1][1]:
                        if int(sPos) < read_reference_end:
                            cigar = cigar[:-1] + [(1,read_reference_end-sPos)] + sCigarTuples[1:]
                        else:
                            cigar = cigar[:-1] + [(2,read_reference_end-sPos)] + sCigarTuples[1:]

                    else:
                        cigar = make_cigar_tuples(align.write_alignment_to_cigar(align.align_optimal(seq.NucleotideSequence(regionSeq[read_reference_start-regionStart:sPos-regionStart+sum([x[1] for x in sCigarTuples if x[0] == 0 or x[0] == 2 or x[0] == 3])]),
                                        seq.NucleotideSequence(read.query_sequence),
                                        matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                                        gap_penalty=(-10,-1),local=False,max_number=1)[0]))

                read_reference_end = read_reference_start + sum([x[1] for x in cigar if x[0] == 0 or x[0] == 2 or x[0] == 3])

            else:
                if rightSoftClip >= minSoftClipLength:
                    clippedSeq = read.query_sequence[-rightSoftClip:]
                    refSeq = regionSeq[read_reference_end-regionStart:read_reference_end-regionStart + window]
                    rightClipPos = refSeq.find(clippedSeq)
                    if rightClipPos != -1:
                        delSeq = refSeq[0:rightClipPos]
                        cigar = cigar[:-1] + [(2,len(delSeq))] + [(0,len(clippedSeq))]
                        read_reference_end = read_reference_start + sum([x[1] for x in cigar if x[0] == 0 or x[0] == 2 or x[0] == 3])
        
        # Check if adjusted read is within distance
        if read_reference_start - leftSoftClip >= end + distance or read_reference_end + rightSoftClip <= start - distance:
            continue

        # Remove leading or trailing soft clips
        if cigar[0][0] >= 4:
            read_query_start += cigar[0][1]
            del cigar[0]
        
        if cigar[-1][0] >= 4:
            read_query_end -= cigar[-1][1]
            del cigar[-1]

        # Process indels if multiple cigar operations
        if len(cigar) > 1:
            # Adjust alignment start for first aligned block
            if cigar[0][0] == 0:
                read_reference_start += cigar[0][1]
                read_query_start = read_query_start + cigar[0][1] 
                del cigar[0]

            # Adjust alignment end for last aligned block
            if cigar[-1][0] == 0:
                read_reference_end = read_reference_end - cigar[-1][1]
                read_query_end = read_query_end - cigar[-1][1]
                del cigar[-1]

            # Process single indel event
            if len(cigar) == 1:

                # theres one event, so record the coordinates
                indelinfo['chrom'] = read.reference_name
                indelinfo['pos'] = read_reference_start
                indelinfo['strands'] = '++'

                # if insertion
                if cigar[0][0] == 1:           
                    indelinfo['ref'] = regionSeq[read_reference_start-regionStart-1:read_reference_start-regionStart] #fasta.fetch(chr,read_reference_start-1-1,read_reference_start-1) # 0 based and position before the isertion
                    indelinfo['alt'] = read.query_sequence[read_query_start-1:read_query_start+cigar[0][1]]
                    indelinfo['type'] = 'INDEL'

                # if deletion
                elif cigar[0][0] == 2 or cigar[0][0] == 3:
                    indelinfo['ref'] = regionSeq[read_reference_start-regionStart:read_reference_start-regionStart + cigar[0][1] + 1] #fasta.fetch(chr,read_reference_start-1,read_reference_start+cigar[0][1])
                    # Check if ref is not empty before accessing first character
                    if indelinfo['ref']:
                        indelinfo['alt'] = indelinfo['ref'][0]
                    else:
                        # Handle empty ref case - could skip this indel or use a default value
                        indelinfo['alt'] = 'N'  # or skip this indel entirely
                    indelinfo['type'] = 'INDEL'

            # if its a complex indel
            else:
                indelinfo['chrom'] = read.reference_name
                indelinfo['pos'] = read_reference_start
                indelinfo['strands'] = '++'
                indelinfo['type'] = 'INDEL'

                # iterate through cigar and get indels
                varlen = sum([ x[1] for x in cigar ])
                indelinfo['ref'] = regionSeq[read_reference_start-regionStart:read_reference_start-regionStart+sum([ x[1] for x in cigar ])+1] #fasta.fetch(chr,read_reference_start-1,read_reference_start+sum([ x[1] for x in cigar ]))
                indelinfo['alt'] = read.query_sequence[read_query_start+1:read_query_start+sum([x[1] if x[0] == 1 or x[0] == 0 else 0 for x in cigar])]

        # Process chimeric reads with supplementary alignments
        elif read.has_tag('SA') is True and (leftSoftClip > 0 or rightSoftClip > 0):
            sChr, sPos, sStrand, sCigar, sMq, sNm = read.get_tag('SA').split(';')[0].split(',') if read.has_tag('SA') else [None,None,None,None,None,None]

            sCigarTuples = make_cigar_tuples(sCigar)
            sPos = int(sPos)
            sEnd = sPos + sum([x[1] for x in sCigarTuples if x[0] == 0 or x[0] == 2 or x[0] == 3]) - 1
            sLeftSoftClip = sCigarTuples[0][1] if sCigarTuples[0][0] >= 4 else 0
            sRightSoftClip = sCigarTuples[-1][1] if sCigarTuples[-1][0] >= 4 else 0

            pSeq = read.query_alignment_sequence
            sSeq = read.query_sequence[sLeftSoftClip:sLeftSoftClip+sum([x[1] for x in sCigarTuples if x[0] == 0 or x[0] == 1])] if sStrand == read_strand else revcomp(read.query_sequence[sLeftSoftClip:sLeftSoftClip+sum([x[1] for x in sCigarTuples if x[0] == 0 or x[0] == 1])])

            # Correct for overlapping alignments
            if len(pSeq)+len(sSeq) > len(read.query_sequence):
                trimLen = len(pSeq)+len(sSeq) - len(read.query_sequence)
                if sLeftSoftClip > sRightSoftClip:
                    sPos += trimLen
                    sSeq = sSeq[trimLen:]
                elif sRightSoftClip > sLeftSoftClip:
                    sEnd -= trimLen
                    sSeq = sSeq[:-trimLen]

            if (int(sMq) >= minSecMapQual and int(sNm) <= maxNM and int(sMq)>0 and
                (sChr != read.reference_name or abs(int(sPos) - read_reference_start) >= svDistanceThreshold) or sStrand != read_strand or 
                (rightSoftClip > 0 and sEnd < read_reference_start) or (leftSoftClip > 0 and sPos > read_reference_end)):

                indelinfo['chrom'] = read.reference_name
                indelinfo['chrom2'] = sChr
                indelinfo['strands'] = '++'if sStrand == read_strand else '+-'
                indelinfo['type'] = 'BND'

                if indelinfo['strands'] == '++':
                    if rightSoftClip > leftSoftClip and sLeftSoftClip > sRightSoftClip:
                        indelinfo['pos'] = read_reference_end
                        indelinfo['pos2'] = sPos
                    elif leftSoftClip > rightSoftClip and sRightSoftClip > sLeftSoftClip:
                        indelinfo['pos'] = read_reference_start - 1 
                        indelinfo['pos2'] = sEnd + 1
                    else:
                        print(f"Error: no proper soft clip orientation found: {read.query_name}: {str(read)}", file=sys.stderr)
                        if strict:
                            exit(1)
                        else:
                            continue

                elif indelinfo['strands'] == '+-':
                    if rightSoftClip > leftSoftClip and sRightSoftClip > sLeftSoftClip:
                        indelinfo['pos'] = read_reference_end
                        indelinfo['pos2'] = sEnd
                    elif leftSoftClip > rightSoftClip and sLeftSoftClip > sRightSoftClip:
                        indelinfo['pos'] = read_reference_start - 1
                        indelinfo['pos2'] = sPos
                    else:
                        print(f"Error: no proper soft clip orientation found: {read.query_name}: {str(read)}", file=sys.stderr)
                        if strict:
                            exit(1)
                        else:
                            continue
                else:
                    print(f"Error: no proper soft clip orientation found: {read.query_name}: {str(read)}", file=sys.stderr)
                    if strict:
                        exit(1)
                    else:
                        continue
            
            elif read_reference_start < end and read_reference_end > start:
                indelinfo['chrom'] = chr
                indelinfo['pos'] = start
                indelinfo['ref'] = '.'
                indelinfo['alt'] = '.'
                indelinfo['type'] = 'REF'
            else:
                continue
                    
        elif read_reference_start < end and read_reference_end > start:
            indelinfo['chrom'] = chr
            indelinfo['pos'] = start
            indelinfo['ref'] = '.'
            indelinfo['alt'] = '.'
            indelinfo['type'] = 'REF'
        else:
           continue

        # Add indel info to dataframe
        readaln = pd.concat([readaln, pd.DataFrame([indelinfo])], ignore_index=True)

    if verbose:
        print("\tSorting indels", file=sys.stderr)

    # Group by read and process results
    readaln = readaln.sort_values(by=['read','chrom','pos','chrom2','pos2','strands','ref','alt','type'],key=lambda col: col != '',ascending=False).groupby('read').first().reset_index()
    indelcounts = readaln.groupby(['chrom','pos','chrom2','pos2','strands','ref','alt','type'],dropna=False).size().reset_index(name='counts')

    # Recast as int type, allowing for NA values
    indelcounts['pos'] = indelcounts['pos'].astype(pd.Int64Dtype())
    indelcounts['pos2'] = indelcounts['pos2'].astype(pd.Int64Dtype())

    if verbose:
        print("\tGetting control counts", file=sys.stderr)

    # Get counts from control BAM
    if len(indelcounts) > 0:
        # Add pam positions to the df
        if pam_positions is not pd.NA:
            indelcounts['Positions'] = [pam_positions] * len(indelcounts)
            indelcounts['Distance'] = indelcounts.apply(
                lambda r: min(
                    [abs(r['pos'] - int(x)) for x in r['Positions']] + 
                    [abs(r['pos'] + len(r['ref']) - 1 - int(x)) for x in r['Positions']]
                ) if pd.notna(r['pos']) and r['Positions'] else pd.NA,
                axis=1
            )

        else:
            indelcounts['Positions'] = pd.NA
            indelcounts['Distance'] = pd.NA

        indelcounts = add_normal_counts(indelcounts, [ x for x in controlbam.fetch(contig=chr,start=start-window,end=end) ], fasta, window=distance)
    else:
        indelcounts['control_alt_counts'] = 0
        indelcounts['control_total_counts'] = 0
        indelcounts['Positions'] = pd.NA
        indelcounts['Distance'] = pd.NA

    # Apply filters
    indelcounts = indelcounts[(indelcounts['counts'] >= minreads) | (indelcounts['ref']=='.')]
    indelcounts = indelcounts[(indelcounts['control_alt_counts'] <= maxcontrol) | (indelcounts['ref']=='.')]
    indelcounts = indelcounts[(indelcounts['Distance'] <= distance) | (indelcounts['ref']=='.')]
    indelcounts = indelcounts[~indelcounts['alt'].str.contains('N')]

    return(indelcounts)

# ============================================================================
# SECTION 2: CRISPR PREDICTION FUNCTIONS
# ============================================================================

def detect_deletion_from_softclips(read, min_softclip_size=3, max_gap=50):
    """Detect potential large deletions from soft-clips and supplementary alignments."""
    if not read.cigartuples:
        return 0
    
    total_deletion = 0
    
    # Parse CIGAR operations
    cigar_ops = []
    for op, length in read.cigartuples:
        cigar_ops.append((op, length))
    
    # Count soft-clips as potential deletions
    for op, length in cigar_ops:
        if op == 4:  # Soft-clip
            if length >= min_softclip_size:
                total_deletion += length
    
    # Look for patterns: softclip -> deletion -> softclip
    for i in range(len(cigar_ops) - 2):
        op1, len1 = cigar_ops[i]
        op2, len2 = cigar_ops[i + 1]
        op3, len3 = cigar_ops[i + 2]
        
        if op1 == 4 and op2 == 2 and op3 == 4:  # Softclip -> Deletion -> Softclip
            if len1 >= min_softclip_size and len2 <= max_gap and len3 >= min_softclip_size:
                total_deletion += len2
        elif op1 == 2 and op2 == 4 and op3 == 2:  # Deletion -> Softclip -> Deletion
            if len2 >= min_softclip_size and len1 <= max_gap and len3 <= max_gap:
                total_deletion += (len1 + len3)
    
    # Look for single softclip patterns with nearby deletions
    for i in range(len(cigar_ops) - 1):
        op1, len1 = cigar_ops[i]
        op2, len2 = cigar_ops[i + 1]
        
        if op1 == 4 and op2 == 2:  # Softclip -> Deletion
            if len1 >= min_softclip_size and len2 <= max_gap:
                total_deletion += len2
        elif op1 == 2 and op2 == 4:  # Deletion -> Softclip
            if len2 >= min_softclip_size and len1 <= max_gap:
                total_deletion += len1
    
    # Look for end soft-clips (left/right) representing large deletions
    if len(cigar_ops) >= 1:
        # Left soft-clip
        if cigar_ops[0][0] == 4:
            left_softclip = cigar_ops[0][1]
            if left_softclip >= min_softclip_size:
                middle_align = sum(length for op, length in cigar_ops[1:] if op == 0)
                if middle_align > 5:
                    total_deletion += left_softclip
        
        # Right soft-clip
        if cigar_ops[-1][0] == 4:
            right_softclip = cigar_ops[-1][1]
            if right_softclip >= min_softclip_size:
                middle_align = sum(length for op, length in cigar_ops[:-1] if op == 0)
                if middle_align > 5:
                    total_deletion += right_softclip
    
    # Look for internal soft-clips
    for i, (op, length) in enumerate(cigar_ops):
        if op == 4 and length >= min_softclip_size:
            left_align = sum(cigar_ops[j][1] for j in range(i) if cigar_ops[j][0] == 0)
            right_align = sum(cigar_ops[j][1] for j in range(i+1, len(cigar_ops)) if cigar_ops[j][0] == 0)
            
            if left_align > 2 and right_align > 2:
                total_deletion += length
    
    # Check supplementary alignments (SA tags) for soft-clip deletions
    if read.has_tag('SA'):
        sa_tag = read.get_tag('SA')
        sa_entries = sa_tag.split(';')
        
        for sa_entry in sa_entries:
            if not sa_entry:
                continue
            
            sa_parts = sa_entry.split(',')
            if len(sa_parts) >= 4:
                sa_cigar = sa_parts[3]
                sa_softclips = parse_cigar_for_softclips(sa_cigar)
                for softclip_size in sa_softclips:
                    if softclip_size >= min_softclip_size:
                        total_deletion += softclip_size
    
    return total_deletion

def parse_cigar_for_softclips(cigar_string):
    """Parse CIGAR string to extract soft-clip sizes."""
    softclips = []
    import re
    
    cigar_parts = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)
    
    for length_str, operation in cigar_parts:
        if operation == 'S':  # Soft-clip
            softclips.append(int(length_str))
    
    return softclips

def parse_cigar_once(read):
    """Parse CIGAR operations to extract insertion, deletion, and soft-clip counts."""
    cigar_summary = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
    if not hasattr(read, 'cigartuples') or read.cigartuples is None:
        return {'insertions': 0, 'deletions': 0, 'softclips': 0}
    
    for op, length in read.cigartuples:
        if op == 1:
            cigar_summary['I'] += length
        elif op == 2:
            cigar_summary['D'] += length
        elif op == 4:
            cigar_summary['S'] += length
    
    # Total deletions = explicit deletions + deletions from soft-clips
    total_deletions = cigar_summary['D'] + detect_deletion_from_softclips(read)
    
    return {
        'insertions': cigar_summary['I'],
        'deletions': total_deletions,
        'softclips': cigar_summary['S']
    }

def count_mismatches_fast(read):
    """Count mismatches from MD tag."""
    if not read.has_tag('MD'):
        return 0
    md = read.get_tag('MD')
    return sum(1 for c in md if c.isalpha())

def calculate_control_fractions(control_bam, chrom, position, window=50):
    """Calculate fractions of control reads with different variant types."""
    fractions = {
        'fraction_control_reads_del': 0.0,
        'fraction_control_reads_ins': 0.0,
        'fraction_control_reads_mismatch': 0.0,
        'fraction_control_reads_softclip': 0.0
    }
    if not control_bam:
        return fractions
        
    total = del_count = ins_count = mismatch_count = softclip_count = 0
    
    for read in control_bam.fetch(chrom, position - window, position + window):
        if read.is_unmapped or read.is_duplicate:
            continue
        total += 1
        if read.cigartuples:
            if any(op == 2 for op, _ in read.cigartuples): 
                del_count += 1
            if any(op == 1 for op, _ in read.cigartuples): 
                ins_count += 1
            if any(op == 4 for op, _ in read.cigartuples): 
                softclip_count += 1
        if read.has_tag('MD') and any(c.isalpha() for c in read.get_tag('MD')):
            mismatch_count += 1
    
    if total > 0:
        fractions.update({
            'fraction_control_reads_del': del_count / total,
            'fraction_control_reads_ins': ins_count / total,
            'fraction_control_reads_mismatch': mismatch_count / total,
            'fraction_control_reads_softclip': softclip_count / total
        })
    return fractions

def calculate_exclusivity_features(read, control_bam, chrom, position, window=100):
    """Calculate whether variants are exclusive to edited samples."""
    exclusivity = {
        'variant_exclusive_to_edited': 0,
        'insertion_exclusive_to_edited': 0,
        'deletion_exclusive_to_edited': 0,
        'control_has_same_variant': 0,
    }
    if not control_bam or read.is_unmapped or read.is_duplicate or not read.cigartuples:
        return exclusivity

    has_insertion = any(op == 1 for op, _ in read.cigartuples)
    has_deletion = any(op == 2 for op, _ in read.cigartuples)
    
    read_start, read_end = read.reference_start, read.reference_end

    control_has_same_ins = False
    control_has_same_del = False
    
    for control_read in control_bam.fetch(chrom, position - window, position + window):
        if control_read.is_unmapped or control_read.is_duplicate or not control_read.cigartuples:
            continue
        ctrl_start, ctrl_end = control_read.reference_start, control_read.reference_end
        
        # Check for overlapping reads
        if not (read_end < ctrl_start or ctrl_end < read_start):
            if has_insertion and any(op == 1 for op, _ in control_read.cigartuples):
                control_has_same_ins = True
            if has_deletion and any(op == 2 for op, _ in control_read.cigartuples):
                control_has_same_del = True

    insertion_exclusive = has_insertion and not control_has_same_ins
    deletion_exclusive = has_deletion and not control_has_same_del
    
    any_exclusive = insertion_exclusive or deletion_exclusive
    control_has_similar = control_has_same_ins or control_has_same_del

    exclusivity['deletion_exclusive_to_edited'] = 1 if deletion_exclusive else 0
    exclusivity['control_has_same_variant'] = 1 if control_has_similar else 0
    return exclusivity

def calculate_distance_to_closest_pam(position, targets_df, chrom):
    """Calculate distance to closest PAM site."""
    if targets_df is None or len(targets_df) == 0:
        return -1
    distances = np.abs(targets_df['Start'] - position)
    return distances.min() if len(distances) > 0 else -1

def predict_reads_at_position(bam_file, chrom, start, end, pampos, model, fasta, is_on_target=0, control_bam=None, threshold=0.80):

    reads = []
    for read in bam_file.fetch(chrom, start, end):
        if read.is_unmapped or read.is_duplicate:
            continue
        reads.append(read)
    if not reads:
        return 0, 0.0

    # Pre-compute site-level reference features once (fasta is a pysam.FastaFile handle)
    ref_site_feats = _precompute_site_ref_features(fasta, chrom, start)

    # Pre-compute site-level control features once
    control_fractions = calculate_control_fractions(control_bam, chrom, start, window=50) if control_bam else {
        'fraction_control_reads_del': 0.0,
        'fraction_control_reads_ins': 0.0,
        'fraction_control_reads_mismatch': 0.0,
        'fraction_control_reads_softclip': 0.0
    }

    # Extract features for each read
    features_list = []
    for read in reads:
        cigar_data = parse_cigar_once(read)

        read_start = read.reference_start
        read_end = read.reference_end

        # Check if read overlaps with target sites
        is_at_any_target = 0
        if pampos is not None:
            for p in pampos:
                if read_start <= p + 25 and read_end >= p - 25:
                    is_at_any_target = 1

        # Calculate exclusivity features
        exclusivity = calculate_exclusivity_features(read, control_bam, chrom, start, window=100) if control_bam else {
            'deletion_exclusive_to_edited': 0,
            'control_has_same_variant': 0,
        }

        # Convert read features to binary (still needed for vs-control deltas)
        read_has_deletion = 1 if cigar_data['deletions'] > 0 else 0
        read_has_insertion = 1 if cigar_data['insertions'] > 0 else 0
        read_has_mismatch = 1 if count_mismatches_fast(read) > 0 else 0

        # Calculate indel characteristics
        total_indel_size = cigar_data['insertions'] + cigar_data['deletions']

        # Indel size category
        indel_size_category = 0
        if total_indel_size > 0:
            if total_indel_size <= 3:
                indel_size_category = 1
            elif total_indel_size <= 10:
                indel_size_category = 2
            else:
                indel_size_category = 3

        # Indel complexity score
        indel_complexity_score = 0.0
        if read.cigartuples:
            indel_operations = [op for op, length in read.cigartuples if op in [1, 2]]
            complexity = len(indel_operations) + (total_indel_size / 10.0)
            indel_complexity_score = min(complexity, 10.0)

        # Create features dictionary (14 model features + is_on_target_site metadata)
        features = {
            'read_pair_gap': abs(read.template_length) if hasattr(read, 'template_length') and read.is_paired and read.is_proper_pair else -1,
            'read_softclip': cigar_data['softclips'],
            'read_del_vs_control': read_has_deletion - control_fractions['fraction_control_reads_del'],
            'read_ins_vs_control': read_has_insertion - control_fractions['fraction_control_reads_ins'],
            'read_mismatch_vs_control': read_has_mismatch - control_fractions['fraction_control_reads_mismatch'],
            'deletion_exclusive_to_edited': exclusivity['deletion_exclusive_to_edited'],
            'distance_to_closest_pam': min(abs(x - start) for x in pampos) if pampos is not None else -1,
            'is_on_target_site': is_on_target,
            'is_at_any_target_site': is_at_any_target,
            'total_indel_size': total_indel_size,
            'indel_size_category': indel_size_category,
            'indel_complexity_score': indel_complexity_score,
            'softclip_dist_to_PAM': _softclip_pam_distance(read, pampos if pampos is not None else []),
            'is_in_homopolymer': ref_site_feats.get('is_in_homopolymer', 0),
            'dist_to_microsatellite': ref_site_feats.get('dist_to_microsatellite', 999),
        }

        features_list.append(features)

    if not features_list:
        return 0, 0.0

    # Get expected features for the model
    if hasattr(model, 'feature_names_in_'):
        expected_features = list(model.feature_names_in_)
    else:
        expected_features = [
            'read_pair_gap', 'read_softclip',
            'read_del_vs_control', 'read_ins_vs_control', 'read_mismatch_vs_control',
            'deletion_exclusive_to_edited',
            'is_on_target_site', 'is_at_any_target_site', 'distance_to_closest_pam',
            'total_indel_size', 'indel_size_category', 'indel_complexity_score',
            'softclip_dist_to_PAM', 'is_in_homopolymer', 'dist_to_microsatellite',
        ]

    features_df = pd.DataFrame(features_list)
    features_df = features_df.reindex(columns=expected_features, fill_value=0)

    # Get model predictions
    preds = model.predict_proba(features_df)[:, 1]

    # Calculate results
    avg_probability = float(preds.mean() * 100)

    return int((preds >= threshold).sum()), avg_probability

# ============================================================================
# SECTION 3: MAIN FUNCTION
# ============================================================================

def main():
    
    parser = argparse.ArgumentParser(description='Extract variant reads and predict CRISPR reads')
    parser.add_argument('-f','--fasta',type=str,default="/storage2/fs1/dspencer/Active/clinseq/projects/scge/data/refdata/singh_v4.3.6/hg38_PLVM_CD19_CARv4_cd34.fa",help='Reference fasta file')
    parser.add_argument('-w','--target-window',type=int,default=100,help='Window size')
    parser.add_argument('-d','--distance',type=int,default=25,help='Distance')
    parser.add_argument('-m','--minreads',type=int,default=1,help='Minimum reads')
    parser.add_argument('-c','--chromosome',type=str,default=None,help='Chromosome to process')
    parser.add_argument('-x','--maxcontrol',type=int,default=0,help='Maximum supporting reads to report an indel/bnd event.')
    parser.add_argument('-s','--strict',action='store_true',help='Exit if a read cannot be properly parsed.')
    parser.add_argument('-v','--verbose',action='store_true',help='Print verbose output')
    parser.add_argument('-o','--outfile',type=str,help='Output file (optional)')
    
    # CRISPR prediction arguments
    parser.add_argument('--crispr-model',type=str,default='models/site14_site5_combined_model.pkl',help='Trained CRISPR ML model file (.pkl) for read prediction (default: site14_site5_combined_model.pkl)')
    parser.add_argument('--crispr-threshold',type=float,default=0.70,help='Probability threshold for CRISPR prediction (default: 0.70)')
    parser.add_argument('--enable-crispr-prediction',action='store_true',help='Enable CRISPR read prediction (uses default model and target-file as targets)')
    parser.add_argument('--targets-csv',type=str,help='Alternative target sites CSV file for feature extraction (optional - uses target-file if not specified)')
    parser.add_argument('--filter-off-target-fp', action='store_true', help='Filter off-target sites that are likely false positives (indel reads > 0 but no predicted CRISPR reads).')
    parser.add_argument('--fp-log', type=str, help='Log file for filtered false positive off-target sites.')
    
    # Required BAM files and target file
    parser.add_argument('--edited-bam',type=str,required=True,help='Edited/experimental BAM/CRAM file')
    parser.add_argument('--control-bam',type=str,required=True,help='Control BAM/CRAM file')
    parser.add_argument('--target-file',type=str,required=True,help='Target sites coordinate file (CSV/BED format) - used for both bed coordinates and CRISPR targets')

    args = parser.parse_args()

    if args.verbose:
        print("Processing input file", file=sys.stderr)

    # ========================================================================
    # STEP 1: Process input BED file and create genomic intervals
    # ========================================================================
    
    # target regions, in VCF format.
    vcf_data = []

    try:
        # Open VCF file with pysam
        vcf_in = pysam.VariantFile(args.target_file)
        
        for rec in vcf_in:
            start = rec.pos - 1
            end = rec.pos
            info_dict = dict(rec.info)
            is_target = 1 if info_dict.get('TARGET', False) else 0
            
            vcf_data.append({
                'Chromosome': rec.chrom,
                'Start': start,
                'End': end,
                'Pos': rec.pos,  # Keep 1-based POS for reference/output
                'Info': info_dict,
                'Ontarget': is_target
            })

    except Exception as e:
        print(f"Error reading VCF input: {e}", file=sys.stderr)

    # Create DataFrame from list of dicts
    bedDf = pd.DataFrame(vcf_data)

    # Check for empty input data
    if len(bedDf) == 0:
        print("Warning: Input target file has no data rows. Writing empty output.", file=sys.stderr)
        # Write empty output file with header only
        output_columns = ['Cluster', 'Chromosome', 'Start', 'End', 'Ontarget', 'Gene', 'indel_type', 
                        'indel_fraction', 'indel_allele_fraction', 'indel_size', 'indel_bases', 
                        'num_edited', 'num_control', 'total_edited', 'total_control', 'significance',
                        'prediction', 'probability', 'model_info']
        empty_df = pd.DataFrame(columns=output_columns)
        if args.outfile:
            empty_df.to_csv(args.outfile, sep='\t', index=False)
        else:
            empty_df.to_csv(sys.stdout, sep='\t', index=False)
        sys.exit(0)

    # Create PyRanges object and cluster intervals
    bedPr = pr.PyRanges(bedDf[['Chromosome', 'Start', 'End', 'Pos', 'Info', 'Ontarget']])
    bedPr = bedPr.cluster(slack=args.target_window)
    mergedBedPr = bedPr.merge(by='Cluster', strand=False, slack=args.target_window)

    # Join aggregated data back to the merged intervals
    mergedBedDf = mergedBedPr.df.join(
        bedPr.df.groupby('Cluster')['Pos'].agg(list).reset_index().set_index('Cluster'),
        on='Cluster', how='left'
    )
    mergedBedDf = mergedBedDf.join(
        bedPr.df.groupby('Cluster')['Info'].agg(list).reset_index().set_index('Cluster'),
        on='Cluster', how='left'
    )
    mergedBedDf = mergedBedDf.join(
        bedPr.df.groupby('Cluster')['Ontarget'].agg('max').reset_index().set_index('Cluster'),
        on='Cluster', how='left'
    )

    # Filter by chromosome if specified
    if args.chromosome is not None:
        print("Processing chromosome", args.chromosome, file=sys.stderr)
        if args.chromosome.startswith('chr'):
            mergedBedDf = mergedBedDf[mergedBedDf['Chromosome'] == args.chromosome]
        else:
            mergedBedDf = mergedBedDf[(mergedBedDf['Chromosome'] == args.chromosome) | 
                                    (mergedBedDf['Chromosome'] == f"chr{args.chromosome}")]
        
        if len(mergedBedDf) == 0:
            print(f"Warning: No data found for chromosome {args.chromosome}", file=sys.stderr)
            all_chromosomes = sorted(bedPr.df['Chromosome'].unique())
            print(f"Available chromosomes: {all_chromosomes}", file=sys.stderr)
            sys.exit(1)
        else:
            print(f"Found {len(mergedBedDf)} intervals for chromosome {args.chromosome}", file=sys.stderr)

    if args.verbose:
        print("Done processing input file", file=sys.stderr)

    # ========================================================================
    # STEP 2: Open BAM files and reference
    # ========================================================================
    
    expsamfile = pysam.AlignmentFile(args.edited_bam,"rc",reference_filename=args.fasta)
    consamfile = pysam.AlignmentFile(args.control_bam,"rc",reference_filename=args.fasta)
    refFasta = pysam.FastaFile(args.fasta)

    # Redirect output if specified
    if args.outfile:
        sys.stdout = open(args.outfile, 'w')

    # ========================================================================
    # STEP 3: Load CRISPR prediction model and targets (if enabled)
    # ========================================================================    

    crispr_model = None
    if args.enable_crispr_prediction:
        try:
            print(f"Loading CRISPR ML model: {args.crispr_model}", file=sys.stderr)
            crispr_model = joblib.load(args.crispr_model)
            print(f"CRISPR model loaded successfully", file=sys.stderr)
            print(f"CRISPR prediction enabled with threshold {args.crispr_threshold}", file=sys.stderr)
        except Exception as e:
            print(f"Error loading CRISPR model '{args.crispr_model}': {e}", file=sys.stderr)
            print(f"Make sure the model file exists in the current directory or provide full path", file=sys.stderr)
            sys.exit(1)

    # ========================================================================
    # STEP 4: Create output header
    # ========================================================================
    
    header_columns = ['chrom', 'start', 'end', 'pam_positions', 'total_reads', 'indel_reads', 'indel_fraction', 
                     'control_reads', 'control_indel_reads', 'control_indel_fraction', 'indel_count', 'indel_info', 
                     'bnd_count', 'bnd_info', 'target_info', 'is_target']
    if args.enable_crispr_prediction:
        header_columns.extend(['crispr_predicted_reads', 'crispr_prediction_fraction', 'crispr_prediction_probability'])
    print("\t".join(header_columns), flush=True)

    fp_log = None
    if args.fp_log:
        fp_log = open(args.fp_log, 'w')
        fp_log.write("\t".join(header_columns) + "\n")

    # ========================================================================
    # STEP 5: Process each genomic interval
    # ========================================================================
    
    for index, row in mergedBedDf.iterrows():
        if args.verbose:
            print(f"Processing interval {row['Chromosome']}:{row['Start']}-{row['End']}", file=sys.stderr)

        # ORIGINAL INDEL CALCULATION
        indels = get_indels(bam=expsamfile,controlbam=consamfile,chr=row['Chromosome'],start=row['Start'],end=row['End'],
                            fasta=refFasta,window=args.target_window,distance=args.distance,pam_positions=row['Pos'],minreads=args.minreads,maxcontrol=args.maxcontrol,verbose=args.verbose)
        
        # Process indel results
        total_reads, indel_reads, control_total_reads, control_indel_reads = 0, 0, 0, 0
        indel_fraction, control_indel_fraction = 0, 0
        indel_keys, bnd_keys = '.', '.'
        bnds = []

        if len(indels) > 0:
            total_reads = sum(indels['counts'])
            indel_reads = sum(indels[indels['type']!='REF']['counts'])
            control_indel_reads = int(indels['control_alt_counts'].mean())
            control_total_reads = int(indels['control_total_counts'].mean())

            # Separate BNDs and indels
            bnds = indels[indels['type']=='BND'].copy()
            indels = indels[indels['type']=='INDEL'].copy()

            if len(indels) > 0:
                indels['Key'] = indels.apply(lambda r: f"{r['chrom']}:{r['pos']}:{r['ref']}:{r['alt']}:{r['counts']}:{r['control_alt_counts']}:{r['Distance']}", axis=1)

            if len(bnds) > 0:
                bnds['Key'] = bnds.apply(lambda r: f"{r['chrom']}:{r['pos']}:{r['chrom2']}:{r['pos2']}:{r['strands']}:{r['counts']}:{r['control_alt_counts']}:{r['Distance']}", axis=1)

        # Calculate fractions and prepare output
        indel_fraction = round(indel_reads/total_reads,4) if total_reads > 0 else 0
        control_indel_fraction = round(control_indel_reads/control_total_reads,4) if control_total_reads > 0 else 0
        indel_keys = ';'.join(indels['Key'].tolist()) if len(indels) > 0 else '.'
        bnd_keys = ';'.join(bnds['Key'].tolist()) if len(bnds) > 0 else '.'
        ontarget = row['Ontarget']
        positions = row['Pos']
        offtargetsites = make_info_string(row['Info']) if row['Info'] else '.'

        # CRISPR PREDICTION (if enabled)
        crispr_predicted_reads = 0
        crispr_prediction_fraction = 0.0
        crispr_prediction_probability = 0.0
        if args.enable_crispr_prediction and crispr_model:
            if args.verbose:
                print(f"  Predicting CRISPR reads for {row['Chromosome']}:{row['Start']}-{row['End']} (total_reads: {total_reads})", file=sys.stderr)
            
            crispr_predicted_reads, crispr_prediction_probability = predict_reads_at_position(
                expsamfile, row['Chromosome'], row['Start'], row['End'], positions, 
                crispr_model, refFasta, is_on_target=ontarget, control_bam=consamfile, threshold=args.crispr_threshold
            )
            crispr_prediction_fraction = round(crispr_predicted_reads/total_reads, 4) if total_reads > 0 else 0.0
            
            if args.verbose:
                print(f"    Predicted {crispr_predicted_reads}/{total_reads} reads as CRISPR-related (≥{args.crispr_threshold})", file=sys.stderr)
                print(f"    Average prediction probability: {crispr_prediction_probability:.1f}%", file=sys.stderr)

        # Filter off-target false positives if enabled
        if args.filter_off_target_fp:
            if ontarget == 0 and indel_reads > 0 and crispr_predicted_reads == 0:
                if args.verbose:
                    print(f"  Filtering FP off-target site {row['Chromosome']}:{row['Start']}-{row['End']} (indel_reads: {indel_reads}, predicted_reads: {crispr_predicted_reads})", file=sys.stderr)
                
                if fp_log:
                    log_fields = [
                        row['Chromosome'], row['Start'], row['End'], 
                        ';'.join([str(x) for x in row['Pos']]), total_reads, indel_reads, indel_fraction,
                        control_total_reads, control_indel_reads, control_indel_fraction,
                        len(indels), indel_keys, len(bnds), bnd_keys, offtargetsites, ontarget
                    ]
                    if args.enable_crispr_prediction:
                        log_fields.extend([crispr_predicted_reads, crispr_prediction_fraction, round(crispr_prediction_probability, 1)])
                    fp_log.write("\t".join([str(field) for field in log_fields]) + "\n")
                continue
        
        # Output results
        output_fields = [
            row['Chromosome'], row['Start'], row['End'], 
            ';'.join([str(x) for x in row['Pos']]), total_reads, indel_reads, indel_fraction,
            control_total_reads, control_indel_reads, control_indel_fraction,
            len(indels), indel_keys, len(bnds), bnd_keys, offtargetsites, ontarget
        ]
        if args.enable_crispr_prediction:
            output_fields.extend([crispr_predicted_reads, crispr_prediction_fraction, round(crispr_prediction_probability, 1)])
        print("\t".join([str(field) for field in output_fields]), flush=True)

    # ========================================================================
    # STEP 6: Cleanup
    # ========================================================================
    
    if args.outfile:
        sys.stdout.close()

    if fp_log:
        fp_log.close()

    expsamfile.close()
    consamfile.close()
    refFasta.close()

if __name__ == "__main__":
    main()
