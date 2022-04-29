from __future__ import division
import pysam
from Bio import pairwise2
import argparse
#from string import maketrans

def revcomp(seq):
    return seq.translate(maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]

parser = argparse.ArgumentParser(description='Add amplicon data to VCF file with variants identified from HaloplexHS data.')
parser.add_argument('vcffile',help='VCF file')
parser.add_argument('cramfile',help='CRAM file')
parser.add_argument('samplename',help='Sample name')
parser.add_argument('-r',"--reference",help='reference genome fasta file')
parser.add_argument('-w',"--window",type=int,default=150,
                    help='window for creating ref and alt sequences')

args = parser.parse_args()

mysample = args.samplename

# open vcf file
vcffile = pysam.VariantFile(args.vcffile)
# open bam file
samfile = pysam.AlignmentFile(args.cramfile,"rc")
# open reffasta
fa = pysam.FastaFile(args.reference)

# window on either side of a variant for indel annotation 
window = args.window

vcffile.header.formats.add("NR", 1, 'Integer', 'Number of reference allele reads')
vcffile.header.formats.add("NV", 1, 'Integer', 'Number of variant allele reads')
vcffile.header.formats.add("VAF", 1, 'Float', 'VAF, Variant Allele Fraction, which is the NV field divided by the NR field')

hdr = str(vcffile.header).rstrip().split("\n")
hd = hdr.pop()
print("\n".join(hdr) + '\n' + "\t".join((hd.split("\t"))[0:9]) + '\t' + mysample)

for rec in vcffile.fetch():

    cnts = {'ref':0,'alt':0}

    # if the variant is a substitution
    if len(rec.ref) == len(rec.alts[0]) and len(rec.ref) == 1:
        
        for pileup in samfile.pileup(rec.contig, rec.pos-1, rec.pos):            
            if pileup.pos == rec.pos-1: # only consider the variant position

                for read in pileup.pileups:

                    # skip positions with indels
                    if not read.is_del and not read.is_refskip:

                        # count alleles by amplicon
                        if read.alignment.query_sequence[read.query_position] == rec.ref:
                            cnts['ref'] += 1
                        elif read.alignment.query_sequence[read.query_position] == rec.alts[0]:
                            cnts['alt'] += 1
                        
    else: # if indel

        refseqstart = rec.pos-window-1
        refseqend = rec.pos+window
        refseq = fa.fetch(rec.contig,refseqstart,refseqend)
        altseq = fa.fetch(rec.contig,refseqstart,rec.pos-1) # 0-based

        # if deletion
        if len(rec.ref) > len(rec.alts[0]):
            altseq = altseq + fa.fetch(rec.contig,rec.pos+len(rec.ref)-1,rec.pos+window+len(rec.ref))

        else: # if insertion
            altseq = altseq + rec.alts[0] + fa.fetch(rec.contig,rec.pos,rec.pos+window)

        for pileup in samfile.pileup(rec.contig, rec.pos-1, rec.pos):
            if pileup.pos == rec.pos-1:
                for read in pileup.pileups:

                    # if the read has fewer mismatches (NM tag) than the indel length and has no softclipped bases
                    # then assign to reference allele
                    if len(read.alignment.cigartuples) == 1 and read.alignment.cigartuples[0][0] == 0 and read.alignment.cigartuples[0][1] == read.alignment.query_length and read.alignment.get_cigar_stats()[0][10] < abs(len(rec.ref) - len(rec.alts[0])):
                        cnts['ref'] += 1

                    # if the record is a deletion and the read has a deletion at this position
                    # with the proper length then assign to the alt allele
                    elif len(rec.ref) > len(rec.alts[0]) and read.indel < 0 and abs(read.indel) == len(rec.ref)-len(rec.alts[0]):
                        cnts['alt'] += 1
                        
                    # if the record is an insertion and the read has an insertion
                    # with the proper bases and length then assign to the alt allele
                    elif len(rec.ref) < len(rec.alts[0]) and read.indel > 0 and rec.alts[0] == read.alignment.seq[read.query_position:read.query_position+read.indel+1]:
                        cnts['alt'] += 1
                    # if none of the above are satisified, use the alignment method
                    # to assign reads to alleles
                    else:
                        rdseq = read.alignment.query_sequence
                            
                        alnref = pairwise2.align.localms(rdseq,refseq, 2, -1, -2, -1,score_only=1,one_alignment_only=1)
                        alnalt = pairwise2.align.localms(rdseq,altseq, 2, -1, -2, -1,score_only=1,one_alignment_only=1)
                    
                        if alnref >= alnalt:
                            cnts['ref'] += 1
                        else:
                            cnts['alt'] += 1

    myvaf = ((cnts['alt']) / (cnts['ref'] + cnts['alt']))

    mysample=0
    mygt = (0,1)
    for s in rec.samples:
        if rec.samples[s]['GT'] == (1,1) and myvaf > .99:
            mygt = (1,1)
 
    rec.samples[mysample]['GT'] = mygt
    rec.samples[mysample]['NR'] = cnts['alt'] + cnts['ref']
    rec.samples[mysample]['NV'] = cnts['alt']
    rec.samples[mysample]['VAF'] = myvaf

    if cnts['alt'] > 0:
        print ("\t".join(str(rec).rstrip().split("\t")[0:10]))


# end vcf fetch

fa.close()
vcffile.close()
samfile.close()
