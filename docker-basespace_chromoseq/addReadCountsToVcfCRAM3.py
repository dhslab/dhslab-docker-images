from __future__ import division
import pysam
from Bio import pairwise2
import argparse
#from string import maketrans

def revcomp(seq):
    return seq.translate(maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]

MinReads = 3
MinVAF = 0.05
PrintAllVariants = False

parser = argparse.ArgumentParser(description='Add amplicon data to VCF file with variants identified from HaloplexHS data.')
parser.add_argument('vcffile',help='VCF file')
parser.add_argument('cramfile',help='CRAM file')
parser.add_argument('samplename',help='Sample name')
parser.add_argument('-r',"--reference",help='reference genome fasta file')
parser.add_argument('-n',"--minreads",type=int,help='minimum variant supporting reads')
parser.add_argument('-v',"--minvaf",type=float,help='minimum vaf')
parser.add_argument('-f',"--includefiltered",action='store_true',help='Excluded filtered variants [default is false]')
parser.add_argument('-w',"--window",type=int,default=150,
                    help='window for creating ref and alt sequences')

args = parser.parse_args()

mysample = args.samplename

if args.minreads:
    MinReads = args.minreads
    
if args.minvaf:
    MinVAF = args.minvaf

if args.includefiltered is True:
    PrintAllVariants = True

# open vcf file
vcffile = pysam.VariantFile(args.vcffile)
#vcffile2 = pysam.VariantFile(args.vcffile)
# open bam file
samfile = pysam.AlignmentFile(args.cramfile,"rc",reference_filename=args.reference)
# open reffasta
fa = pysam.FastaFile(args.reference)

# window on either side of a variant for indel annotation 
window = args.window

vcffile.header.filters.add("LowReads",None,None,'Fails requirement of having >='+str(MinReads)+' reads supporting the variant')
vcffile.header.filters.add("LowVAF",None,None,'Fails requirement of having >='+str(MinVAF)+' VAF')

vcffile.header.formats.add("NR", 1, 'Integer', 'Number of reference allele reads')
vcffile.header.formats.add("NV", 1, 'Integer', 'Number of variant allele reads')
vcffile.header.formats.add("VAF", 1, 'Float', 'VAF, Variant Allele Fraction, which is the NV field divided by the NR field')

hdr = str(vcffile.header).rstrip().split("\n")
hd = hdr.pop()
print("\n".join(hdr) + '\n' + "\t".join((hd.split("\t"))[0:9]) + '\t' + mysample)

for rec in vcffile.fetch(reopen=True):

    cnts = {'ref':0,'alt':0}

    # if the variant is a substitution
    if len(rec.ref) == len(rec.alts[0]) and len(rec.ref) == 1:

        for pileup in samfile.pileup(rec.contig, rec.pos-1, rec.pos, multiple_iterators = False):

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


        # now get reference for other records within window, if any
        overlappingAlts = []
        for ovrec in vcffile.fetch(contig=rec.contig,start=rec.pos-20,end=rec.pos+20):
            if ovrec.alts[0] != rec.alts[0]:
                ovalt = fa.fetch(ovrec.contig,ovrec.pos-window-1,ovrec.pos-1)
                # if deletion
                if len(ovrec.ref) > len(ovrec.alts[0]):
                    ovalt = ovalt + fa.fetch(ovrec.contig,ovrec.pos+len(ovrec.ref)-1,ovrec.pos+window+len(ovrec.ref))
                    
                elif len(ovrec.ref) < len(ovrec.alts[0]): # if insertion
                    ovalt = ovalt + ovrec.alts[0] + fa.fetch(ovrec.contig,ovrec.pos,ovrec.pos+window)

                else: # if substitution
                    ovalt = fa.fetch(ovrec.contig,ovrec.pos-window-1,ovrec.pos-2) + ovrec.alts[0] + fa.fetch(ovrec.contig,ovrec.pos,ovrec.pos+window)
                    
                overlappingAlts.append(ovalt)
                                                                        
        pu = samfile.pileup(rec.contig, rec.pos-1, rec.pos, multiple_iterators = False)
        for pileup in pu:
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
                            maxscore = 0
                            for s in overlappingAlts:
                                sa = pairwise2.align.localms(rdseq,s, 2, -1, -2, -1,score_only=1,one_alignment_only=1)
                                if maxscore < sa:
                                    maxscore = sa

                            if alnalt > maxscore:                            
                                cnts['alt'] += 1
                            else:
                                cnts['ref'] += 1                                

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

    ispass = False

    if cnts['alt'] < MinReads:
        rec.filter.add("LowReads")

    elif myvaf < MinVAF:
        rec.filter.add("LowVAF")

    else:
        rec.filter.add("PASS")
        ispass = True
        
    if PrintAllVariants is True or ispass is True:
        print ("\t".join(str(rec).rstrip().split("\t")[0:10]))

# end vcf fetch

fa.close()
vcffile.close()
samfile.close()
