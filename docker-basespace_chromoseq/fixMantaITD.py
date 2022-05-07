from __future__ import division
import pysam
import argparse

parser = argparse.ArgumentParser(description='Add amplicon data to VCF file with variants identified from HaloplexHS data.')
parser.add_argument('vcffile',help='VCF file')
parser.add_argument('-r',"--reference",help='reference genome fasta file')

args = parser.parse_args()

# open vcf file
vcffile = pysam.VariantFile(args.vcffile)

# open reffasta
fa = pysam.FastaFile(args.reference)

vcffile.header.formats.add("GT", 1, 'String', 'Sample genotype')

hdr = str(vcffile.header).rstrip().split("\n")
hd = hdr.pop()
print("\n".join(hdr) + '\n' + "\t".join((hd.split("\t"))[0:9]) + '\t' + vcffile.header.samples[0])

for rec in vcffile.fetch(reopen=True):

    if rec.alts[0] == "<DUP>" or rec.alts[0] == "<DUP:TANDEM>":
        rec.alts = (str(rec.ref + fa.fetch(rec.contig,rec.pos,rec.stop) + fa.fetch(rec.contig,rec.pos,rec.stop)),)
        rec.info['SVTYPE']  = 'INS'
        rec.samples[0]['GT'] = (0,1)

    elif rec.info['SVTYPE'] == 'INS':
        rec.samples[0]['GT'] = (0,1)        

    else:
        continue
    
    print ("\t".join(str(rec).rstrip().split("\t")[0:10]))

# end vcf fetch

fa.close()
vcffile.close()

