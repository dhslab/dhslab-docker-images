#!/usr/bin/env python

import sys
from time import gmtime, strftime
import os
import re
from cyvcf2 import VCF
import tempfile
import csv
import binascii

def parse_csq_header(vcf_file):
    for header in vcf_file.header_iter():
        info = header.info(extra=True)
        if b'ID' in info.keys() and info[b'ID'] == b'CSQ':
            format_pattern = re.compile('Format: (.*)"')
            match = format_pattern.search(info[b'Description'].decode())
            return match.group(1).split('|')

def parse_csq_entries(csq_entries, csq_fields):
    transcripts = {}
    for entry in csq_entries:
        values = entry.split('|')
        transcript = {}
        for key, value in zip(csq_fields, values):
            transcript[key] = value
        if transcript['Allele'] not in transcripts.keys():
            transcripts[transcript['Allele']] = []
        transcripts[transcript['Allele']].append(transcript)
    return transcripts

def get_csq_entries_bygene(csq_entries):
    genes = {}
    for entry in csq_entries:
        genes[entry['SYMBOL']] = entry
        
    return genes

def resolve_alleles(entry, csq_alleles):
    alleles = {}
    if entry.is_indel:
        for alt in entry.ALT:
            alt = str(alt)
            if alt[0:1] != entry.REF[0:1]:
                csq_allele = alt
            elif alt[1:] == "":
                csq_allele = '-'
            else:
                csq_allele = alt[1:]
            alleles[alt] = csq_allele
    elif entry.is_sv:
        for alt in alts:
            if len(alt) > len(entry.REF) and 'insertion' in csq_alleles:
                alleles[alt] = 'insertion'
            elif len(alt) < len(entry.REF) and 'deletion' in csq_alleles:
                alleles[alt] = 'deletion'
            elif len(csq_alleles) == 1:
                alleles[alt] = list(csq_alleles)[0]
    else:
        for alt in entry.ALT:
            alt = str(alt)
            alleles[alt] = alt
    return alleles

def get_genes(transcripts):
    genes = set()
    for transcript in transcripts:
        genes.add(transcript['SYMBOL'])

    return genes

def decode_hex(string):
    hex_string = string.group(0).replace('%', '')
    return binascii.unhexlify(hex_string).decode('utf-8')

def convert_aa(codon):
    three = ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Ter"]
    one  = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*"]
    
    for i in range(0,len(three)):
        p = re.compile(three[i])
        codon = p.sub(one[i],codon)

    return codon


#
# Script
#

mincnbins = 10

(script, Name, genevcf_file, svvcf_file, genelist) = sys.argv

# variants to print out
vars = {}
vars['knownsv'] = []
vars['cna'] = []
vars['knowngenes'] = []
vars['genelevel'] = []
vars['novelsv'] = []

# read genelist
knowngenelist = set()

with open(genelist, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        l = row[0]
        p = re.compile("_exon_\d+")
        l = p.sub("",l)
        knowngenelist.add(l)

# get gene-level variants

genevcf = VCF(genevcf_file)
genecsq_fields = parse_csq_header(genevcf)

for variant in genevcf:

    csq = variant.INFO.get('CSQ')
    
    if csq is None:
        sys.exit("No VEP fields")
                    
    transcripts = list(parse_csq_entries(csq.split(','), genecsq_fields).items())[0][1] # just get the first allele in the list.
    genes = get_csq_entries_bygene(transcripts)

    vartype = ''
    if len(variant.REF) == len(variant.ALT[0]):
        vartype = 'SNV'
    else:
        vartype = 'INDEL'

    filter = 'PASS'
    if variant.FILTER is not None:
        filter = variant.FILTER

    chr1 = str(variant.CHROM)
    pos1 = variant.POS
    svlen = len(variant.ALT) - len(variant.REF)
        
    gene = list(knowngenelist.intersection(get_genes(transcripts)))[0]
    
    if re.match("synonymous|UTR|stream",genes[gene]['Consequence']):
        continue
    
    abundance = variant.format("VAF")[0][0] * 100
        
    csyntax = 'NA'
    if genes[gene]['HGVSc'] is not None and genes[gene]['HGVSc'] is not '':
       csyntax = genes[gene]['HGVSc'].split(":")[1]
       
    psyntax = 'NA'
    if genes[gene]['HGVSp'] is not None and genes[gene]['HGVSp'] is not '':
        psyntax = genes[gene]['HGVSp'].split(":")[1]
    
    pmaf = genes[gene]['MAX_AF']
    if pmaf is None or pmaf == '':
        pmaf = 'NA'
    else:
        pmaf = str(float(pmaf) * 100) + '%'
        
    psyntax = convert_aa(psyntax)
    
    vars['genelevel'].append([vartype,chr1,str(pos1),variant.REF,variant.ALT[0],gene,genes[gene]['Consequence'],csyntax,psyntax,str(genes[gene]['EXON']),filter,
                              str(variant.ID),str(round(abundance,1))+"%",str(variant.format("NV")[0][0]),str(variant.format("NR")[0][0]),pmaf])
    
    
# done getting gene variants

# get SVs. have to go through twice because BNDs have 2 entries per

passedvars = {} # the ones that passed all filters
alreadydone = set() # so that BND mates can be skipped appropriately 

svvcf = VCF(svvcf_file)
svcsq_fields = parse_csq_header(svvcf)

for variant in svvcf:

    vartype = variant.INFO.get('SVTYPE')

    filter = 'PASS'
    if variant.FILTER is not None:
        filter = variant.FILTER

    # the only filtering done in this script--return only known or novel SVs
    if variant.INFO.get('KNOWNSV') is not None or (filter is 'PASS' and (variant.INFO.get('BLACKLIST_AF') == 0 or variant.INFO.get('BLACKLIST_AF') is None)):
#        if variant.INFO.get('CN') is None or (variant.INFO.get('CN') is not None and variant.INFO.get('CNBINS') >= mincnbins):
        passedvars[variant.ID] = variant

# done with first pass

for v in passedvars.items():

    out = []
    
    variant = v[1]
    
    isknown = 0
    knowngenes = []
    mate = '';

    if variant.INFO.get('KNOWNSV') is not None:
        isknown = 1
        
    vartype = variant.INFO.get('SVTYPE')

    csq = variant.INFO.get('CSQ')

    if csq is None:
        sys.exit("No VEP fields")

    transcripts = list(parse_csq_entries(csq.split(','), svcsq_fields).items())[0][1] # just get the first allele in the list.
    genes = get_csq_entries_bygene(transcripts)
                                    
    filter = 'PASS'
    if variant.FILTER is not None:
        filter = variant.FILTER

    if vartype in ('DEL','DUP'):
        
        chr1 = str(variant.CHROM)
        pos1 = variant.POS
        
        chr2 = chr1
        pos2 = int(variant.INFO.get('END'))
        svlen = pos2 - pos1 + 1
        
        gs = get_genes(transcripts)
        if '' in gs:
            gs.remove('')

        if len(gs) == 0:
            gs.add('None')
            
        bands = transcripts[0]['cytobands'].split("&")
        bandstr = bands[0]
        if len(bands) > 1:
            bandstr = bands[0] + bands[-1]

        knowngenes = gs.intersection(knowngenelist)
        knowngenestring = ",".join(knowngenes)
        if len(knowngenes) == 0:
            knowngenestring = 'None'
            
        genestring = ",".join(gs)
        if len(genestring) > 10:
            genestring = str(len(gs)) + " genes"
            
        # For example:  seq[GRCh38] del(X)(q21.31q22.1)
        #          chrX:g.89555676_100352080del

        csyntax = '.'
        psyntax = '.'
        if vartype == 'DEL':
            csyntax = chr1 + ":g." + str(pos1) + "_" + str(pos2) + "del"
            if bands[0].find('p') > -1 and bands[-1].find('q') > -1: # if the CNA spans the centromere then the whole chromosome is lost/gained
                psyntax = "-" + chr1.replace('chr','')
                
            elif bands[0].find('p') > -1:
                bands.reverse()
                psyntax = "del(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
                
            else:
                psyntax = "del(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
            
        elif vartype == 'DUP':
            csyntax = chr1 + ":g." + str(pos1) + "_" + str(pos2) + "gain"
            if bands[0].find('p') > -1 and bands[-1].find('q') > -1:
                psyntax = "+" + chr1.replace('chr','')
                
            elif bands[0].find('p') > -1:
                bands.reverse()
                psyntax = "del(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
                
            else:
                psyntax = "gain(" + chr1.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
            
        # abundance
        abundance = 0.0
        pr = (0,0)
        sr = (0,0)
        if variant.INFO.get('LOG2RATIO') is not None:
            CN = variant.INFO.get('CN')
            
            if variant.INFO.get('LOG2RATIO') < 0 and CN > 2:
                CN = 1
            
            abundance = ((2**variant.INFO.get('LOG2RATIO') - 1.0) / ((CN/2.0 - 1.0)))*100;

            out = [vartype,chr1,str(pos1),chr2,str(pos2),str(svlen),bandstr,knowngenestring,csyntax,psyntax,genestring,filter,str(variant.ID),str(round(abundance,1))+"%",str(variant.INFO.get('CN')),str(round(variant.INFO.get('LOG2RATIO'),3))]
                
        elif variant.format("SR") is not None and variant.format("PR")[0] is not None:
            sr = variant.format("SR")[0]
            pr =  variant.format("PR")[0]
                
            abundance = (sr[1] + pr[1]) / (pr[0] + pr[1] + sr[0] + sr[1]) * 100
            
            numhits = '.'
            if (variant.INFO.get('CONTIGHITS') is not None):
                numhits = variant.INFO.get('CONTIGHITS')
            
            out = [vartype,chr1,str(pos1),chr2,str(pos2),str(svlen),bandstr,knowngenestring,csyntax,psyntax,genestring,filter,str(numhits),str(variant.ID),str(round(abundance,1))+"%",str(pr[1]) + '/' + str(pr[0]+pr[1]),str(sr[1]) + '/' + str(sr[0]+sr[1]),str(variant.INFO.get('CONTIG'))]
                
    elif vartype == 'BND':

        chr1 = str(variant.CHROM)
        pos1 = variant.POS
                        
        # skip if this is the mate
        if variant.INFO.get('MATEID') in alreadydone:
            continue
            
        # get the mate
        mate = passedvars[variant.INFO.get('MATEID')]

        if mate.INFO.get('KNOWNSV') is not None:
            isknown = 1
            
        matecsq = mate.INFO.get('CSQ')

        if matecsq is None:
            sys.exit("No VEP fields")
                
        matetranscripts = list(parse_csq_entries(matecsq.split(','), svcsq_fields).items())[0][1] # just get the first allele in the list.
        mategenes = get_csq_entries_bygene(matetranscripts)
            
        chr2 = mate.CHROM
        pos2 = mate.POS
        
        gs1 = get_genes(transcripts)
        if '' in gs1:
            gs1.remove('')
            
        if len(gs1) == 0:
            gs1.add('INTERGENIC')
        
        bands1 = transcripts[0]['cytobands'].split("&")
        
        gs2 = get_genes(matetranscripts)
        if '' in gs2:
            gs2.remove('')
            
        if len(gs2) == 0:
            gs2.add('INTERGENIC')
            
        bands2 = matetranscripts[0]['cytobands'].split("&")
        bandstr = bands1[0] + bands2[0]
        
        genestring = ",".join(gs1) + "--" + ",".join(gs2)
                                
        # abundance
        abundance = 0.0
        pr = (0,0)
        sr = (0,0)            
        if variant.format("SR") is not None:
            sr = variant.format("SR")[0]
            
        if variant.format("PR")[0] is not None:                
            pr =  variant.format("PR")[0]

        abundance = (sr[1] + pr[1]) / (pr[0] + pr[1] + sr[0] + sr[1]) * 100

        numhits = '.'
        if (variant.INFO.get('CONTIGHITS') is not None):
            numhits = variant.INFO.get('CONTIGHITS')

        knowngenes1 = list(gs1.intersection(knowngenelist))
        knowngenes2 = list(gs2.intersection(knowngenelist))
        knowngenes = knowngenes1 + knowngenes2

        knowngenestring = ",".join(gs1) + '--' + ",".join(gs2)
        if len(knowngenes) == 0:
            knowngenestring = 'None'
                                
        alt = variant.ALT[0]
        strand = '+'
        if alt.find("[") == 0 or alt.find("]") > 0:
            strand = '-'
            
        csyntax = '';
        psyntax = '';
        if (chr1.find('X') == -1 and chr2.find('X') == -1 and chr1.find('Y') == -1 and chr2.find('Y') == -1 and int(chr1.replace('chr','')) < int(chr2.replace('chr',''))) or chr1.find('X') > -1 or chr1.find('Y') > -1: # this isnt working. Want to list lower chromosome first in these strings. If X is involved, then X first.
            csyntax = chr1 + ":g." + str(pos1) + "(+)::" + chr2 + ":g." + str(pos2) + "(" + strand + ")"
            psyntax = 't(' + chr1.replace('chr','') + ';' + chr2.replace('chr','') + ')(' + bands1[0] + ';' + bands2[0] + ')'
        else:
            csyntax = chr2 + ":g." + str(pos2) + "(+)::" + chr1 + ":g." + str(pos1) + "(" + strand + ")"
            psyntax = 't(' + chr2.replace('chr','') + ';' + chr1.replace('chr','') + ')(' + bands2[0] + ';' + bands1[0] + ')'
                
        out = [vartype,chr1,str(pos1),chr2,str(pos2),"NA",bandstr,knowngenestring,csyntax,psyntax,genestring,filter,str(numhits),str(variant.ID) + ";" + str(mate.ID),str(round(abundance,1))+"%",str(pr[1]) + '/' + str(pr[0]+pr[1]),str(sr[1]) + '/' + str(sr[0]+sr[1]),str(variant.INFO.get('CONTIG'))]

        alreadydone.add(variant.ID)
        
    if isknown == 1:
        vars['knownsv'].append(out)
        
    elif variant.INFO.get('LOG2RATIO') is not None:
        vars['cna'].append(out)

    elif len(knowngenes) > 0:
        vars['knowngenes'].append(out)
        
    else:
        vars['novelsv'].append(out)


print("ChromoSeq Report for " + Name + " ---- Generated on: " + strftime("%Y-%m-%d %H:%M:%S", gmtime()) + "\n")

print("*** COPY NUMBER ALTERATIONS ***\n")
if len(vars['cna']) > 0:
    print("\t".join(('TYPE','CHR1','POS1','CHR2','POS2','LENGTH','BANDS','KNOWN_GENES','HGVS-LIKE','ISCN-LIKE','TOTAL_GENES','FILTERS','ID','ABUNDANCE','PLOIDY','LOG2RATIO')))
    for i in vars['cna']:
        print("\t".join(i))
    print()        
else:
    print("\tNONE DETECTED\n")    

print("*** RECURRENT TRANSLOCATIONS ***\n")
if len(vars['knownsv']) > 0:
    print("\t".join(('TYPE','CHR1','POS1','CHR2','POS2','LENGTH','BANDS','KNOWN_GENES','HGVS-LIKE','ISCN-LIKE','TOTAL_GENES','FILTERS','CONTIG_HITS','ID','ABUNDANCE','SR_READS','PR_READS','CONTIG')))
    for i in vars['knownsv']:
        print("\t".join(i))
    print()
else:
    print("\tNONE DETECTED\n")    

print("*** GENE MUTATIONS ***\n")
if len(vars['genelevel']) > 0:
    print("\t".join(('TYPE','CHR','POS','REF','ALT','GENE','CONSEQUENCE','HGVSc','HGVSp','EXON','FILTERS','ID','VAF','VariantReads','TotalReads','PopulationMAF')))
    for r in vars['genelevel']:
        print("\t".join(r))
    print()
else:
    print("\t"+"NONE DETECTED\n")

print("*** NOVEL STRUCTURAL VARIANTS INVOLVING KNOWN GENES ***\n")
if len(vars['knowngenes']) > 0:
    print("\t".join(('TYPE','CHR1','POS1','CHR2','POS2','LENGTH','BANDS','KNOWN_GENES','HGVS-LIKE','ISCN-LIKE','TOTAL_GENES','FILTERS','CONTIG_HITS','ID','ABUNDANCE','PR_READS','SR_READS','CONTIG')))
    for r in vars['knowngenes']:
        print("\t".join(r))
    print()
else:
    print("\t"+"NONE DETECTED\n")

print("*** OTHER HIGH-CONFIDENCE STRUCTURAL VARIANTS ***\n")
if len(vars['novelsv']) > 0:
    print("\t".join(('TYPE','CHR1','POS1','CHR2','POS2','LENGTH','BANDS','KNOWN_GENES','HGVS-LIKE','ISCN-LIKE','TOTAL_GENES','FILTERS','CONTIG_HITS','ID','ABUNDANCE','PR_READS','SR_READS','CONTIG')))
    for r in vars['novelsv']:
        print("\t".join(r))
    print()
else:
    print("\t"+"NONE DETECTED\n")

