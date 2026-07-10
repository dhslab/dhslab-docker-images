#!/usr/bin/env python3

from io import StringIO
import sys, os, re, tempfile, csv, pysam, json, binascii, argparse, subprocess, gzip
import pandas as pd
import pyranges as pr
import numpy as np
from time import gmtime, strftime
from cyvcf2 import VCF
from pathlib import Path

__version__ = '2.0.0'

ACROCENTRICS = ['13','14','15','21','22']

def checkfile(file_path):
    """Check if a file exists at the given path."""
    if not os.path.exists(file_path):
        raise argparse.ArgumentTypeError(f"The file {file_path} does not exist.")
    return file_path

def fileexists(file_path):
    """Check if a file exists at the given path."""
    if os.path.exists(file_path):
        raise argparse.ArgumentTypeError(f"The outfile {file_path} exists!")
    return file_path

def revcomp(dna):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}  # DNA complement pairs
    reverse_complement = "".join(complement.get(base, base) for base in reversed(dna))
    return reverse_complement

def pos2codon(exonstart,exonend,cpos,pos,strand):
    if pos<=exonend and pos>=exonstart:
        if strand == '+':
            return(int((pos-exonstart + cpos) / 3 + .99))

        elif strand=='-':
             return(int((exonend-pos + cpos) / 3 + .99))

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

def sort_chrompos(row,chrom='Chromosome',pos='Start'):
    # Extract the numeric part from the 'Chromosome' value
    if row[chrom].startswith('chr'):
        num_part = row[chrom][3:]
        if num_part.isdigit():
            chromosome_value = int(num_part)
        else:
            chromosome_value = float('inf')
    else:
        chromosome_value = float('inf')
    
    return (chromosome_value, row[pos])

def checkQc(value, minmax):
    if len(minmax)==1:
        if value < float(minmax[0]):
            return ['>' + str(minmax[0]),'(!)']
        else: 
            return ['>' + str(minmax[0]),'']
    elif len(minmax)==2:
        if value > float(minmax[0]) and value < float(minmax[1]):
            return [str(minmax[0])+'-'+str(minmax[1]),'']
        else:
            return [str(minmax[0])+'-'+str(minmax[1]), '(!)']
        
    else:
        sys.exit("Reference range must have minimum and or maximum")

def get_runinfo(cramfile,reference=None):

    runinfo = {'runid':[],'instrument':[],'flowcell':[],'flowcell_type':[],'flowcell_lot':[],'reagent_lot':[],'recipe':[]}

    # Open the BAM/CRAM file
    with pysam.AlignmentFile(cramfile, "rb", reference_filename=reference) as samfile:
        # Access the header
        header = samfile.header

        # Iterate through read groups in the header
        for read_group in header.get('RG', []):
            if 'PL' in read_group:
                rgpl = read_group['PL'].split('.')
                if len(rgpl) == 7:
                    runid, instrument, flowcellid, flowcelltype, flowcelllot, reagentlot, recipe = rgpl
                    runinfo['runid'].append(runid)
                    runinfo['instrument'].append(instrument)
                    runinfo['flowcell'].append(flowcellid)
                    runinfo['flowcell_type'].append(flowcelltype)
                    runinfo['flowcell_lot'].append(flowcelllot)
                    runinfo['reagent_lot'].append(reagentlot)
                    runinfo['recipe'].append(recipe)

        if len(runinfo['runid']) == 0:
            for read in samfile.fetch():
                readName = read.query_name
                parts = readName.split(':')
                runinfo['runid'].append(f"RUN_{parts[0]}_{str(int(parts[1])).zfill(4)}_{parts[2]}")
                runinfo['instrument'].append(parts[0])
                runinfo['flowcell'].append(parts[2])
                break

        for key in runinfo:
            if len(runinfo[key]) > 0 and runinfo[key]!='?':
                runinfo[key] = ','.join(set(runinfo[key]))
            else:
                runinfo[key] = 'NONE'

    return runinfo

def vepToTable(csq,header):
    fields = header['Description'].strip('"').split("|")
    data = []
    for i in csq.split(','):
        data.append(dict(zip(fields, i.split("|"))))
    df = pd.DataFrame(data)

    # if no symbol, use Gene ID
    df.loc[df['SYMBOL']=='','SYMBOL'] = df.loc[df['SYMBOL']=='','Gene']

    df.loc[df['STRAND']=='1','STRAND'] = '+'
    df.loc[df['STRAND']=='-1','STRAND'] = '-'

    if 'DISTANCE' in df.columns:
        df['DISTANCE'] = df['DISTANCE'].apply(lambda x: 0 if x=='' else int(x))

    if 'PICK' in df.columns:
        df['PICK'] = df['PICK'].apply(lambda x: 0 if x=='' else 1)

    return df

def getVepFields(vcf):
    vep = {}
    i = 0
    for j in vcf.get_header_type('CSQ')['Description'].strip('"').split("|"):
        vep[j] = i
        i+=1
    return(vep)

def vepGeneEffect(row):
    if row['EXON']:
        return('(exon' + row['EXON'].split('/')[0] + ')')
    elif row['INTRON']:
        return('(intron' + row['INTRON'].split('/')[0] + ')')
    elif row['DISTANCE']:
        if 'upstream' in row['Consequence']:
            return('(' + str(row['DISTANCE']) + 'bp upstream)')
        elif 'downstream' in row['Consequence']:
            return('(' + str(row['DISTANCE']) + 'bp downstream)')
        else:
            return('')
    else:
        return('')
    

def df_to_dict_nan_to_none(df, index=False):
  return df.replace({np.nan: None}).to_dict('split', index=index)

def parse_small_variants(vcffile,individual=[0]):
    
    #########################################
    #
    # Get small variants
    #
    #########################################

    variants_data = []
    variants_columns = ['type','filters','chrom','pos','ref','alt','gene','gene_id','transcript','consequence','csyntax','psyntax','exon','intron','pop_af','annotations']

    vcf = VCF(vcffile)

    all_sample_names = vcf.samples
    sample_names = [all_sample_names[i] for i in individual]
    variants_columns += [
        f"{sample}_{suffix}" 
        for sample in sample_names
        for suffix in ["coverage", "altreads", "vaf"]
    ]

    # get VEP fields
    vepFields = getVepFields(vcf)

    # get variants
    for variant in vcf:

        vartype = ''
        if len(variant.REF) == len(variant.ALT[0]):
            vartype = 'SNV'
        else:
            vartype = 'INDEL'

        varfilter = 'PASS'
        if variant.FILTER is not None:
            continue
            #varfilter = variant.FILTER
    
        # adding 241101, skip variants if components of MNV via MNV_tag
        if variant.INFO.get('MNVTAG'):
            # get MNV_tag value
            mnv_tag = variant.INFO['MNVTAG']
            # format is chrom:pos_ref->alt
            mnv_alt = mnv_tag.split(">")[-1]
            # continue if this alt doesn't match the current alt
            if mnv_alt != variant.ALT[0]:
                continue

        sample_data = ['NA','NA','NA']
        # multiply sample_data by the number of samples
        sample_data = sample_data * len(sample_names)

        for i, sample in enumerate(sample_names):
            sample_data[i*3] = variant.format("DP")[i][0]
            sample_data[i*3+1] = variant.format("AD")[i][1]
            sample_data[i*3+2] = round(variant.format("AF")[i][0] * 100,2)
            
        # get VEP annotation
        csq = variant.INFO['CSQ']
        
        if csq is None:
            sys.exit("No VEP fields")
        
        gene='NA'
        gene_id='NA'
        transcript='NA'
        csyntax='NA'
        psyntax='NA'
        consequence='NA'
        exon='NA'
        intron='NA'
        popmaf = 'NA'
        customannotation = 'NA'
        for i in variant.INFO['CSQ'].split(','):
            csq = i.split("|")

            if csq[vepFields['PICK']] == '':
                continue

            # get pop allele frequency. This is present for each transcript annotation, but is always the same 
            if csq[vepFields['MAX_AF']] != '':
                popmaf = float(csq[vepFields['MAX_AF']])

            transcript = csq[vepFields['Feature']]
            gene = csq[vepFields['SYMBOL']]
            gene_id = csq[vepFields['Gene']]
            consequence = csq[vepFields['Consequence']].split("&")[0]

            csyntax = csq[vepFields['HGVSc']].split(":")
            if len(csyntax) > 1:
                csyntax = csyntax[1]
            else:
                csyntax = consequence

            psyntax = csq[vepFields['HGVSp']].split(":")
            if len(psyntax) > 1:
                psyntax = convert_aa(psyntax[1])
                psyntax = re.sub("\%3D","=",psyntax)
            else:
                psyntax = csyntax
            
            impact = csq[vepFields['IMPACT']]
            exon = csq[vepFields['EXON']] or 'NA'
            intron = csq[vepFields['INTRON']] or 'NA'
            customannotation = csq[vepFields['Existing_variation']] or 'NA'    

        # convert pop maf to percent
        if popmaf!='NA':
            popmaf = round(float(popmaf)*100,3)

        # only include all variants <=0.1% and ns or specific noncoding variants 
        variants_data.append(dict(zip(variants_columns,[vartype,varfilter,str(variant.CHROM),variant.POS,variant.REF,variant.ALT[0],gene,gene_id,transcript,consequence,csyntax,psyntax,exon,intron,popmaf,customannotation]+sample_data)))
    
    variants = pd.DataFrame(variants_data, columns=variants_columns)
    return variants

def parse_svs(svvcffile,individual=0):

    nonSynon = ["splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant","stop_lost","start_lost","transcript_ablation","transcript_amplification","inframe_insertion","inframe_deletion","missense_variant","protein_altering_variant"]

    svs_data = []
    svs_columns = ['type','chrom1','pos1','chrom2','pos2','length','csyntax','psyntax','bands','known_genes','known_gene_detail','total_genes','filters','id','abundance','info']

    ########################
    #
    # Get SVs
    # 
    ########################

    print("Gathering SVs...",file=sys.stderr)

    svvcf = VCF(svvcffile)

    passedvars = {}

    for variant in svvcf:

        # save BNDs because we need both records for these
        if variant.INFO['SVTYPE']=='BND' or variant.INFO['SVTYPE']=='INV':
            passedvars[variant.ID] = variant

        else:

            filter = []
            if variant.FILTER is not None:
                filter = variant.FILTER.split(';')

            vartype = variant.INFO.get('SVTYPE')

            # get gene 1 info
            chr = str(variant.CHROM)
            pos1 = variant.POS
            pos2 = variant.INFO.get('END')
            svlen = 'NA'
            svlen_info = variant.INFO.get('SVLEN')
            if svlen_info is not None:
                if isinstance(svlen_info, tuple) and len(svlen_info) > 0:
                    svlen = abs(svlen_info[0])
                else:
                    svlen = abs(svlen_info)

            csyntax = ''
            psyntax = ''

            bandstring = 'None'
            genestring = None
            genedetail = None
            total_genes = 'None'

            # get genes for this variant
            # Skip if CSQ annotation is missing
            try:
                csq_data = variant.INFO['CSQ']
            except KeyError:
                # No CSQ annotation, skip this variant
                continue

            vepCsq = vepToTable(csq_data,svvcf.get_header_type('CSQ'))

            total_genes = f"{len(vepCsq['SYMBOL'].unique())} genes"

            bands = vepCsq['cytobands'][0].split("&")
            bandstring = bands[0]
            if len(bands) > 1:
                bandstring = bands[0] + bands[-1]

            vepCsq = vepCsq.sort_values(by=['PICK'], ascending=[False]).drop_duplicates(subset='SYMBOL',keep='first')
            vepCsq['GeneEffect'] = vepCsq.apply(lambda r: vepGeneEffect(r),axis=1)
            vepCsq['GeneImpact'] = vepCsq['Consequence'].apply(lambda r: int(len(set(r.split('&')) & set(nonSynon))>0))

            if len(vepCsq) > 10:
                genestring = f'{len(vepCsq)} genes'
                genedetail = genestring
            elif len(vepCsq) > 0:
                genestring = ','.join(vepCsq['SYMBOL'].unique().tolist())
                genedetail = vartype.lower() + '[' + ','.join(vepCsq.apply(lambda x: x['SYMBOL'] + x['GeneEffect'],axis=1).tolist()) + ']'

            if vartype == 'DEL':
                csyntax = chr + ":g." + str(pos1) + "_" + str(pos2) + "del"
                if bands[0].find('p') > -1 and bands[-1].find('q') > -1: # if the CNA spans the centromere then the whole chromosome is lost/gained
                    psyntax = "seq[GRCh38] -" + chr.replace('chr','')
                    
                elif 'q11' in bands and 'qter' in bands and chr1 in ACROCENTRICS:
                    psyntax = "seq[GRCh38] -" + chr.replace('chr','')

                elif bands[0].find('p') > -1:
                    bands.reverse()
                    psyntax = "seq[GRCh38] del(" + chr.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
                    
                else:
                    psyntax = "seq[GRCh38] del(" + chr.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
            
            elif vartype == 'DUP' or vartype == 'INS':
                csyntax = chr + ":g." + str(pos1) + "_" + str(pos2) + vartype.lower()
                if bands[0].find('p') > -1 and bands[-1].find('q') > -1:
                    psyntax = "seq[GRCh38] +" + chr.replace('chr','')

                elif 'q11' in bands and 'qter' in bands and chr in ACROCENTRICS:
                    psyntax = "seq[GRCh38] +" + chr.replace('chr','')

                elif bands[0].find('p') > -1:
                    bands.reverse()
                    psyntax = "seq[GRCh38] " + vartype.lower() + "(" + chr.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"
                    
                else:
                    psyntax = "seq[GRCh38] " + vartype.lower() + "(" + chr.replace('chr','') + ")(" + bands[0] + bands[-1] + ")"

            abundance = 0.0
            CN='.'
            PR = (0,0)
            SR = (0,0)            
            if "SR" in variant.FORMAT:
                SR = variant.format("SR")[individual]
                
            if "PR" in variant.FORMAT:
                PR =  variant.format("PR")[individual]

            if "MAF" in variant.FORMAT and "CN" in variant.FORMAT:
                if len(variant.format("MAF"))==1:
                    abundance = round(variant.format("MAF")[0][0]* 100,2)
                    CN = variant.format("CN")[0][0]
                else:
                    abundance = round(variant.format("MAF")[individual][0]* 100,2)
                    CN = variant.format("CN")[individual][0]

            elif PR[0] + SR[0] + PR[1] + SR[1] > 0:
                abundance = round((SR[1] + PR[1]) / (PR[0] + PR[1] + SR[0] + SR[1])*100,1)

            infostring = 'CN=' + str(CN) + ';PR_READS=' + str(PR[1]) + '/' + str(PR[0]+PR[1]) + ';SR_READS=' + str(SR[1]) + '/' + str(SR[0]+SR[1]) + ';CONTIG=' + str(variant.INFO.get('CONTIG'))

            if len(filter) == 0:
                filter = 'PASS'
            else:
                filter = ';'.join(filter)

            svs_data.append(dict(zip(svs_columns,[vartype,chr,pos1,chr,pos2,svlen,csyntax,psyntax,bandstring,genestring,genedetail,total_genes,filter,str(variant.ID),abundance,infostring])))

    # now handle BNDs, which each have 2 entries.
    # this includes translocations and inversions
    alreadydone = set()

    i  = 0
    for v in passedvars.items():

        # get the first record
        variant = v[1]

        # skip if this is the mate and we already processed the record pair
        if variant.INFO.get('MATEID') in alreadydone or variant.INFO.get('MATEID') not in passedvars:
            continue

        # get the mate
        mate = passedvars[variant.INFO.get('MATEID')]

        vartype = variant.INFO.get('SVTYPE')
        
        # get filter info
        filter = []
        if variant.FILTER is not None:
            filter = variant.FILTER.split(';')
        
        # Info for variant (first end of BND)
        chr1 = str(variant.CHROM)
        pos1 = variant.POS
        genes1='NA'
        transcript1='NA'
        region1 = 'NA'
        strand1 = '+'
        bands1 = 'NA'

        # Info for mate (second end of BND)
        chr2 = mate.CHROM
        pos2 = mate.POS     
        genes2='NA'
        transcript2='NA'
        region2 = 'NA'
        strand2 = '+'
        bands2 = 'NA'

        # variant info.
        csyntax = 'None'
        psyntax = 'None'
        bandstring = 'None'
        genestring = 'None'
        genedetail = 'None'
        total_genes = '2 genes'

        # get gene info for VARIANT
        # Skip if CSQ annotation is missing (e.g., stub data without VEP annotations)
        try:
            csq_data = variant.INFO['CSQ']
        except KeyError:
            # No CSQ annotation, skip this BND pair
            alreadydone.add(variant.ID)
            continue

        vepCsq = vepToTable(csq_data,svvcf.get_header_type('CSQ'))

        bands1 = vepCsq['cytobands'][0].split("&")[0]

        vepCsq = vepCsq.sort_values(by=['PICK'], ascending=[False]).drop_duplicates(subset='SYMBOL',keep='first')
        vepCsq['GeneEffect'] = vepCsq.apply(lambda r: vepGeneEffect(r),axis=1)
        gene1Df= vepCsq[['SYMBOL','DISTANCE','STRAND','GeneEffect']]

        # sort by known transcript and then distance and get only the first transcript.
        gene1Df = gene1Df.sort_values(by=['DISTANCE','GeneEffect'], ascending=[True,False],na_position='last').fillna('')

        # get gene info for MATE
        # Skip if mate has no CSQ annotation
        try:
            mate_csq_data = mate.INFO['CSQ']
        except KeyError:
            # No CSQ annotation on mate, skip this BND pair
            alreadydone.add(variant.ID)
            continue

        vepCsq = vepToTable(mate_csq_data,svvcf.get_header_type('CSQ'))    

        bands2 = vepCsq['cytobands'][0].split("&")[0]

        vepCsq = vepCsq.sort_values(by=['DISTANCE','PICK'], ascending=[True,False]).drop_duplicates(subset='SYMBOL',keep='first')
        vepCsq['GeneEffect'] = vepCsq.apply(lambda r: vepGeneEffect(r),axis=1)
        gene2Df= vepCsq[['SYMBOL','DISTANCE','STRAND','GeneEffect']]

        # sort by known transcript and then distance and get only the first transcript.
        gene2Df = gene2Df.sort_values(by=['DISTANCE','GeneEffect'], ascending=[True,False],na_position='last').fillna('')
        
        if gene1Df.shape[0] > 0:
            genes1 = gene1Df['SYMBOL'].tolist()[0]
            region1 = gene1Df['GeneEffect'].tolist()[0]
            strand1 = gene1Df['STRAND'].tolist()[0]

        if gene2Df.shape[0] > 0:
            genes2 = gene2Df['SYMBOL'].tolist()[0]
            region2 = gene2Df['GeneEffect'].tolist()[0]
            strand2 = gene2Df['STRAND'].tolist()[0]

        if genes1=='NA' or genes1=='':
            genes1 = 'INTERGENIC'

        if genes2=='NA' or genes2=='':
            genes2 = 'INTERGENIC'

        orientation = '+'
        if variant.ALT[0].find("[") == 0 or variant.ALT[0].find("]") > 0:
            orientation = '-'    

        # abundance
        abundance = 0.0
        PR = (0,0)
        SR = (0,0)            
        if variant.format("SR") is not None:
            SR = variant.format("SR")[0]
                
        if variant.format("PR")[0] is not None:                
            PR =  variant.format("PR")[0]

        if 'DUX4:IGH' in variant.ID:
            SR[0] = 0
            PR[0] = 0

        supporting_reads = SR[1] + PR[1]
        abundance = round((SR[1] + PR[1]) / (PR[0] + PR[1] + SR[0] + SR[1]) * 100,1)

        alt = variant.ALT[0]
        orientation = '+'
        if alt.find("[") == 0 or alt.find("]") > 0:
            orientation = '-'

        svtype = 't'
        if vartype == 'INV':
            svtype = 'inv'

        if (chr1.find('M') == -1 and chr2.find('M') == -1 and chr1.find('X') == -1 and chr2.find('X') == -1 and chr1.find('Y') == -1 and chr2.find('Y') == -1 and int(chr1.replace('chr','')) < int(chr2.replace('chr',''))) or chr1.find('X') > -1 or chr1.find('Y') > -1:
            csyntax = chr1 + ":g." + str(pos1) + "(+)::" + chr2 + ":g." + str(pos2) + "(" + orientation + ")"
            psyntax = f'seq[GRCh38] {svtype}(' + chr1.replace('chr','') + ';' + chr2.replace('chr','') + ')(' + bands1 + ';' + bands2 + ')'
            bandstring = bands1 + ';' + bands2
            genestring = genes1+'::'+genes2

            if genes1 == 'INTERGENIC':
                genedetail = genes1+region1+'::'+genes2+'('+strand2+')'+region2
            elif genes2 == 'INTERGENIC':
                genedetail = genes1+'('+strand1+')'+region1+'::'+genes2+region2
            else:
                genedetail = genes1+'('+strand1+')'+region1+'::'+genes2+'('+strand2+')'+region2

        else:
            csyntax = chr2 + ":g." + str(pos2) + "(+)::" + chr1 + ":g." + str(pos1) + "(" + orientation + ")"
            psyntax = f'seq[GRCh38] {svtype}(' + chr2.replace('chr','') + ';' + chr1.replace('chr','') + ')(' + bands2 + ';' + bands1 + ')'
            bandstring = bands2 + ';' + bands1
            genestring = genes2+'::'+genes1

            if genes1 == 'INTERGENIC':
                genedetail = genes2+'('+strand2+')'+region2+'::'+genes1+region1
            elif genes2 == 'INTERGENIC':
                genedetail = genes2+region2+'::'+genes1+'('+strand1+')'+region1
            else:
                genedetail = genes2+'('+strand2+')'+region2+'::'+genes1+'('+strand1+')'+region1

        infostring = 'PR_READS=' + str(PR[1]) + '/' + str(PR[0]+PR[1]) + ';SR_READS=' + str(SR[1]) + '/' + str(SR[0]+SR[1]) + ';CONTIG=' + str(variant.INFO.get('CONTIG'))

        if len(filter) == 0:
            filter = 'PASS'
        else:
            filter = ';'.join(filter)

        svs_data.append(dict(zip(svs_columns,[vartype,chr1,pos1,chr2,pos2,svlen,csyntax,psyntax,bandstring,genestring,genedetail,total_genes,filter,str(variant.ID) + ";" + str(mate.ID),abundance,infostring])))

        alreadydone.add(variant.ID)
    
    svs = pd.DataFrame(svs_data, columns=svs_columns)
    return svs

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Script
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parser = argparse.ArgumentParser(description='Vep to table')
parser.add_argument('-v','--vcf',required=False,type=str,help='Small variant vcf')
parser.add_argument('-s','--svvcf',required=False,type=str,help='SV vcf')
parser.add_argument('-i','--individual',required=False,default=[0],type=int,nargs='+',help='Individual in the VCF to return format fields')
parser.add_argument('-o','--outfile',required=True,type=str,help='Output file')

args = parser.parse_args()

sample = list(args.individual)[0] if len(args.individual) == 1 else args.individual

if args.vcf is not None:
    variants = parse_small_variants(args.vcf,individual=sample)
    variants.to_csv(args.outfile,sep='\t',index=False)

if args.svvcf is not None:
    svs = parse_svs(args.svvcf,individual=sample)
    svs.to_csv(args.outfile,sep='\t',index=False)

    
