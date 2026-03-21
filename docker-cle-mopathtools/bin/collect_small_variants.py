#!/usr/bin/env python

import argparse
import binascii
import json
import os
import re
import sys
from math import ceil as ceiling
from io import StringIO
from pathlib import Path
from time import gmtime, strftime

import numpy as np
import pandas as pd
import pysam
import natsort

__version__ = "1.0.0"

def df_to_dict_nan_to_none(df, index=False):
    return df.replace({np.nan: None}).to_dict("split", index=index)

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

def decode_hex(string):
    hex_string = string.group(0).replace("%", "")
    return binascii.unhexlify(hex_string).decode("utf-8")

def convert_aa(codon):
    three = [
        "Ala",
        "Arg",
        "Asn",
        "Asp",
        "Cys",
        "Glu",
        "Gln",
        "Gly",
        "His",
        "Ile",
        "Leu",
        "Lys",
        "Met",
        "Phe",
        "Pro",
        "Ser",
        "Thr",
        "Trp",
        "Tyr",
        "Val",
        "Ter",
    ]
    one = [
        "A",
        "R",
        "N",
        "D",
        "C",
        "E",
        "Q",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
        "*",
    ]

    for i in range(0, len(three)):
        p = re.compile(three[i])
        codon = p.sub(one[i], codon)

    return codon

# write a function to parse the VEP CSQ VCF header
def vepHeader(header):
    fields = (
        header.info.get("CSQ")
        .record.get("Description")
        .split(":")[1]
        .replace('"', "")
        .strip()
        .split("|")
    )
    return dict(zip(fields, range(0, len(fields))))

def _get_fields_from_header(description: str) -> list[str]:
    """
    Parses a VCF header description string to extract field names.
    This version removes all whitespace, splits by ':', and takes the second element.

    Example format: "Description: FIELD1|FIELD2|FIELD3"
    """
    try:
        # Remove all whitespace, split by ':', take the second part, and split by '|'
        return description.strip('"').strip('.').strip(' ').split(':')[1].strip(' ').split('|')
    except IndexError:
        # Return an empty list if the format is unexpected (e.g., no ':')
        return []

def vep_csq_to_dataframe(csq_tag: str | None, header: dict) -> pd.DataFrame:
    """
    Parses a VEP CSQ string from a VCF INFO tag into a DataFrame.

    Args:
        csq_tag: The CSQ string from the VCF INFO field.
        header: The VCF header dictionary for the CSQ field.

    Returns:
        A DataFrame of VEP annotations.
    """
    fields = _get_fields_from_header(header.get("Description", ""))
    
    # Return an empty DataFrame with the correct columns if the tag is missing
    if not csq_tag or not fields:
        return pd.DataFrame(columns=fields + ["KnownGene", "KnownTrx"])
        
    # 1. Create a list of dictionaries from the CSQ string
    annotations = [dict(zip(fields, entry.split("|"))) for entry in csq_tag.split(",")]
    
    # 2. Convert the list to a DataFrame
    df = pd.DataFrame(annotations)

    # --- 3. Data Cleaning and Type Conversion ---

    # remove leading or trailing '.' from Allele column 
    df["Allele"] = df["Allele"].str.strip('.')
    
    # If a SYMBOL is missing, use the Gene ID instead
    df["SYMBOL"] = df.apply(lambda row: row['Gene'] if row['SYMBOL'] == '' else row['SYMBOL'], axis=1)
    
    # Standardize STRAND representation
    df["STRAND"] = df["STRAND"].replace({"1": "+", "-1": "-"})

    # if TranscriptCoordinates (START and END) are in columns
    if "START" in df.columns and "END" in df.columns:
        # make START and END numeric
        df["START"] = pd.to_numeric(df["START"])
        df["END"] = pd.to_numeric(df["END"])
    
    # Convert DISTANCE to a numeric type, defaulting to 0 for empty strings
    if "DISTANCE" in df.columns:
        df["DISTANCE"] = pd.to_numeric(df["DISTANCE"], errors='coerce').fillna(0).astype(int)

    # Convert PICK to a binary indicator
    if "PICK" in df.columns:
        df["PICK"] = df["PICK"].apply(lambda x: 1 if x == "1" else 0)
    
    # Add custom boolean flags
    df["KnownGene"] = 0
    df["KnownTrx"] = 0

    return df

def known_sv_genes_tag_to_dataframe(variant, header: dict) -> pd.DataFrame:
    """
    Parses the KnownSvGenes INFO tag into a DataFrame with VEP-like columns.

    Args:
        variant: A variant object with POS, INFO attributes.
        header: The VCF header dictionary for the KnownSvGenes field.

    Returns:
        A DataFrame of known gene overlaps with calculated consequences.
    """
    fields = _get_fields_from_header(header.get("Description", ""))
    known_sv_tag = variant.INFO.get("KnownSvGenes")
    
    # Return an empty DataFrame if the tag is missing or variant positions are invalid
    pos1 = variant.POS
    pos2 = variant.INFO.get("END", pos1)
    if not known_sv_tag or not fields or pos1 is None or pos2 is None:
        return pd.DataFrame(columns=fields + ["DISTANCE", "Consequence","KnownGene","KnownTrx"])

    # 1. Create a list of dictionaries from the KnownSvGenes string
    parsed_genes = [dict(zip(fields, entry.split("|"))) for entry in known_sv_tag.split(",")]

    # 2. Convert the list to a DataFrame
    df = pd.DataFrame(parsed_genes)

    # --- 3. Data Cleaning and Type Conversion ---
    # Convert position columns to numeric types
    df["START"] = pd.to_numeric(df["START"])
    df["END"] = pd.to_numeric(df["END"])
    
    # Standardize STRAND representation
    df["STRAND"] = df["STRAND"].replace({"1": "+", "-1": "-"})

    # Add custom boolean flags
    df["KnownGene"] = 1
    df["KnownTrx"] = 1
    
    # --- 4. Custom Calculations ---
    # Calculate DISTANCE to the variant
    df["DISTANCE"] = df.apply(
        lambda row: 0 if (pos1 <= row["END"] and row["START"] <= pos2)
        else min(abs(pos1 - row["END"]), abs(row["START"] - pos2), abs(pos1 - row["START"]), abs(pos2 - row["END"])),
        axis=1,
    )
    
    # Determine the Consequence based on overlap and strand
    def get_consequence(row):
        if row["DISTANCE"] == 0 and (row["END"] > pos1 or row["START"] < pos2): # breakpoint is within gene
            return "transcript_ablation"
        if row["DISTANCE"] == 0 and row["START"] > pos1 and row["END"] < pos2: # SV spans gene
            # if variant.INFO.get("SVTYPE") == DEL return deletion if DUP return transcript_amplification
            if variant.INFO.get("SVTYPE") == "DEL":
                return "deletion"
            else:
                return "transcript_amplification"
        if (row["STRAND"] == "+" and row["START"] > pos2) or (row["STRAND"] == "-" and row["END"] < pos1):
            return "upstream_gene_variant"
        return "downstream_gene_variant"
    df["Consequence"] = df.apply(get_consequence, axis=1)

    return df
 
def getVepFields(vcf):
    vep = {}
    i = 0
    for j in vcf.get_header_type("CSQ")["Description"].strip('"').split("|"):
        vep[j] = i
        i += 1
    return vep

def vepGeneEffect(row):
    # if row is all None, return None
    if row.isna().all():
        return None
    
    if row["EXON"] and '/' in row["EXON"]:
        return f"exon{row['EXON'].split('/')[0]}"
    elif row["INTRON"] and '/' in row["INTRON"]:
        return f"intron{row['INTRON'].split('/')[0]}"
    elif row["DISTANCE"] and row["DISTANCE"]>0:
        if "upstream" in row["Consequence"]:
            return f"{row['DISTANCE']}bp upstream"
        elif "downstream" in row["Consequence"]:
            return f"{row['DISTANCE']}bp downstream"
        else:
            return ""
    else:
        return "intragenic"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Script
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def main():
    parser = argparse.ArgumentParser(description="Collect small variants and return a tab-delimited file")
    parser.add_argument("--vcf-file", required=True, type=str, help="Small variant VCF file")
    parser.add_argument("--target-bed", "-t", required=True, type=str, help="Mopath-style Target BED file")
    parser.add_argument("--qc-file", "-q", required=True, type=str, help="Mopath-style QC JSON file")
    
    parser.add_argument("--outfile", "-o", required=True, type=str, help="Output file")

    # Add argument for VEP blacklist tag
    parser.add_argument(
        "--blacklist-tag", "-b", required=False, type=str, default=None, help="VEP blacklist tag"
    )

    # Add argument for VEP custom annotation tag
    parser.add_argument(
        "--custom-annotation-tag", required=False, type=str, default=None, help="Info tag for custom annotation used for VEP."
    )

    # Add argument for Variant database tag
    parser.add_argument(
        "--variant-db-vcf-tag", required=False, type=str, default=None, help="Info tag for variant database."
    )

    # population allele frequency threshold for categorizing variants as SNPs
    parser.add_argument(
        "--pop-af", "-p", required=False, type=float, default=0.01, help="Population allele frequency threshold"
    )

    parser.add_argument(
        "--version", "-v", action="version", version="%(prog)s: " + __version__
    )

    args = parser.parse_args()


    #########################################
    #
    # Get thresholds and reference values for QC Metrics
    #
    #########################################

    qcranges = {}
    with open(args.qc_file, "r") as json_file:
        qcranges = json.load(json_file)

    # get nonsynonymous annotations
    nonSynon = qcranges["PARAMETERS"]["NONSYNONYMOUS_ANNOTATIONS"]
    nonCodingGeneList = qcranges["PARAMETERS"]['NONCODINGVARIANTLIST']

    # dataframe with small variants
    variants = pd.DataFrame(
        columns=[
            "category",
            "type",
            "filters",
            "chrom",
            "pos",
            "ref",
            "alt",
            "gene",
            "gene_id",
            "transcript",
            "consequence",
            "csyntax",
            "psyntax",
            "exon",
            "intron",
            "pop_af",
            "annotations",
            "coverage",
            "altreads",
            "vaf",
            "priorvariants"
        ]
    )

    # Read in target file to get gene/transcript info
    targetDf = pd.read_csv(
        args.target_bed,
        header=None,
        names=["Chromosome", "Start", "End", "Gene", "Info"],
        sep="\t"
    )
    geneCovDf = targetDf[targetDf["Info"].str.contains(r"^gene", regex=True)]
    svCovDf = targetDf[targetDf["Info"].str.contains(r"^sv", regex=True)]

    #
    # get transcript info for genes and sv genes
    #
    ids = (
        geneCovDf["Info"]
        .str.split("|", expand=True)
    )
    ids.columns = ["RegionType", "Region", "Gene", "GeneId", "Transcript", "cdsStart", "cdsEnd", "strand"]

    geneTrx = (
        pd.concat([geneCovDf, ids[["Region", "Transcript", "cdsStart", "cdsEnd", "strand"]]], axis=1)
        .drop_duplicates()
        .reset_index()
    )
    geneTrx["exonStart"] = geneTrx.apply(lambda r: r["Start"]+3 if r["strand"]=='1' else r["End"]-2, axis=1)
    geneTrx["exonEnd"] = geneTrx.apply(lambda r: r["End"]-2 if r["strand"]=='1' else r["Start"]+3, axis=1)
    geneTrx["cdsStart"] = geneTrx["cdsStart"].astype(int)
    geneTrx["cdsEnd"] = geneTrx["cdsEnd"].astype(int)

    exonTrx = geneTrx[~geneTrx["Region"].str.contains("region")].copy().reset_index(drop=True)
    exonTrx["Exon"] = exonTrx["Region"].str.replace(r"^exon_", "", regex=True).astype(int)
    exonTrx = exonTrx.sort_values(by=["Transcript", "Exon"]).reset_index(drop=True)
    exonTrx['Exon'] = exonTrx['Exon'].astype(str)

    geneTrx = (
        geneTrx[["Gene", "Transcript", "cdsStart", "cdsEnd", "strand"]]
                .drop_duplicates()
                .reset_index()
        )

    ids = (
        svCovDf["Info"]
        .str.split("|", expand=True)
    )
    ids.columns = ["RegionType", "Region", "Gene", "GeneId", "Transcript", "cdsStart", "cdsEnd", "strand"]
    svTrx = pd.concat([svCovDf, ids[["Region", "Transcript", "cdsStart", "cdsEnd", "strand"]]], axis=1)[
        ["Gene", "Transcript", "cdsStart", "cdsEnd", "strand"]
    ].drop_duplicates()
    svTrx["cdsStart"] = svTrx["cdsStart"].astype(int)
    svTrx["cdsEnd"] = svTrx["cdsEnd"].astype(int)

    # list of all genes and transcripts
    knownTrx = pd.concat(
        [svTrx[["Gene", "Transcript"]], geneTrx[["Gene", "Transcript"]]],
        axis=0,
        ignore_index=True,
    ).drop_duplicates()

    targetDf["Region"] = targetDf["Info"].str.split("|").str[1]

    #########################################
    #
    # Get small variants
    #
    #########################################

    print("Gathering gene variants...", file=sys.stderr)

    vcf_in = pysam.VariantFile(args.vcf_file, "r")
    vep = vepHeader(vcf_in.header)

    # Set variant db header fields if variant_db_tag is provided
    variant_db = None
    if args.variant_db_vcf_tag:
        variant_db = {}
        for hdr in vcf_in.header.records:
            if hdr.key == args.variant_db_vcf_tag:
                i = 0
                for j in hdr.value.split('|'):
                    variant_db[j] = i
                    i+=1
    
    # get variants
    variant_list = []

    for variant in vcf_in:
        vartype = ""
        if len(variant.ref) == len(variant.alts[0]):
            vartype = "SNV"
        else:
            vartype = "INDEL"

        filter_keys = list(variant.filter.keys())
        varfilter = "PASS" if not filter_keys else ";".join(filter_keys)

        # adding 241101, skip variants if components of MNV via MNV_tag
        if "MNVTAG" in variant.info:
            # get MNV_tag value
            mnv_tag = variant.info.get("MNVTAG")
            # format is chrom:pos_ref->alt
            mnv_alt = mnv_tag.split(">")[-1]
            # continue if this alt doesn't match the current alt
            if mnv_alt != variant.alts[0]:
                continue

        abundance = "N/A"
        totalReads = "N/A"
        variantReads = "N/A"
        abundance = round(variant.samples[0]["AF"][0] * 100, 2)
        totalReads = variant.samples[0]["DP"]
        variantReads = variant.samples[0]["AD"][1]

        # get VEP annotation
        csq = variant.info.get("CSQ")

        if csq is None:
            sys.exit("No VEP fields")

        gene = "N/A"
        gene_id = "N/A"
        transcript = "N/A"
        csyntax = "N/A"
        psyntax = "N/A"
        consequence = "N/A"
        exon = "N/A"
        intron = "N/A"
        popmaf = "N/A"
        customannotation = "N/A"
        in_custom_db = False

        for i in variant.info.get("CSQ", ""):
            csq = i.split("|")

            # skip blacklisted variants, if blacklist tag is provided and present in VEP annotation
            if args.blacklist_tag and vep.get(args.blacklist_tag) is not None:
                if csq[vep[args.blacklist_tag]]:
                    continue

            # get pop allele frequency. This is present for each transcript annotation, but is always the same
            if csq[vep["MAX_AF"]] != "":
                popmaf = float(csq[vep["MAX_AF"]])

            # check if this is in the list of transcripts. only variants annotated with a known transcript will be reported
            if geneTrx["Transcript"].str.contains(csq[vep["Feature"]]).any():
                transcript = csq[vep["Feature"]]
                gene = csq[vep["SYMBOL"]]
                gene_id = csq[vep["Gene"]]

                consequences = csq[vep["Consequence"]].split("&")
                
                if gene in nonCodingGeneList.keys() and any(c in nonCodingGeneList[gene] for c in consequences):
                    consequence = next(c for c in consequences if c in nonCodingGeneList[gene])

                elif any(c in nonSynon for c in consequences):
                    consequence = next(c for c in consequences if c in nonSynon)

                else:
                    continue

                csyntax = csq[vep["HGVSc"]].split(":")
                if len(csyntax) > 1:
                    csyntax = csyntax[1]
                else:
                    if (
                        csq[vep["STRAND"]]
                        != geneTrx[geneTrx["Transcript"] == transcript]["strand"].tolist()[
                            0
                        ]
                    ):
                        sys.exit("strands dont match")

                    if "upstream" in consequence:
                        if csq[vep["STRAND"]] == '1':
                            # *  ---->
                            distance = (
                                geneTrx[geneTrx["Transcript"] == transcript][
                                    "cdsStart"
                                ].min()
                                - variant.POS
                            )
                            csyntax = "c.-" + str(distance) + variant.REF + ">" + csq[0]
                        else:
                            # <---- *
                            distance = (
                                variant.POS
                                - geneTrx[geneTrx["Transcript"] == transcript][
                                    "cdsEnd"
                                ].max()
                            )
                            csyntax = (
                                "c.-"
                                + str(distance)
                                + revcomp(variant.REF)
                                + ">"
                                + revcomp(csq[0])
                            )

                    elif "downstream" in consequence:
                        if csq[vep["STRAND"]] == '1':
                            # ---->  *
                            distance = (
                                variant.POS
                                - geneTrx[geneTrx["Transcript"] == transcript][
                                    "cdsEnd"
                                ].max()
                            )
                            csyntax = (
                                "c.+"
                                + str(csq[vep["DISTANCE"]])
                                + variant.REF
                                + ">"
                                + csq[0]
                            )
                        else:
                            # *  <----
                            distance = (
                                geneTrx[geneTrx["Transcript"] == transcript]["cdsEnd"].min()
                                - variant.POS
                            )
                            csyntax = (
                                "c.+"
                                + str(csq[vep["DISTANCE"]])
                                + revcomp(variant.REF)
                                + ">"
                                + revcomp(csq[0])
                            )
                    elif 'deletion' in csq[0]: # this is a deletion that affects an exon
                        if csq[vep["cDNA_position"]].startswith('?'): # deletion is at start of exon
                            exonStart = exonTrx.loc[
                                            (exonTrx["Transcript"] == transcript) & 
                                            (exonTrx["Exon"].isin(csq[vep["EXON"]].split('/')[0].split('-'))),
                                            'exonStart'
                                        ].iloc[0]

                            delcDNAEnd = int(csq[vep["cDNA_position"]][2:])
                            delStart = variant.POS
                            delEnd = variant.INFO["END"]
                            
                            if csq[vep["STRAND"]] == '1': # deletion is like this -->intr[on/exo]n-->
                                delcDNAStart = delcDNAEnd - (delEnd - exonStart)
                                csyntax = f"c.{delcDNAStart}-{exonStart-delStart}_{delcDNAEnd}del"
                            else: # deletion is like this <--ex[on/intr]on<--
                                delcDNAStart = delcDNAEnd - (exonStart - delStart)
                                csyntax = f"c.{delcDNAStart}-{delEnd-exonStart}_{delcDNAEnd}del"

                        elif csq[vep["cDNA_position"]].endswith('?'): # deletion is at end of exon
                            exonEnd = exonTrx.loc[
                                            (exonTrx["Transcript"] == transcript) & 
                                            (exonTrx["Exon"].isin(csq[vep["EXON"]].split('/')[0].split('-'))), 
                                            'exonEnd'
                                        ].iloc[-1]

                            delcDNAStart = int(csq[vep["cDNA_position"]][:-2])
                            delStart = variant.POS
                            delEnd = variant.INFO["END"]
                            
                            if csq[vep["STRAND"]] == '1': # deletion is like this -->ex[on/intr]on-->
                                delcDNAEnd = delcDNAStart + (exonEnd - delStart)
                                csyntax = f"c.{delcDNAStart}_{delcDNAEnd}+{exonEnd-delStart}del"
                            else: # deletion is like this <--intr[on/ex]on<--
                                delcDNAEnd = delcDNAStart + (delEnd - exonEnd)
                                csyntax = f"c.{delcDNAStart}_{delcDNAEnd}+{exonEnd-delStart}del"

                        elif '-' in csq[vep["cDNA_position"]]:
                            delcDNAStart = int(csq[vep["cDNA_position"]].split('-')[0])
                            delcDNAEnd = int(csq[vep["cDNA_position"]].split('-')[1])
                            csyntax = f"c.{delcDNAStart}_{delcDNAEnd}del"

                        else :
                            csyntax = consequence

                    else:
                        csyntax = consequence

                psyntax = csq[vep["HGVSp"]].split(":")
                if len(psyntax) > 1:
                    psyntax = convert_aa(psyntax[1])
                    psyntax = re.sub("\%3D", "=", psyntax)

                elif csq[vep["Protein_position"]] and 'deletion' in csq[0]:                
                    if csq[vep["cDNA_position"]].startswith('?'): # deletion is at start of exon
                            exonStart = exonTrx.loc[
                                            (exonTrx["Transcript"] == transcript) & 
                                            (exonTrx["Exon"].isin(csq[vep["EXON"]].split('/')[0].split('-'))), 
                                            'exonStart'
                                        ].iloc[0]

                            delDNAEnd = int(csq[vep["cDNA_position"]][2:])
                            delCodonEnd = int(csq[vep["Protein_position"]][2:])
                            delStart = variant.POS
                            delEnd = variant.INFO["END"]
                            
                            if csq[vep["STRAND"]] == '1': # deletion is like this -->intr[on/exo]n-->
                                delCodonStart = delCodonEnd-ceiling((delEnd - exonStart)/3)
                                psyntax = f"p.({delCodonStart}_{delCodonEnd})del"
                            else: # deletion is like this <--ex[on/intr]on<--
                                delCodonStart = delCodonEnd-ceiling((exonStart-delStart)/3)
                                psyntax = f"p.({delCodonStart}_{delCodonEnd})del"

                    elif csq[vep["cDNA_position"]].endswith('?'): # deletion is at end of exon
                        exonEnd = exonTrx.loc[
                                        (exonTrx["Transcript"] == transcript) & 
                                        (exonTrx["Exon"].isin(csq[vep["EXON"]].split('/')[0].split('-'))), 
                                        'exonEnd'
                                    ].iloc[-1]

                        delcDNAStart = int(csq[vep["cDNA_position"]][:-2])
                        delCodonStart = int(csq[vep["Protein_position"]][:-2])
                        delStart = variant.POS
                        delEnd = variant.INFO["END"]
                        
                        if csq[vep["STRAND"]] == '1': # deletion is like this -->ex[on/intr]on-->
                            delCodonEnd = delcDNAStart + ceiling((exonEnd - delStart)/3)
                            psyntax = f"p.({delCodonStart}_{delCodonEnd})del"
                        else: # deletion is like this <--intr[on/ex]on<--
                            delCodonEnd = delcDNAStart + ceiling((delStart - exonEnd)/3)
                            psyntax = f"p.({delCodonStart}_{delCodonEnd})del"

                    elif '-' in csq[vep["Protein_position"]]:
                        delCodonStart = int(csq[vep["Protein_position"]].split('-')[0])
                        delCodonEnd = int(csq[vep["Protein_position"]].split('-')[1])
                        psyntax = f"p.({delCodonStart}_{delCodonEnd})del"

                    elif csq[vep["Protein_position"]]:
                        psyntax = "p.?"

                    else:
                        psyntax = consequence

                elif csq[vep["Protein_position"]]:
                    psyntax = "p.?"

                else:
                    psyntax = consequence

                impact = csq[vep["IMPACT"]]
                exon = csq[vep["EXON"]] or "N/A"
                intron = csq[vep["INTRON"]] or "N/A"
                customannotation = [csq[vep["Existing_variation"]]] or []
                priorvariants = None

                if args.custom_annotation_tag and vep.get(args.custom_annotation_tag) is not None:
                    if csq[vep[args.custom_annotation_tag]]:
                        customannotation.append(f"{args.custom_annotation_tag}=" + str(csq[vep[args.custom_annotation_tag]] or 0))
                        in_custom_db = True

                # check for forced genotype tag and update annotation field.
                if 'FGT' in variant.info:
                    customannotation.append(f'FORCEGT={variant.id}')

                if variant_db and variant.info.get(args.variant_db_vcf_tag) is not None:
                    dbrecords = []
                    for v in variant.info.get(args.variant_db_vcf_tag):
                        pvnfo = dict(zip(variant_db.keys(),v.split('|')))
                        dbrecords.append('|'.join([pvnfo['accession'],pvnfo['date'],pvnfo['filter'],str(round(float(pvnfo['vaf'])*100,2))+'%']))
                        
                    if len(dbrecords) > 0:
                        priorvariants = ','.join(dbrecords)


        # convert pop maf to percent
        if popmaf != "N/A":
            popmaf = round(float(popmaf) * 100, 3)

        # categories
        category = "None"
        if varfilter != "PASS":
            category = "FILTERED"

        elif (
            popmaf != "N/A"
            and popmaf >= float(args.pop_af)
            and not in_custom_db
        ):
            category = "SNP"

        else:
            category = "PASS"

        # only include all variants <=0.1% and ns or specific noncoding variants
        if consequence != "N/A":
            # Set popmaf NA to None for JSON output of null to conform to spec
            if popmaf == "N/A":
                popmaf = None

            variant_list.append(dict(zip(variants.columns, [
                    category,
                    vartype,
                    varfilter,
                    str(variant.contig),
                    variant.pos,
                    variant.ref,
                    variant.alts[0],
                    gene,
                    gene_id,
                    transcript,
                    consequence,
                    csyntax,
                    psyntax,
                    exon,
                    intron,
                    popmaf,
                    ';'.join(customannotation),
                    totalReads,
                    variantReads,
                    abundance,
                    priorvariants
                ])))

    if variant_list:
        variants = pd.DataFrame(variant_list)

        # print variants to tab-delimited file or stdout if no file is given
        if args.outfile:
            variants.to_csv(args.outfile, sep="\t", index=False)
        else:
            print(variants.to_csv(sep="\t", index=False))
    else:
        print("No variants found that meet criteria.", file=sys.stderr)

    





if __name__ == "__main__":
    main()