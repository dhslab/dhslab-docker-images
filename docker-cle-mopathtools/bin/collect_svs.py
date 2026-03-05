#!/usr/bin/env python

import argparse
import json
import os
import re
import sys

import numpy as np
import pandas as pd
import pysam
import natsort

__version__ = "1.0.0"

ACROCENTRICS = ["chr13", "chr14", "chr15", "chr21", "chr22"]

def checkfile(file_path):
    """Check if a file exists at the given path."""
    if not os.path.exists(file_path):
        raise argparse.ArgumentTypeError(f"The file {file_path} does not exist.")
    return file_path

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

def vep_csq_to_dataframe(csq_tag: str | tuple | None, header_description: str) -> pd.DataFrame:
    """
    Parses a VEP CSQ string from a VCF INFO tag into a DataFrame.

    Args:
        csq_tag: The CSQ string from the VCF INFO field. Can be a string or a tuple of strings.
        header_description: The VCF header description for the CSQ field. 

    Returns:
        A DataFrame of VEP annotations.
    """
    fields = _get_fields_from_header(header_description)
    
    # Return an empty DataFrame with the correct columns if the tag is missing
    if not csq_tag or not fields:
        return pd.DataFrame(columns=fields + ["KnownGene", "KnownTrx"])

    # Pysam may return a tuple of strings for CSQ, whereas cyvcf2 returns a single comma-separated string.
    # This handles both cases.
    if isinstance(csq_tag, str):
        csq_entries = csq_tag.split(',')
    else:
        csq_entries = csq_tag

    # 1. Create a list of dictionaries from the CSQ string
    annotations = [dict(zip(fields, entry.split("|"))) for entry in csq_entries]
    
    # 2. Convert the list to a DataFrame
    df = pd.DataFrame(annotations, columns=fields)

    # --- 3. Data Cleaning and Type Conversion ---

    # remove leading or trailing '.' from Allele column 
    if "Allele" in df.columns:
        df["Allele"] = df["Allele"].str.strip('.')
    
    # If a SYMBOL is missing, use the Gene ID instead
    if "SYMBOL" in df.columns and "Gene" in df.columns:
        df["SYMBOL"] = df.apply(lambda row: row['Gene'] if pd.isna(row['SYMBOL']) or row['SYMBOL'] == '' else row['SYMBOL'], axis=1)
    
    # Standardize STRAND representation
    if "STRAND" in df.columns:
        df["STRAND"] = pd.to_numeric(df["STRAND"], errors='coerce').fillna(0).astype(int)

    # if TranscriptCoordinates (START and END) are in columns
    if "START" in df.columns and "END" in df.columns:
        # make START and END numeric
        df["START"] = pd.to_numeric(df["START"], errors='coerce')
        df["END"] = pd.to_numeric(df["END"], errors='coerce')
    
    else:
        df['START'] = None
        df['END'] = None

    # Convert DISTANCE to a numeric type, defaulting to 0 for empty strings
    if "DISTANCE" in df.columns:
        df["DISTANCE"] = pd.to_numeric(df["DISTANCE"], errors='coerce').fillna(0).astype(int)

    else:
        df['DISTANCE'] = None

    # Convert PICK to a binary indicator
    if "PICK" in df.columns:
        df["PICK"] = df["PICK"].apply(lambda x: 1 if x == "1" else 0)

    else:
        df['PICK'] = None
    
    # Add custom boolean flags
    df["KnownGene"] = 0
    df["KnownTrx"] = 0

    return df

def known_sv_genes_tag_to_dataframe(variant, header_description: str) -> pd.DataFrame:
    """
    Parses the KnownSvGenes INFO tag into a DataFrame with VEP-like columns.

    Args:
        variant: A pysam variant object with pos, stop, info attributes.
        header_description: The VCF header description for the KnownSvGenes field.

    Returns:
        A DataFrame of known gene overlaps with calculated consequences.
    """
    fields = _get_fields_from_header(header_description)
    
    # Safely get the tag, return empty df if not present in header or record
    if "KnownSvGenes" not in variant.info:
        return pd.DataFrame(columns=fields + ["DISTANCE", "Consequence","KnownGene","KnownTrx"])
    known_sv_tag = variant.info["KnownSvGenes"]
    
    # Return an empty DataFrame if the tag is missing or variant positions are invalid
    pos1 = variant.pos
    pos2 = variant.stop
    if not known_sv_tag or not fields or pos1 is None or pos2 is None:
        return pd.DataFrame(columns=fields + ["DISTANCE", "Consequence","KnownGene","KnownTrx"])

    # Pysam may return a tuple of strings, whereas cyvcf2 returns a single comma-separated string.
    if isinstance(known_sv_tag, str):
        known_sv_entries = known_sv_tag.split(',')
    else:
        known_sv_entries = known_sv_tag

    # 1. Create a list of dictionaries from the KnownSvGenes string
    parsed_genes = [dict(zip(fields, entry.split("|"))) for entry in known_sv_entries]

    # 2. Convert the list to a DataFrame
    df = pd.DataFrame(parsed_genes, columns=fields)

    # --- 3. Data Cleaning and Type Conversion ---
    # Convert position columns to numeric types
    if "START" in df.columns:
        df["START"] = pd.to_numeric(df["START"], errors='coerce')
    if "END" in df.columns:
        df["END"] = pd.to_numeric(df["END"], errors='coerce')
    
    # Standardize STRAND representation
    if "STRAND" in df.columns:
        df["STRAND"] = pd.to_numeric(df["STRAND"], errors='coerce').fillna(0).astype(int)
 
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
            # if variant.info.get("SVTYPE") == "DEL" return deletion if DUP return transcript_amplification
            if variant.info.get("SVTYPE") == "DEL":
                return "deletion"
            else:
                return "transcript_amplification"
        if (row["STRAND"] == "+" and row["START"] > pos2) or (row["STRAND"] == "-" and row["END"] < pos1):
            return "upstream_gene_variant"
        return "downstream_gene_variant"
    df["Consequence"] = df.apply(get_consequence, axis=1)

    return df

def vepGeneEffect(row):
    exon = row.get("EXON")
    intron = row.get("INTRON")
    dist = row.get("DISTANCE")
    cons = row.get("Consequence")

    if isinstance(exon, str) and "/" in exon:
        return f"exon{exon.split('/')[0]}"
    elif isinstance(intron, str) and "/" in intron:
        return f"intron{intron.split('/')[0]}"
    elif pd.notna(dist) and dist > 0:
        if isinstance(cons, str) and "upstream" in cons:
            return f"{int(dist)}bp upstream"
        elif isinstance(cons, str) and "downstream" in cons:
            return f"{int(dist)}bp downstream"
        return ""
    return "intragenic"

def _tail_or_none(value):
    if isinstance(value, str) and value:
        return value.split("-")[-1]
    return None


def get_gene_syntax(row, bnd_orientation, chr_l, chr_r):
    """Calculates gene-based syntax (genestring and genedetail) for a BND annotation row."""
    gene_l = row.get("SYMBOL_l") or "INTERGENIC"
    region_l = "N/A" if gene_l == "INTERGENIC" else f"{row.get('Feature_l')}:{row.get('GeneEffect_l')}"
    strand_l = row.get("STRAND_l") or 0

    gene_r = row.get("SYMBOL_r") or "INTERGENIC"
    region_r = "N/A" if gene_r == "INTERGENIC" else f"{row.get('Feature_r')}:{row.get('GeneEffect_r')}"
    strand_r = row.get("STRAND_r") or 0

    require_strand = row.get("KNOWNSVREQUIRESTRAND", 1)
    inframe_splice = row.get("InframeSplice", 0)

    # Calculate fusion orientation: 1 is same strand and -1 is opposite.
    fusion_orientation = strand_l * strand_r * bnd_orientation

    chr_l_num_str = chr_l.replace('chr', '').replace('X', '23').replace('Y', '24').replace('M', '25')
    chr_r_num_str = chr_r.replace('chr', '').replace('X', '23').replace('Y', '24').replace('M', '25')

    genestring = "N/A"
    genedetail = "N/A"

    # This is an in-frame fusion in the correct orientation for transcription
    if fusion_orientation >= 0:
        gene_separator = "(::)" if inframe_splice or not require_strand else "(--)"
        gene_separator_detail = "(::)" if inframe_splice else "(--)"

        if strand_l == 1:  # FWD strand, left to right
            genestring = f"{gene_l}{gene_separator}{gene_r}"
            genedetail = f"{gene_l}(+)({region_l}){gene_separator_detail}{gene_r}(+)({region_r})"

        elif strand_l == -1:  # REV strand, right to left
            genestring = f"{gene_r}{gene_separator}{gene_l}"
            genedetail = f"{gene_r}(-)({region_r}){gene_separator_detail}{gene_l}(-)({region_l})"

    # Genes are not in the proper orientation for a fusion
    else:
        gene_separator = "//"
        
        # Order genes by chromosome number
        if chr_l in ['chrX', 'chrY', 'chrM'] or int(chr_l_num_str) < int(chr_r_num_str):
            genestring = f"{gene_l}{gene_separator}{gene_r}"
            if gene_l == "INTERGENIC":
                genedetail = f"{gene_l}({region_l}){gene_separator}{gene_r}({'+' if strand_r == 1 else '-'})({region_r})"
            elif gene_r == "INTERGENIC":
                genedetail = f"{gene_l}({'+' if strand_l == 1 else '-'})({region_l}){gene_separator}{gene_r}({region_r})"
            else:
                genedetail = f"{gene_l}({'+' if strand_l == 1 else '-'})({region_l}){gene_separator}{gene_r}({'+' if strand_r == 1 else '-'})({region_r})"
        else:
            genestring = f"{gene_r}{gene_separator}{gene_l}"
            if gene_l == "INTERGENIC":
                genedetail = f"{gene_r}({'+' if strand_r == 1 else '-'})({region_r}){gene_separator}{gene_l}({region_l})"
            elif gene_r == "INTERGENIC":
                genedetail = f"{gene_r}({region_r}){gene_separator}{gene_l}({'+' if strand_l == 1 else '-'})({region_l})"
            else:
                genedetail = f"{gene_r}({'+' if strand_r == 1 else '-'})({region_r}){gene_separator}{gene_l}({'+' if strand_l == 1 else '-'})({region_l})"

    return genestring, genedetail

def main():
    parser = argparse.ArgumentParser(description="Collect Structural Variants from ChromoSeq output.")
    parser.add_argument("--sv_vcf", required=True, type=checkfile, help="Annotated SV VCF file (*.sv_annotated.vcf.gz)")
    parser.add_argument("--sv_targets", required=True, type=checkfile, help="Recurrent SV target list CSV")
    parser.add_argument("--bed_file", required=True, type=checkfile, help="Mopath formatted bedfile or coverage report")
    parser.add_argument("--outfile", type=str, help="Output tab-delimited file. If not provided, output is written to stdout.")
    parser.add_argument("--version", "-v", action="version", version="%(prog)s: " + __version__)
    args = parser.parse_args()

    # Load accessory files
    nonSynon = {
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "frameshift_variant",
        "stop_lost",
        "start_lost",
        "transcript_ablation",
        "transcript_amplification",
        "inframe_insertion",
        "inframe_deletion",
        "missense_variant",
        "protein_altering_variant",
        "feature_truncation",
        "coding_sequence_variant",
    }

    recurrentSvs = pd.read_csv(args.sv_targets, sep=",", header=None)
    swapped = recurrentSvs.copy()
    swapped.iloc[:, [0, 1]] = swapped.iloc[:, [1, 0]].values
    recurrentSvs = pd.concat([recurrentSvs, swapped], axis=0, ignore_index=True)
    recurrentSvs.columns = ["KNOWNSVGENE1", "KNOWNSVGENE2", "KNOWNSVREQUIRESTRAND", "KNOWNSVTYPE"]
    recurrentSvs = recurrentSvs.where(pd.notna(recurrentSvs), None)

    covDf = pd.read_csv(
        args.bed_file,
        header=None,
        skiprows=1,
        names=["Chromosome", "Start", "End", "Gene", "Info"],
        sep="\t",
    )
    geneCovDf = covDf[covDf["Info"].str.contains("exon")]
    svCovDf = covDf[covDf["Info"].str.contains(r"transcript|gene", regex=True)]

    ids = (
        geneCovDf["Info"]
        .str.split("|", expand=True)
        .replace("\s\+\s\S+", "", regex=True)
        .loc[:, 2:]
    )
    ids.columns = ["Transcript", "Region", "cdsStart", "cdsEnd", "strand"]
    geneTrx = (
        pd.concat([geneCovDf, ids], axis=1)
        .drop_duplicates()
        .reset_index()
    )
    geneTrx = (
        geneTrx[["Gene", "Transcript", "cdsStart", "cdsEnd", "strand"]]
                .drop_duplicates()
                .reset_index(drop=True)
        )

    ids = (
        svCovDf["Info"]
        .str.split("|", expand=True)
        .replace("\s\+\s\S+", "", regex=True)
        .loc[:, 2:]
    )
    ids.columns = ["Transcript", "Region", "cdsStart", "cdsEnd", "strand"]
    svTrx = pd.concat([svCovDf, ids], axis=1)[
        ["Gene", "Transcript", "cdsStart", "cdsEnd", "strand"]
    ].drop_duplicates()

    knownTrx = pd.concat(
        [svTrx[["Gene", "Transcript"]], geneTrx[["Gene", "Transcript"]]],
        axis=0,
        ignore_index=True,
    ).drop_duplicates()

    reportableCnvGeneList = geneTrx["Gene"].drop_duplicates().tolist()

    # Initialize dataframe for all SVs
    svs_df_columns = [
        "category", "type", "chrom1", "pos1", "chrom2", "pos2", "length",
        "csyntax", "psyntax", "bands", "known_genes", "known_gene_detail",
        "total_genes", "filters", "id", "abundance", "info",
    ]
    svs = pd.DataFrame(columns=svs_df_columns)

    # Process SVs
    print("Gathering SVs...", file=sys.stderr)
    svvcf = pysam.VariantFile(args.sv_vcf)
    csq_header_desc = svvcf.header.info['CSQ'].description
    known_sv_genes_header_desc = ""
    if 'KnownSvGenes' in svvcf.header.info:
        known_sv_genes_header_desc = svvcf.header.info['KnownSvGenes'].description

    passedvars = {}
    sv_deldup_list = []

    for variant in svvcf:
        
        if variant.info["SVTYPE"] in ["BND", "INV"]:
            passedvars[variant.id] = variant
        
        else:
            # Hard filter DEL/DUPs that do not involve genes
            if "KnownSvGenes" not in variant.info and "CSQ" not in variant.info:
                continue

            filter_keys = list(variant.filter.keys())
            filter = "PASS" if not filter_keys else ";".join(filter_keys)

            vartype = variant.info.get("SVTYPE")
            chr1 = str(variant.chrom)
            pos1 = variant.pos
            svlen = None
            if variant.info.get("SVLEN") is not None:
                svlen = int(abs(variant.info.get("SVLEN")[0]))

            pos2 = variant.stop if variant.stop is not None else "N/A"
            if pos2 == "N/A" and svlen != None:
                pos2 = pos1 + svlen
            
            # Get abundance information
            abundance = 0.0
            PR = (0, 0)
            SR = (0, 0)
            if "SR" in variant.samples[0] and variant.samples[0]["SR"][0] is not None:
                SR = variant.samples[0]["SR"]
            if "PR" in variant.samples[0] and variant.samples[0]["PR"][0] is not None:
                PR = variant.samples[0]["PR"]

            denominator = PR[0] + PR[1] + SR[0] + SR[1]
            abundance = round((SR[1] + PR[1]) / denominator * 100, 1) if denominator > 0 else 0.0

            info_dict = { "PR_READS": f"{PR[1]}/{PR[0] + PR[1]}",
                          "SR_READS": f"{SR[1]}/{SR[0] + SR[1]}",
                          "CONTIG": str(variant.info.get("CONTIG")) }

            csyntax = ""
            psyntax = ""
            bandstring = "None"
            genestring = "None"
            genedetail = "None"
            total_genes = "None"

            bands = []
            if "Cytobands" in variant.info and variant.info.get("Cytobands") is not None:
                bands = variant.info.get("Cytobands")
                bands = [b.replace("acen_", "") for b in bands]
            
            bandstring = bands[0] if bands else "None"
            if len(bands) > 1:
                bandstring = bands[0] + bands[-1]

            vepCsq = vep_csq_to_dataframe(variant.info.get("CSQ"), csq_header_desc)
            vepCsq["KnownGene"] = vepCsq["SYMBOL"].isin(knownTrx["Gene"]).astype(int)
            vepCsq.loc[vepCsq["SYMBOL"].isin(reportableCnvGeneList), "KnownGene"] = 1
            vepCsq["KnownTrx"] = vepCsq["Feature"].isin(knownTrx["Transcript"]).astype(int)
            vepCsq = vepCsq.sort_values(
                by=["SYMBOL", "KnownTrx", "PICK"], ascending=[True, False, False]
            ).drop_duplicates(subset="SYMBOL", keep="first")
        
            knownSvGeneDf = known_sv_genes_tag_to_dataframe(variant, known_sv_genes_header_desc)
            if not knownSvGeneDf.empty:                
                vepCsq = pd.concat([vepCsq, knownSvGeneDf[~knownSvGeneDf['SYMBOL'].isin(vepCsq['SYMBOL'])]], ignore_index=True)
                vepCsq = vepCsq.where(pd.notna(vepCsq), None)

            if not vepCsq.empty:
                total_genes = f"{len(vepCsq['SYMBOL'].unique())} genes"
                vepCsq = vepCsq.sort_values(by=["START"], ascending=[True])
                vepCsq["GeneEffect"] = vepCsq.apply(lambda r: vepGeneEffect(r), axis=1)
                vepCsq["GeneImpact"] = vepCsq["Consequence"].apply(
                    lambda r: int(len(set(r.split("&")) & set(nonSynon)) > 0) if isinstance(r, str) and r else None
                )

            # Format DEL/DUP events
            if vartype == "DEL":
                csyntax = chr1 + ":g." + str(pos1) + "_" + str(pos2) + "del"
                if (bands and bands[0].find("p") > -1 and bands[-1].find("q") > -1):
                    psyntax = "seq[GRCh38] -" + chr1.replace("chr", "")
                elif (bands and any("q11" in b for b in bands) and any("qter" in b for b in bands) and chr1 in ACROCENTRICS):
                    psyntax = "seq[GRCh38] -" + chr1.replace("chr", "")
                elif (bands and bands[0].find("p") > -1):
                    psyntax = ("seq[GRCh38] del(" + chr1.replace("chr", "") + ")(" + bands[0] + bands[-1] + ")")
                elif bands:
                    psyntax = ("seq[GRCh38] del(" + chr1.replace("chr", "") + ")(" + bands[0] + bands[-1] + ")")
            elif vartype == "DUP" or vartype == "INS":
                csyntax = chr1 + ":g." + str(pos1) + "_" + str(pos2) + vartype.lower()
                if (bands and bands[0].find("p") > -1 and bands[-1].find("q") > -1):
                    psyntax = "seq[GRCh38] +" + chr1.replace("chr", "")
                elif (bands and any("q11" in b for b in bands) and any("qter" in b for b in bands) and chr1 in ACROCENTRICS):
                    psyntax = "seq[GRCh38] +" + chr1.replace("chr", "")
                elif (bands and bands[0].find("p") > -1):
                    psyntax = ("seq[GRCh38] " + vartype.lower() + "(" + chr1.replace("chr", "") + ")(" + bands[0] + bands[-1] + ")")
                elif bands:
                    psyntax = ("seq[GRCh38] " + vartype.lower() + "(" + chr1.replace("chr", "") + ")(" + bands[0] + bands[-1] + ")")

            knownGeneDf = vepCsq[(vepCsq["KnownTrx"] == 1) & (vepCsq["DISTANCE"] == 0)].sort_values(by=["START"])[["SYMBOL", "GeneImpact", "GeneEffect"]]
            
            if len(knownGeneDf) > 15:
                genestring = f"{len(knownGeneDf)} genes"
                genedetail = genestring
            elif len(knownGeneDf) > 0:
                genestring = ",".join(knownGeneDf["SYMBOL"].unique().tolist())
                genedetail = (
                    vartype.lower() + "[" + ",".join(
                        knownGeneDf.apply(
                            lambda x: f"{x['SYMBOL']}({x['GeneEffect']})" if x['GeneEffect']!="intragenic" else f"{x['SYMBOL']}", axis=1
                        ).tolist()
                    ) + "]"
                )

            isRecurrentSv = ((recurrentSvs["KNOWNSVGENE1"].isin(knownGeneDf['SYMBOL']))
                            & (recurrentSvs["KNOWNSVGENE2"].isin(knownGeneDf['SYMBOL']))
                            & (recurrentSvs["KNOWNSVREQUIRESTRAND"] == 1)
                            & (recurrentSvs["KNOWNSVTYPE"].isin([vartype,"*"]))).any()
            
            if isRecurrentSv:
                filter_list = [f for f in filter.split(";") if f not in ["MaxDepth", "MinSomaticScore"]]
                filter = "PASS" if not filter_list else ';'.join(filter_list)

            category = ""
            if filter == "PASS" and isRecurrentSv:
                category = "CNV"
            elif filter != "PASS" and isRecurrentSv:
                category = "OTHERSV"
            elif filter == "PASS" and genestring != "None":
                category = "OTHERSV"
            else:
                category = None

            # Add DEL/DUP info to list
            sv_deldup_list.append(dict(zip(svs_df_columns, [
                    category, vartype, chr1, int(pos1), chr1, int(pos2), svlen,
                    csyntax, psyntax, bandstring, genestring, genedetail, total_genes,
                    filter, str(variant.id), abundance, ';'.join([f"{k}={v}" for k, v in info_dict.items()]),
                ])))
            
            # Get fusion genes formed by a deletion if the deletion juxtaposes 2 genes
            if vartype == "DEL" and len(vepCsq[(vepCsq['START']<=pos1)]) > 1 and len(vepCsq[(vepCsq['END']>=pos2)]) > 1:
                # get cross product of left and right candidate fusions (those that extend beyond the DEL interval)
                candidate_fusions = pd.merge(vepCsq[(vepCsq['START']<=pos1)], 
                                             vepCsq[(vepCsq['END']>=pos2)], how='cross', suffixes=('_l', '_r'))

                is_inframe_splice = (
                    (candidate_fusions['IntronFrame_l'] != -1) &
                    (candidate_fusions['IntronFrame_l'] == candidate_fusions['IntronFrame_r']) &
                    (candidate_fusions['INTRON_l'].notna()) & (candidate_fusions['INTRON_l'] != '') &
                    (candidate_fusions['INTRON_r'].notna()) & (candidate_fusions['INTRON_r'] != '')
                )
                candidate_fusions['InframeSplice'] = np.where(is_inframe_splice, 1, 0)

                candidate_fusions = (candidate_fusions.sort_values(by=["KnownGene_l", "KnownGene_r", "InframeSplice", "KnownTrx_l", "KnownTrx_r", "PICK_l", "PICK_r"], 
                                                                    ascending=[False, False, False, False, False, False, False])
                                    .drop_duplicates(subset=["SYMBOL_l","SYMBOL_r"], keep="first")
                                    .query("SYMBOL_l != 'INTERGENIC' and SYMBOL_r != 'INTERGENIC'")
                                    .reset_index(drop=True))

                candidate_fusions = pd.merge(candidate_fusions,recurrentSvs,how='left', left_on=['SYMBOL_l','SYMBOL_r'], right_on=['KNOWNSVGENE1','KNOWNSVGENE2'])
                candidate_fusions['RecurrentSV'] = candidate_fusions.apply(lambda r: 1 if r['KNOWNSVGENE1'] is not None and (r['KNOWNSVTYPE'] is None or r['KNOWNSVTYPE'] == vartype) else 0, axis=1)

                candidate_fusions = candidate_fusions[(candidate_fusions['RecurrentSV'] == 1) |
                                                        ((candidate_fusions['SYMBOL_l'] != candidate_fusions['SYMBOL_r']) & 
                                                        ((candidate_fusions['BIOTYPE_l']=='protein_coding') | (candidate_fusions['BIOTYPE_r']=='protein_coding')))]

                if not candidate_fusions.empty:
                    candidate_fusions['EXON_l'] = candidate_fusions['EXON_l'].apply(_tail_or_none)
                    candidate_fusions['INTRON_l'] = candidate_fusions['INTRON_l'].apply(_tail_or_none)
                    candidate_fusions['EXON_r'] = candidate_fusions['EXON_r'].apply(_tail_or_none)
                    candidate_fusions['INTRON_r'] = candidate_fusions['INTRON_r'].apply(_tail_or_none)
                    candidate_fusions['SYMBOL_l'] = candidate_fusions['SYMBOL_l'].fillna('INTERGENIC')
                    candidate_fusions['SYMBOL_r'] = candidate_fusions['SYMBOL_r'].fillna('INTERGENIC')

                    candidate_fusions[['genestring', 'genedetail']] = candidate_fusions.apply(
                        lambda row: get_gene_syntax(row, 1, chr1, chr1),
                        axis=1,
                        result_type='expand'
                    )

                    candidate_fusions = candidate_fusions.sort_values(by=['RecurrentSV', 'InframeSplice'],ascending=[False,False])

                    category = ""
                    if filter == "PASS" and candidate_fusions.iloc[0]['RecurrentSV']:
                        category = "CNV"
                    elif filter != "PASS" and candidate_fusions.iloc[0]['RecurrentSV']:
                        category = "OTHERSV"
                    else:
                        category = None

                    if len(candidate_fusions) > 1:
                        info_dict['OTHERANNOTATIONS'] = '|'.join(candidate_fusions['genedetail'].tolist()[1:])

                    sv_deldup_list.append(dict(zip(svs_df_columns, [
                        category, vartype, chr1, int(pos1), chr1, int(pos2), svlen,
                        csyntax, psyntax, bandstring, candidate_fusions.iloc[0]['genestring'], candidate_fusions.iloc[0]['genedetail'], '2 genes',
                        filter, str(variant.id), abundance, ';'.join([f"{k}={v}" for k, v in info_dict.items()]),
                    ])))

    # Concatenate only if there are records to add
    if sv_deldup_list:
        svs = pd.concat([svs.dropna(axis=1, how='all'), pd.DataFrame(sv_deldup_list)], ignore_index=True)

    # Handle BNDs
    alreadydone = set()
    sv_bnd_list = []
    for v_id, variant in passedvars.items():
        
        mate_id_tuple = variant.info.get("MATEID")
        if not mate_id_tuple:
            continue

        mate_id = mate_id_tuple[0]
        
        if variant.id in alreadydone or mate_id in alreadydone or mate_id not in passedvars:
            continue
        
        mate = passedvars[mate_id]

        if ("KnownSvGenes" not in variant.info and "KnownSvGenes" not in mate.info) and \
           ("CSQ" not in variant.info or "CSQ" not in mate.info):
            continue

        # swap variant and mate if variant is downstream of mate
        if variant.alts[0].find("[") == 0 or variant.alts[0].find("]") == 0:
            variant, mate = mate, variant

        bnd_orientation = -1 if variant.alts[0].find("[") == 0 or variant.alts[0].find("]") > 0 else 1

        vartype = variant.info.get("SVTYPE")
        filter_keys = list(variant.filter.keys())
        filter = "PASS" if not filter_keys else ";".join(filter_keys)

        chr_l, pos_l = str(variant.chrom), variant.pos
        chr_r, pos_r = str(mate.chrom), mate.pos
        svlen = None
        genestring = "N/A"
        genedetail = "N/A"
        total_genes = "N/A"
        bandstring = "N/A"

        bands_l = "N/A"
        if "Cytobands" in variant.info and variant.info.get("Cytobands") is not None:
            bands_l = variant.info.get("Cytobands")[0].replace("acen_", "")

        bands_r = "N/A"
        if "Cytobands" in mate.info and mate.info.get("Cytobands") is not None:
            bands_r = mate.info.get("Cytobands")[0].replace("acen_", "")

                # make ISCN syntax
        svtype = None
        chr_l_num_str = chr_l.replace('chr', '').replace('X', '-1').replace('Y', '-2').replace('M', '-3')
        chr_r_num_str = chr_r.replace('chr', '').replace('X', '-1').replace('Y', '-2').replace('M', '3')

        if vartype == "INV":
            svtype = "inv"
            svlen = int(abs(pos_l - pos_r))
            csyntax = f"{chr_l}:g.{pos_l}(+)::{chr_r}:g.{pos_r}({'+' if bnd_orientation == 1 else '-'})"
            psyntax = f"seq[GRCh38] {svtype}({chr_l.replace('chr', '')})({bands_l};{bands_r})"
            bandstring = f"{bands_l};{bands_r}"

        elif vartype == "BND":
            svtype = "t"
            if chr_l in ['chrX', 'chrY', 'chrM'] or int(chr_l_num_str) < int(chr_r_num_str):
                csyntax = f"{chr_l}:g.{pos_l}(+)::{chr_r}:g.{pos_r}({'+' if bnd_orientation == 1 else '-'})"
                psyntax = f"seq[GRCh38] {svtype}({chr_l.replace('chr', '')};{chr_r.replace('chr', '')})({bands_l};{bands_r})"
                bandstring = f"{bands_l};{bands_r}"
            else:
                csyntax = f"{chr_r}:g.{pos_r}(+)::{chr_l}:g.{pos_l}({'+' if bnd_orientation == 1 else '-'})"
                psyntax = f"seq[GRCh38] {svtype}({chr_r.replace('chr', '')};{chr_l.replace('chr', '')})({bands_r};{bands_l})"
                bandstring = f"{bands_r};{bands_l}"

        # Calculate abundance from read counts
        PR = (0, 0)
        SR = (0, 0)
        if "SR" in variant.samples[0] and variant.samples[0]["SR"][0] is not None:
            SR = variant.samples[0]["SR"]
        if "PR" in variant.samples[0] and variant.samples[0]["PR"][0] is not None:
            PR = variant.samples[0]["PR"]
        
        if "DUX" in variant.id:
            SR = (0, SR[1])
            PR = (0, PR[1])

        supporting_reads = SR[1] + PR[1]
        total_reads = PR[0] + PR[1] + SR[0] + SR[1]
        abundance = round(supporting_reads / total_reads * 100, 1) if total_reads > 0 else 0.0

        info_dict = { "PR_READS": f"{PR[1]}/{PR[0] + PR[1]}",
                      "SR_READS": f"{SR[1]}/{SR[0] + SR[1]}",
                      "CONTIG": str(variant.info.get("CONTIG")) }
        
        vepCsq1 = vep_csq_to_dataframe(variant.info.get("CSQ"), csq_header_desc)
        vepCsq1 = vepCsq1[(vepCsq1['Allele']==variant.ref) | (vepCsq1['Allele']=="BND")].reset_index(drop=True)
        if not vepCsq1.empty:
            vepCsq1["KnownTrx"] = vepCsq1["Feature"].isin(knownTrx["Transcript"]).astype(int)
            vepCsq1["KnownGene"] = vepCsq1["SYMBOL"].isin(knownTrx["Gene"]).astype(int)

            if "IntronFrame" not in vepCsq1.columns:
                vepCsq1["IntronFrame"] = None
                    
            knownSvGene1Df = known_sv_genes_tag_to_dataframe(variant, known_sv_genes_header_desc)
            knownSvGene1Df['IntronFrame'] = 0
            vepCsq1 = pd.concat([vepCsq1, knownSvGene1Df[~knownSvGene1Df['SYMBOL'].isin(vepCsq1['SYMBOL'])]], ignore_index=True)
            vepCsq1 = vepCsq1.where(pd.notna(vepCsq1), None)

            vepCsq1["GeneEffect"] = vepCsq1.apply(lambda r: vepGeneEffect(r), axis=1)

        else:
            vepCsq1 = pd.DataFrame([{}], columns=vepCsq1.columns)

        vepCsq2 = vep_csq_to_dataframe(mate.info.get("CSQ"), csq_header_desc)
        vepCsq2 = vepCsq2[(vepCsq2['Allele']==mate.ref) | (vepCsq2['Allele']=="BND")].reset_index(drop=True)
        if not vepCsq2.empty:
            vepCsq2["KnownTrx"] = vepCsq2["Feature"].isin(knownTrx["Transcript"]).astype(int)
            vepCsq2["KnownGene"] = vepCsq2["SYMBOL"].isin(knownTrx["Gene"]).astype(int)

            if "IntronFrame" not in vepCsq2.columns:
                vepCsq2["IntronFrame"] = None

            knownSvGene2Df = known_sv_genes_tag_to_dataframe(mate, known_sv_genes_header_desc)
            knownSvGene2Df['IntronFrame'] = 0
            vepCsq2 = pd.concat([vepCsq2, knownSvGene2Df[~knownSvGene2Df['SYMBOL'].isin(vepCsq2['SYMBOL'])]], ignore_index=True)
            vepCsq2 = vepCsq2.where(pd.notna(vepCsq2), None)

            vepCsq2["GeneEffect"] = vepCsq2.apply(lambda r: vepGeneEffect(r), axis=1)

        else:
            vepCsq2 = pd.DataFrame([{}], columns=vepCsq2.columns)

        # cross vepCsq1 (left) with vepCsq2 (right) for all possible combinations
        bnd_annot = pd.merge(vepCsq1, vepCsq2, how='cross', suffixes=('_l', '_r'))
        # fillna SYMBOL w/ INTERGENIC and the rest with 0
        bnd_annot[['SYMBOL_l','SYMBOL_r']] = bnd_annot[['SYMBOL_l','SYMBOL_r']].fillna('INTERGENIC')

        if not bnd_annot.empty:
            is_inframe_splice = (
                (bnd_annot['IntronFrame_l'] != -1) &
                (bnd_annot['IntronFrame_l'] == bnd_annot['IntronFrame_r']) &
                (bnd_annot['INTRON_l'].notna()) & (bnd_annot['INTRON_l'] != '') &
                (bnd_annot['INTRON_r'].notna()) & (bnd_annot['INTRON_r'] != '')
            )
            bnd_annot['InframeSplice'] = np.where(is_inframe_splice, 1, 0)
        else:
            bnd_annot['InframeSplice'] = 0

        # sort annotated BND pairs: 
        # 1. KnownGene (descending)
        # 2. InframeSplice (descending)
        # 3. KnownTrx (descending)
        # 4. PICK (descending), SYMBOL L/R, InFrame, Distance L/R
        bnd_annot = bnd_annot.sort_values(by=["KnownGene_l", "KnownGene_r", "InframeSplice", "KnownTrx_l", "KnownTrx_r", "PICK_l", "PICK_r"], 
                              ascending=[False, False, False, False, False, False, False]).drop_duplicates(subset=["SYMBOL_l","SYMBOL_r"], keep="first")

        bnd_annot = pd.merge(bnd_annot,recurrentSvs,how='left', left_on=['SYMBOL_l','SYMBOL_r'], right_on=['KNOWNSVGENE1','KNOWNSVGENE2'])
        bnd_annot['RecurrentSV'] = bnd_annot.apply(lambda r: 1 if r['KNOWNSVGENE1'] is not None and (r['KNOWNSVTYPE'] is None or r['KNOWNSVTYPE'] == vartype) else 0, axis=1)

        # exclude pairs where neither annotation is protein coding
        bnd_annot = bnd_annot[(bnd_annot['KnownGene_l'] == 1) |
                                (bnd_annot['KnownGene_r'] == 1) |
                                ((bnd_annot['SYMBOL_l'] != bnd_annot['SYMBOL_r']) & 
                                ((bnd_annot['BIOTYPE_l']=='protein_coding') | (bnd_annot['BIOTYPE_r']=='protein_coding')))]

        if not bnd_annot.empty:
            bnd_annot['EXON_l'] = bnd_annot.apply(lambda r: _tail_or_none(r['EXON_l']) if r['SYMBOL_l']!='INTERGENIC' else None,axis=1)
            bnd_annot['INTRON_l'] = bnd_annot.apply(lambda r: _tail_or_none(r['INTRON_l']) if r['SYMBOL_l']!='INTERGENIC' else None,axis=1)
            bnd_annot['EXON_r'] = bnd_annot.apply(lambda r: _tail_or_none(r['EXON_r']) if r['SYMBOL_r']!='INTERGENIC' else None,axis=1)
            bnd_annot['INTRON_r'] = bnd_annot.apply(lambda r: _tail_or_none(r['INTRON_r']) if r['SYMBOL_r']!='INTERGENIC' else None,axis=1)
                
            bnd_annot[['genestring', 'genedetail']] = bnd_annot.apply(
                lambda row: get_gene_syntax(row, bnd_orientation, chr_l, chr_r),
                axis=1,
                result_type='expand'
            )
            
            bnd_annot = bnd_annot.sort_values(by=['RecurrentSV', 'InframeSplice'],ascending=[False,False])

            top_annot = bnd_annot.iloc[0]
            genestring = top_annot['genestring']
            genedetail = top_annot['genedetail']
            gene_l = top_annot.get('SYMBOL_l')
            gene_r = top_annot.get('SYMBOL_r')
            region_l = top_annot.get('GeneEffect_l') or "N/A"
            region_r = top_annot.get('GeneEffect_r') or "N/A"

            total_genes = len(set([g for g in [gene_l, gene_r] if g != "INTERGENIC"]))
            total_genes = f"{total_genes} gene" + ("s" if total_genes != 1 else "")
            
            if svtype == "inv" and gene_l == gene_r and region_l == region_r and ("intron" in region_l or "intragenic" in region_l):
                continue
            if svtype == "inv" and not knownTrx["Gene"].isin([gene_l, gene_r]).any():
                continue

            isRecurrentSv = (bnd_annot.iloc[0]['RecurrentSV'] and 
                            (bnd_annot.iloc[0]['InframeSplice'] or bnd_annot.iloc[0]['KNOWNSVREQUIRESTRAND'] == 0))

            category = ""
            if isRecurrentSv:
                filter_list = [f for f in filter.split(";") if f not in ["MaxDepth", "MinSomaticScore"]]
                filter = "PASS" if not filter_list else ';'.join(filter_list)

            if filter == "PASS" and isRecurrentSv:
                category = "RECURRENTSV"
            elif filter != "PASS" and isRecurrentSv:
                category = "OTHERSV"
            elif filter == "PASS" and (knownTrx["Gene"].isin([gene_l, gene_r]).any() or (gene_l != "INTERGENIC" and "ENS" not in gene_l and gene_l != "INTERGENIC" and "ENS" not in gene_r)):
                category = "OTHERSV"
            else:
                category = None

            if len(bnd_annot) > 1:
                info_dict['OTHERANNOTATIONS'] = '|'.join(bnd_annot['genedetail'].tolist()[1:])

            sv_bnd_list.append(dict(zip(svs_df_columns, [
                    category, vartype, chr_l, int(pos_l), chr_r, int(pos_r), svlen,
                    csyntax, psyntax, bandstring, genestring, genedetail, total_genes,
                    filter, str(variant.id) + ";" + str(mate.id), abundance, ';'.join([f"{k}={v}" for k, v in info_dict.items()])
                ])))
        
        alreadydone.add(variant.id)
        alreadydone.add(mate.id)

    # Concatenate only if there are records to add
    if sv_bnd_list:
        svs = pd.concat([svs.dropna(axis=1, how='all'), pd.DataFrame(sv_bnd_list)], ignore_index=True)

    # Final processing and output
    svs = svs.sort_values(by=["chrom1", "pos1", "chrom2", "pos2"], key=natsort.natsort_keygen()).reset_index(drop=True)
    svs['length'] = pd.to_numeric(svs['length'], errors='coerce').astype('Int64')
    
    # Use fillna() to fill null length values to "NA" for output
    if args.outfile:
        svs.to_csv(args.outfile, sep="\t", index=False, header=True, na_rep="NA")
        print(f"SV report written to {args.outfile}", file=sys.stderr)
    else:
        svs.to_csv(sys.stdout, sep="\t", index=False, header=True, na_rep="NA")

if __name__ == "__main__":
    main()
