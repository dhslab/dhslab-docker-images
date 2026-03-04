#!/usr/bin/env python

import argparse
import os
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
        missing_symbol = df["SYMBOL"].isna() | (df["SYMBOL"] == "")
        df.loc[missing_symbol, "SYMBOL"] = df.loc[missing_symbol, "Gene"]
    
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

def vep_gene_effect(row):
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

def load_recurrent_svs(path: str) -> pd.DataFrame:
    recurrent_svs = pd.read_csv(path, sep=",", header=None)
    swapped = recurrent_svs.copy()
    swapped.iloc[:, [0, 1]] = swapped.iloc[:, [1, 0]].values
    recurrent_svs = pd.concat([recurrent_svs, swapped], axis=0, ignore_index=True)
    recurrent_svs.columns = ["KNOWNSVGENE1", "KNOWNSVGENE2", "KNOWNSVREQUIRESTRAND", "KNOWNSVTYPE"]
    return recurrent_svs.where(pd.notna(recurrent_svs), None)

def _extract_trx_info(cov_df: pd.DataFrame) -> pd.DataFrame:
    ids = (
        cov_df["Info"]
        .str.split("|", expand=True)
        .replace("\s\+\s\S+", "", regex=True)
        .loc[:, 2:]
    )
    ids.columns = ["Transcript", "Region", "cdsStart", "cdsEnd", "strand"]
    return pd.concat([cov_df, ids], axis=1)[["Gene", "Transcript", "cdsStart", "cdsEnd", "strand"]].drop_duplicates()

def load_cov_gene_transcript_sets(path: str) -> tuple[pd.DataFrame, set[str], set[str], set[str]]:
    with open(path, "r") as f:
        header_line = f.readline().replace("#", "").strip()
    cov_df = pd.read_csv(
        path,
        header=None,
        skiprows=1,
        names=["Chromosome", "Start", "End", "Gene", "Info"] + header_line.split("\t")[5:],
        sep="\t",
    )
    gene_cov_df = cov_df[cov_df["Info"].str.contains("exon")]
    sv_cov_df = cov_df[cov_df["Info"].str.contains(r"transcript|gene", regex=True)]

    gene_trx = _extract_trx_info(gene_cov_df)
    sv_trx = _extract_trx_info(sv_cov_df)
    known_trx = pd.concat(
        [sv_trx[["Gene", "Transcript"]], gene_trx[["Gene", "Transcript"]]],
        axis=0,
        ignore_index=True,
    ).drop_duplicates()

    known_gene_set = set(known_trx["Gene"].dropna().tolist())
    known_trx_set = set(known_trx["Transcript"].dropna().tolist())
    reportable_gene_set = set(gene_trx["Gene"].dropna().drop_duplicates().tolist())
    return known_trx, known_gene_set, known_trx_set, reportable_gene_set

def _base_filter(variant) -> str:
    filter_keys = list(variant.filter.keys())
    return "PASS" if not filter_keys else ";".join(filter_keys)

def _clean_recurrent_filter(filter_value: str) -> str:
    filter_list = [f for f in filter_value.split(";") if f not in ["MaxDepth", "MinSomaticScore"]]
    return "PASS" if not filter_list else ";".join(filter_list)

def _extract_support_reads(variant, zero_ref_for_dux: bool = False) -> tuple[tuple[int, int], tuple[int, int], float]:
    pr = (0, 0)
    sr = (0, 0)
    if "SR" in variant.samples[0] and variant.samples[0]["SR"][0] is not None:
        sr = variant.samples[0]["SR"]
    if "PR" in variant.samples[0] and variant.samples[0]["PR"][0] is not None:
        pr = variant.samples[0]["PR"]

    if zero_ref_for_dux and "DUX" in variant.id:
        sr = (0, sr[1])
        pr = (0, pr[1])

    total_reads = pr[0] + pr[1] + sr[0] + sr[1]
    supporting_reads = sr[1] + pr[1]
    abundance = round(supporting_reads / total_reads * 100, 1) if total_reads > 0 else 0.0
    return pr, sr, abundance

def _build_info_dict(pr: tuple[int, int], sr: tuple[int, int], contig) -> dict[str, str]:
    return {
        "PR_READS": f"{pr[1]}/{pr[0] + pr[1]}",
        "SR_READS": f"{sr[1]}/{sr[0] + sr[1]}",
        "CONTIG": str(contig),
    }

def _serialize_info_dict(info_dict: dict[str, str]) -> str:
    return ";".join([f"{k}={v}" for k, v in info_dict.items()])

def _get_bands(variant) -> list[str]:
    if "Cytobands" in variant.info and variant.info.get("Cytobands") is not None:
        return [b.replace("acen_", "") for b in variant.info.get("Cytobands")]
    return []

def _get_bandstring(bands: list[str], default: str = "None") -> str:
    if not bands:
        return default
    if len(bands) == 1:
        return bands[0]
    return bands[0] + bands[-1]

def _get_first_band(variant, default: str = "N/A") -> str:
    bands = _get_bands(variant)
    return bands[0] if bands else default

def _merge_known_sv_genes(vep_df: pd.DataFrame, known_sv_df: pd.DataFrame) -> pd.DataFrame:
    if known_sv_df.empty:
        return vep_df
    merged = pd.concat([vep_df, known_sv_df[~known_sv_df["SYMBOL"].isin(vep_df["SYMBOL"])]], ignore_index=True)
    return merged.where(pd.notna(merged), None)

def _annotate_known_flags(
    vep_df: pd.DataFrame,
    known_gene_set: set[str],
    known_trx_set: set[str],
    reportable_gene_set: set[str] | None = None,
) -> pd.DataFrame:
    if vep_df.empty:
        return vep_df
    vep_df["KnownGene"] = vep_df["SYMBOL"].isin(known_gene_set).astype(int)
    if reportable_gene_set is not None:
        vep_df.loc[vep_df["SYMBOL"].isin(reportable_gene_set), "KnownGene"] = 1
    vep_df["KnownTrx"] = vep_df["Feature"].isin(known_trx_set).astype(int)
    return vep_df

def _build_breakend_vep(
    variant,
    csq_header_desc: str,
    known_sv_genes_header_desc: str,
    known_gene_set: set[str],
    known_trx_set: set[str],
) -> pd.DataFrame:
    vep_df = vep_csq_to_dataframe(variant.info.get("CSQ"), csq_header_desc)
    vep_df = vep_df[(vep_df["Allele"] == variant.ref) | (vep_df["Allele"] == "BND")].reset_index(drop=True)
    if vep_df.empty:
        return pd.DataFrame([{}], columns=vep_df.columns)

    vep_df = _annotate_known_flags(vep_df, known_gene_set, known_trx_set)
    if "IntronFrame" not in vep_df.columns:
        vep_df["IntronFrame"] = None

    known_sv_df = known_sv_genes_tag_to_dataframe(variant, known_sv_genes_header_desc)
    known_sv_df["IntronFrame"] = 0
    vep_df = _merge_known_sv_genes(vep_df, known_sv_df)
    vep_df["GeneEffect"] = vep_df.apply(lambda r: vep_gene_effect(r), axis=1)
    return vep_df

def _build_deldup_syntax(vartype: str, chrom: str, pos1: int, pos2, bands: list[str]) -> tuple[str, str]:
    csyntax = ""
    psyntax = ""
    if vartype == "DEL":
        csyntax = chrom + ":g." + str(pos1) + "_" + str(pos2) + "del"
        if (bands and bands[0].find("p") > -1 and bands[-1].find("q") > -1):
            psyntax = "seq[GRCh38] -" + chrom.replace("chr", "")
        elif (bands and any("q11" in b for b in bands) and any("qter" in b for b in bands) and chrom in ACROCENTRICS):
            psyntax = "seq[GRCh38] -" + chrom.replace("chr", "")
        elif (bands and bands[0].find("p") > -1):
            psyntax = ("seq[GRCh38] del(" + chrom.replace("chr", "") + ")(" + bands[0] + bands[-1] + ")")
        elif bands:
            psyntax = ("seq[GRCh38] del(" + chrom.replace("chr", "") + ")(" + bands[0] + bands[-1] + ")")
    elif vartype in ["DUP", "INS"]:
        csyntax = chrom + ":g." + str(pos1) + "_" + str(pos2) + vartype.lower()
        if (bands and bands[0].find("p") > -1 and bands[-1].find("q") > -1):
            psyntax = "seq[GRCh38] +" + chrom.replace("chr", "")
        elif (bands and any("q11" in b for b in bands) and any("qter" in b for b in bands) and chrom in ACROCENTRICS):
            psyntax = "seq[GRCh38] +" + chrom.replace("chr", "")
        elif (bands and bands[0].find("p") > -1):
            psyntax = ("seq[GRCh38] " + vartype.lower() + "(" + chrom.replace("chr", "") + ")(" + bands[0] + bands[-1] + ")")
        elif bands:
            psyntax = ("seq[GRCh38] " + vartype.lower() + "(" + chrom.replace("chr", "") + ")(" + bands[0] + bands[-1] + ")")
    return csyntax, psyntax

def _normalize_exon_intron_columns(df: pd.DataFrame) -> pd.DataFrame:
    def _tail_or_none(value):
        if isinstance(value, str) and value:
            return value.split("-")[-1]
        return None

    df["EXON_l"] = df["EXON_l"].apply(_tail_or_none)
    df["INTRON_l"] = df["INTRON_l"].apply(_tail_or_none)
    df["EXON_r"] = df["EXON_r"].apply(_tail_or_none)
    df["INTRON_r"] = df["INTRON_r"].apply(_tail_or_none)
    return df

def main():
    parser = argparse.ArgumentParser(description="Collect Structural Variants from ChromoSeq output.")
    parser.add_argument("--sv_vcf", required=True, type=checkfile, help="Annotated SV VCF file (*.sv_annotated.vcf.gz)")
    parser.add_argument("--sv_targets", required=True, type=checkfile, help="Recurrent SV target list CSV")
    parser.add_argument("--cov_report", required=True, type=checkfile, help="Coverage report TSV (*.coverage_report.tsv)")
    parser.add_argument("--outfile", type=str, help="Output tab-delimited file. If not provided, output is written to stdout.")
    parser.add_argument("--reference", help="Reference fasta (required for CRAM, optional for VCF)")
    parser.add_argument("--version", "-v", action="version", version="%(prog)s: " + __version__)
    args = parser.parse_args()

    # Load accessory files
    nonsynonymous_consequences = {
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

    recurrent_svs = load_recurrent_svs(args.sv_targets)
    known_trx, known_gene_set, known_trx_set, reportable_gene_set = load_cov_gene_transcript_sets(args.cov_report)

    # Initialize dataframe for all SVs
    svs_df_columns = [
        "category", "type", "chrom1", "pos1", "chrom2", "pos2", "length",
        "csyntax", "psyntax", "bands", "known_genes", "known_gene_detail",
        "total_genes", "filters", "id", "abundance", "info",
    ]
    svs = pd.DataFrame(columns=svs_df_columns)

    # Process SVs
    print("Gathering SVs...", file=sys.stderr)
    sv_vcf_handle = pysam.VariantFile(args.sv_vcf)
    csq_header_desc = sv_vcf_handle.header.info['CSQ'].description
    known_sv_genes_header_desc = ""
    if 'KnownSvGenes' in sv_vcf_handle.header.info:
        known_sv_genes_header_desc = sv_vcf_handle.header.info['KnownSvGenes'].description

    pending_breakend_variants = {}
    sv_deldup_list = []

    for variant in sv_vcf_handle:
        
        if variant.info["SVTYPE"] in ["BND", "INV"]:
            pending_breakend_variants[variant.id] = variant
        
        else:
            # Hard filter DEL/DUPs that do not involve genes
            if "KnownSvGenes" not in variant.info and "CSQ" not in variant.info:
                continue

            filter = _base_filter(variant)

            #caseid = variant.info.get("CASE")

            vartype = variant.info.get("SVTYPE")
            chr1 = str(variant.chrom)
            pos1 = variant.pos
            svlen = None
            if variant.info.get("SVLEN") is not None:
                svlen = int(abs(variant.info.get("SVLEN")[0]))

            pos2 = variant.stop if variant.stop is not None else "N/A"
            if pos2 == "N/A" and svlen != None:
                pos2 = pos1 + svlen
            
            pr, sr, abundance = _extract_support_reads(variant)
            info_dict = _build_info_dict(pr, sr, variant.info.get("CONTIG"))

            csyntax = ""
            psyntax = ""
            bandstring = "None"
            genestring = "None"
            genedetail = "None"
            total_genes = "None"

            bands = _get_bands(variant)
            bandstring = _get_bandstring(bands)

            vep_csq = vep_csq_to_dataframe(variant.info.get("CSQ"), csq_header_desc)
            vep_csq = _annotate_known_flags(vep_csq, known_gene_set, known_trx_set, reportable_gene_set)
            vep_csq = vep_csq.sort_values(
                by=["SYMBOL", "KnownTrx", "PICK"], ascending=[True, False, False]
            ).drop_duplicates(subset="SYMBOL", keep="first")
        
            known_sv_gene_df = known_sv_genes_tag_to_dataframe(variant, known_sv_genes_header_desc)
            vep_csq = _merge_known_sv_genes(vep_csq, known_sv_gene_df)

            if not vep_csq.empty:
                total_genes = f"{len(vep_csq['SYMBOL'].unique())} genes"
                vep_csq = vep_csq.sort_values(by=["START"], ascending=[True])
                vep_csq["GeneEffect"] = vep_csq.apply(lambda r: vep_gene_effect(r), axis=1)
                vep_csq["GeneImpact"] = vep_csq["Consequence"].apply(
                    lambda r: int(len(set(r.split("&")) & set(nonsynonymous_consequences)) > 0 if r else None)
                )

            csyntax, psyntax = _build_deldup_syntax(vartype, chr1, pos1, pos2, bands)

            known_gene_df = vep_csq[(vep_csq["KnownTrx"] == 1) & (vep_csq["DISTANCE"] == 0)].sort_values(by=["START"])[["SYMBOL", "GeneImpact", "GeneEffect"]]
            
            if len(known_gene_df) > 15:
                genestring = f"{len(known_gene_df)} genes"
                genedetail = genestring
            elif len(known_gene_df) > 0:
                genestring = ",".join(known_gene_df["SYMBOL"].unique().tolist())
                genedetail = (
                    vartype.lower() + "[" + ",".join(
                        known_gene_df.apply(
                            lambda x: f"{x['SYMBOL']}({x['GeneEffect']})" if x['GeneEffect']!="intragenic" else f"{x['SYMBOL']}", axis=1
                        ).tolist()
                    ) + "]"
                )

            is_recurrent_sv = ((recurrent_svs["KNOWNSVGENE1"].isin(known_gene_df['SYMBOL']))
                            & (recurrent_svs["KNOWNSVGENE2"].isin(known_gene_df['SYMBOL']))
                            & (recurrent_svs["KNOWNSVREQUIRESTRAND"] == 1)
                            & (recurrent_svs["KNOWNSVTYPE"].isin([vartype,"*"]))).any()
            
            if is_recurrent_sv:
                filter = _clean_recurrent_filter(filter)

            category = ""
            if filter == "PASS" and is_recurrent_sv:
                category = "CNV"
            elif filter != "PASS" and is_recurrent_sv:
                category = "OTHERSV"
            elif filter == "PASS" and genestring != "None":
                category = "OTHERSV"
            else:
                category = None

            #info_dict['CASEID'] = caseid

            # Add DEL/DUP info to list
            sv_deldup_list.append(dict(zip(svs_df_columns, [
                    category, vartype, chr1, int(pos1), chr1, int(pos2), svlen,
                    csyntax, psyntax, bandstring, genestring, genedetail, total_genes,
                    filter, str(variant.id), abundance, _serialize_info_dict(info_dict),
                ])))
            
            # Get fusion genes formed by a deletion if the deletion juxtaposes 2 genes
            if vartype == "DEL" and len(vep_csq[(vep_csq['START']<=pos1)]) > 1 and len(vep_csq[(vep_csq['END']>=pos2)]) > 1:
                # Evaluate all left/right transcript pairings spanning the deleted interval.
                candidate_fusions = pd.merge(vep_csq[(vep_csq['START']<=pos1)], 
                                             vep_csq[(vep_csq['END']>=pos2)], how='cross', suffixes=('_l', '_r'))

                is_inframe_splice = (
                    (candidate_fusions['STRAND_l'] == candidate_fusions['STRAND_r']) &
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

                candidate_fusions = pd.merge(candidate_fusions,recurrent_svs,how='left', left_on=['SYMBOL_l','SYMBOL_r'], right_on=['KNOWNSVGENE1','KNOWNSVGENE2'])
                is_recurrent = (candidate_fusions['KNOWNSVGENE1'].notna()) & \
                               (candidate_fusions['KNOWNSVTYPE'].isna() | candidate_fusions['KNOWNSVTYPE'].isin([vartype, '*']))
                candidate_fusions['RecurrentSV'] = np.where(is_recurrent, 1, 0)
                candidate_fusions = candidate_fusions[(candidate_fusions['RecurrentSV'] == 1) |
                                                        ((candidate_fusions['SYMBOL_l'] != candidate_fusions['SYMBOL_r']) & 
                                                        ((candidate_fusions['BIOTYPE_l']=='protein_coding') | (candidate_fusions['BIOTYPE_r']=='protein_coding')))]

                if not candidate_fusions.empty:
                    candidate_fusions = _normalize_exon_intron_columns(candidate_fusions)
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

                    #info_dict['CASEID'] = caseid

                    sv_deldup_list.append(dict(zip(svs_df_columns, [
                        category, vartype, chr1, int(pos1), chr1, int(pos2), svlen,
                        csyntax, psyntax, bandstring, candidate_fusions.iloc[0]['genestring'], candidate_fusions.iloc[0]['genedetail'], '2 genes',
                        filter, str(variant.id), abundance, _serialize_info_dict(info_dict),
                    ])))

    # Concatenate only if there are records to add
    if sv_deldup_list:
        svs = pd.concat([svs.dropna(axis=1, how='all'), pd.DataFrame(sv_deldup_list)], ignore_index=True)

    # Pair breakends once by variant id and mate id.
    processed_variant_ids = set()
    sv_bnd_list = []
    for _variant_id, variant in pending_breakend_variants.items():
        
        mate_id_tuple = variant.info.get("MATEID")
        if not mate_id_tuple:
            continue

        mate_id = mate_id_tuple[0]
        
        if variant.id in processed_variant_ids or mate_id in processed_variant_ids or mate_id not in pending_breakend_variants:
            continue
        
        mate = pending_breakend_variants[mate_id]

        if ("KnownSvGenes" not in variant.info and "KnownSvGenes" not in mate.info) and \
           ("CSQ" not in variant.info or "CSQ" not in mate.info):
            continue

        # Normalize ordering so downstream logic can treat this as left/right consistently.
        if variant.alts[0].find("[") == 0 or variant.alts[0].find("]") == 0:
            variant, mate = mate, variant

        bnd_orientation = -1 if variant.alts[0].find("[") == 0 or variant.alts[0].find("]") > 0 else 1

        #caseid = variant.info.get("CASE")
        
        vartype = variant.info.get("SVTYPE")
        filter = _base_filter(variant)

        chr_l, pos_l = str(variant.chrom), variant.pos
        chr_r, pos_r = str(mate.chrom), mate.pos
        svlen = None
        genestring = "N/A"
        genedetail = "N/A"
        total_genes = "N/A"
        bandstring = "N/A"

        bands_l = _get_first_band(variant, default="N/A")
        bands_r = _get_first_band(mate, default="N/A")

        # Build ISCN/cytogenetic syntax for BND/INV outputs.
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

        pr, sr, abundance = _extract_support_reads(variant, zero_ref_for_dux=True)
        info_dict = _build_info_dict(pr, sr, variant.info.get("CONTIG"))

        left_vep_csq = _build_breakend_vep(
            variant,
            csq_header_desc,
            known_sv_genes_header_desc,
            known_gene_set,
            known_trx_set,
        )
        right_vep_csq = _build_breakend_vep(
            mate,
            csq_header_desc,
            known_sv_genes_header_desc,
            known_gene_set,
            known_trx_set,
        )

        # Cross left/right annotations to score all possible gene pairings.
        bnd_annot = pd.merge(left_vep_csq, right_vep_csq, how='cross', suffixes=('_l', '_r'))
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

        bnd_annot = pd.merge(bnd_annot,recurrent_svs,how='left', left_on=['SYMBOL_l','SYMBOL_r'], right_on=['KNOWNSVGENE1','KNOWNSVGENE2'])
        is_recurrent = (bnd_annot['KNOWNSVGENE1'].notna()) & \
                       (bnd_annot['KNOWNSVTYPE'].isna() | bnd_annot['KNOWNSVTYPE'].isin([vartype, '*']))
        bnd_annot['RecurrentSV'] = np.where(is_recurrent, 1, 0)

        # exclude pairs where neither annotation is protein coding
        bnd_annot = bnd_annot[(bnd_annot['KnownGene_l'] == 1) |
                                (bnd_annot['KnownGene_r'] == 1) |
                                ((bnd_annot['SYMBOL_l'] != bnd_annot['SYMBOL_r']) & 
                                ((bnd_annot['BIOTYPE_l']=='protein_coding') | (bnd_annot['BIOTYPE_r']=='protein_coding')))]

        if not bnd_annot.empty:
            bnd_annot = _normalize_exon_intron_columns(bnd_annot)
            bnd_annot.loc[bnd_annot["SYMBOL_l"] == "INTERGENIC", ["EXON_l", "INTRON_l"]] = None
            bnd_annot.loc[bnd_annot["SYMBOL_r"] == "INTERGENIC", ["EXON_r", "INTRON_r"]] = None
                
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
            if svtype == "inv" and not any(g in known_gene_set for g in [gene_l, gene_r]):
                continue

            is_recurrent_sv = (bnd_annot.iloc[0]['RecurrentSV'] and 
                            (bnd_annot.iloc[0]['InframeSplice'] or bnd_annot.iloc[0]['KNOWNSVREQUIRESTRAND'] == 0))

            category = ""
            if is_recurrent_sv:
                filter = _clean_recurrent_filter(filter)

            if filter == "PASS" and is_recurrent_sv:
                category = "RECURRENTSV"
            elif filter != "PASS" and is_recurrent_sv:
                category = "OTHERSV"
            elif filter == "PASS" and (any(g in known_gene_set for g in [gene_l, gene_r]) or (gene_l != "INTERGENIC" and "ENS" not in gene_l and gene_l != "INTERGENIC" and "ENS" not in gene_r)):
                category = "OTHERSV"
            else:
                category = None

            if len(bnd_annot) > 1:
                info_dict['OTHERANNOTATIONS'] = '|'.join(bnd_annot['genedetail'].tolist()[1:])

            #info_dict['CASEID'] = caseid

            sv_bnd_list.append(dict(zip(svs_df_columns, [
                    category, vartype, chr_l, int(pos_l), chr_r, int(pos_r), svlen,
                    csyntax, psyntax, bandstring, genestring, genedetail, total_genes,
                    filter, str(variant.id) + ";" + str(mate.id), abundance, _serialize_info_dict(info_dict)
                ])))
        
        processed_variant_ids.add(variant.id)
        processed_variant_ids.add(mate.id)

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
