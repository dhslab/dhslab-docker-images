#!/usr/bin/env python

# compile_report_data.py
# This script compiles data from various sources into a single JSON file for the Quarto report.

import argparse
import json
import pandas as pd
import os
from io import StringIO

def parse_offtarget_file(file_path):
    """Parse the off-target file and return a list of dictionaries."""
    records = []
    with open(file_path) as data:
        header = data.readline().strip().split('\t')
        for line in data:
            stripped_line = line.strip().split("\t")
            if len(stripped_line) < len(header):
                continue
            
            row_dict = dict(zip(header, stripped_line))
            records.append(row_dict)
    return records

def parse_coverage_file(file_path):
    """Parse a coverage file and return the coverage value."""
    if not file_path:
        return None
    try:
        with open(file_path, 'r') as f:
            line = f.readline().strip()
            return float(line.split(',')[-1])
    except (IOError, ValueError, IndexError):
        return None

def parse_vcf_file(file_path):
    """Parse a VCF file and return a list of dictionaries."""
    records = []
    if not file_path or not os.path.exists(file_path):
        return records

    with open(file_path, 'r') as f:
        header = []
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                header = line.strip().lstrip('#').split('\t')
                continue

            if not header:
                continue

            stripped_line = line.strip().split('\t')
            if len(stripped_line) < len(header):
                continue

            row_dict = dict(zip(header, stripped_line))
            records.append(row_dict)
    return records



def main():
    """
    Main function to parse arguments and compile data.
    """
    parser = argparse.ArgumentParser(description="Compile data for Quarto report.")
    parser.add_argument("--sample_id", required=True, help="Sample ID.")
    parser.add_argument("--control_id", required=False, default="N/A", help="Control/normal sample identifier.")
    parser.add_argument("--grnas", required=False, default="", help="Comma-separated list of gRNAs.")
    parser.add_argument("--tumor_coverage", required=False, help="Path to tumor coverage metrics file.")
    parser.add_argument("--normal_coverage", required=False, help="Path to normal coverage metrics file.")
    parser.add_argument("--cna_plot", required=False, default=None, help="Path to CNA plot PNG.")
    parser.add_argument("--baf_plot", required=False, default=None, help="Path to BAF plot PNG.")
    parser.add_argument("--circos_plot", required=False, default=None, help="Path to Circos plot PNG.")
    parser.add_argument("--transgene_insertions", required=False, help="Path to VEP-annotated on-target SV and transgene integration TSV.")
    parser.add_argument("--transgene_name", required=False, help="Transgene description string.")
    parser.add_argument("--somatic_variants", required=True, help="Path to VEP-annotated small variant TSV for targeted gene mutations.")
    parser.add_argument("--offtarget_indels", required=True, help="Path to off-target indel analysis file.")
    parser.add_argument("--offtarget_svs", required=False, help="Path to BND VCF file from indels.")
    parser.add_argument("-o", "--output", required=True, help="Output JSON file path.")
    
    args = parser.parse_args()

    # Parse the off-target indel file
    off_target_data = parse_offtarget_file(args.offtarget_indels)

    # Parse the BND VCF file
    bnd_vcf_data = parse_vcf_file(args.offtarget_svs)

    # Parse the on-target SV and transgene data
    transgene_insertions = []
    if args.transgene_insertions and os.path.exists(args.transgene_insertions):
        lines = []
        with open(args.transgene_insertions, 'r') as f:
            for line in f:
                # Skip VEP header comments
                if line.startswith('##'):
                    continue
                lines.append(line)
        
        if lines:
            # Re-join the remaining lines (header + data) into a single string
            data_str = "".join(lines)
            # Use pandas to read from the string, which handles the # in the header
            df = pd.read_csv(StringIO(data_str), sep='\\t', engine='python')
            # Clean the leading '#' from the first column name
            df.rename(columns={df.columns[0]: df.columns[0].lstrip('#')}, inplace=True)
            transgene_insertions = df.to_dict('records')

    # Parse targeted gene mutations (small variants) TSV if present
    targeted_gene_mutations = []
    try:
        sv_df = pd.read_csv(
            args.somatic_variants,
            sep='	',
            comment='#',
            keep_default_na=False,
            na_filter=False,
            dtype=str,
        )
        # Summarize results for canonical list; fall back gracefully
        genes_of_interest = ["TP53", "DNMT3A", "RUNX1", "TET2"]
        for gene in genes_of_interest:
            has_variant = False
            # Try common columns
            gene_col = None
            for cand in ["SYMBOL", "Gene", "gene", "Symbol", "symbol"]:
                if cand in sv_df.columns:
                    gene_col = cand
                    break
            consequence_col = None
            for cand in ["Consequence", "CSQ", "consequence"]:
                if cand in sv_df.columns:
                    consequence_col = cand
                    break
            if gene_col is not None:
                rows = sv_df[sv_df[gene_col] == gene]
                if not rows.empty:
                    has_variant = True
                    if consequence_col is not None:
                        effects = sorted(set(
                            ",".join(rows[consequence_col].astype(str)).split(",")
                        ))
                        targeted_gene_mutations.append({"gene": gene, "result": "; ".join([e for e in effects if e])})
                    else:
                        targeted_gene_mutations.append({"gene": gene, "result": "variant detected"})
            if not has_variant:
                targeted_gene_mutations.append({"gene": gene, "result": "no mutations identified"})
    except Exception:
        # If TSV missing or unreadable, emit default rows
        for gene in ["TP53", "DNMT3A", "RUNX1", "TET2"]:
            targeted_gene_mutations.append({"gene": gene, "result": "no data"})

    # Create a dictionary to hold all the report data.
    report_data = {
        "sample_id": args.sample_id,
        "transgene_description": args.transgene_name,
        "plots": {
            "cna": args.cna_plot,
            "baf": args.baf_plot,
            "circos": args.circos_plot
        },
        "tables": {
            "on_target_sv_transgene": transgene_insertions,
            "off_target_indels": off_target_data,
            "targeted_gene_mutations": targeted_gene_mutations,
            "bnd_vcf": bnd_vcf_data
        },
        "metadata": {
            "drug_product": args.sample_id,
#            "hotspot_file": args.hotspot_file,
            "control_sample": args.control_id,
            "assay": "WGS",
            "grnas": [s.strip() for s in args.grnas.split(",") if s.strip()],
            "mean_coverage": {
                "tumor": parse_coverage_file(args.tumor_coverage),
                "normal": parse_coverage_file(args.normal_coverage)
            }
        }
    }

    # Minimal schema validation
    def require(path, container):
        cur = container
        for key in path:
            if isinstance(cur, dict) and key in cur:
                cur = cur[key]
            else:
                raise ValueError(f"Missing required key in report_data: {'/'.join(path)}")
        return cur

    # Validate critical keys
    require(["sample_id"], report_data)
    if args.cna_plot is not None:
        require(["plots","cna"], report_data)
    if args.baf_plot is not None:
        require(["plots","baf"], report_data)

    require(["tables","on_target_sv_transgene"], report_data)
    require(["tables","off_target_indels"], report_data)
    require(["metadata","drug_product"], report_data)
    require(["metadata","control_sample"], report_data)

    # Write the compiled data to the output JSON file.
    with open(args.output, 'w') as f:
        json.dump(report_data, f, indent=4, allow_nan=False)

if __name__ == "__main__":
    main() 