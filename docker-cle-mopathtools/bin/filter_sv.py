#!/usr/bin/env python3

import argparse
import csv
import os
import sys
import gzip
import pandas as pd
import pysam

__version__ = "2.0.0"

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

# read transcript list
def get_bed_transcripts(file_path):
    df = pd.read_csv(file_path, sep="\t", header=None, names=["chrom", "start", "end", "gene", "info"])
    # extract transcript IDs
    df['transcript'] = (df["info"]
    .str.split("|", expand=True)
    .loc[:, 3].tolist())
    return df

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

# This parses the CSQ field from the VEP annotation and returns a list for a desired field
def parseVepCsq(csq, vep, field=None, allele=None):
    if not isinstance(csq, tuple) and not isinstance(vep, list):
        raise TypeError("csq should be a tuple or list")

    out = []
    for val in csq:        

        # separate the CSQ field into its components
        fields = val.split("|")

        # remove '.' from the reference allele for BNDs, e.g., 'A.' -> 'A'
        fields[0] = fields[0].strip(".")

        if allele is not None:
            if fields[0] != allele:
                continue

        if len(fields) < len(vep):
            raise ValueError("VEP header has more fields than the CSQ record")

        csq_dict = dict(zip(vep, fields))

        out.append(csq_dict)

    if field:
        out = list(set([d[field] for d in out if field in d and d[field] != ""]))

    return out

# MAIN FUNCTION
def main():
    parser = argparse.ArgumentParser(description="Filter VCF records based on criteria")
    parser.add_argument("--vcf", "-vf", type=checkfile, help="Input VCF file")
    parser.add_argument(
        "--bed_file",type=checkfile,required=True,help="Mopath-style transcript bed file",
    )
    parser.add_argument(
        "--min_length", "-l", type=int, default=1000, help="Minimum SV length"
    )
    parser.add_argument(
        "--max_length", "-L", type=int, default=20000000, help="Maximum SV length"
    )
    parser.add_argument(
        "--min_paired_reads",
        "-p",
        type=str,
        default='0',
        help="Filter code for minimum PR alt-supporting reads in the format '[hard filter],soft_filter'. Single value is assumed to be a soft filter and will set MinSvReads. Default 0.",
    )
    parser.add_argument(
        "--min_split_reads",
        "-s",
        type=str,
        default='0',
        help="Filter code for minimum SR alt-supporting reads in the format '[hard filter],soft_filter'. Single value is assumed to be a soft filter and will set MinSvReads. Default 0",
    )
    parser.add_argument(
        "--min_abundance",
        "-a",
        type=float,
        default=5.0,
        help="Minimum SV abundance in percent",
    )
    parser.add_argument("--outfile", "-o", type=fileexists, help="Outfile")
    parser.add_argument(
        "--version", "-v", action="version", version="%(prog)s: " + __version__
    )

    # List of VEP consequences that affect the protein coding sequence for DEL/DUP variants
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
        "coding_sequence_variant"
    }

    args = parser.parse_args()

    # Parse filter codes [hard,soft]. If one value is passed, assume its a soft filter.
    pr_vals = [int(x) for x in args.min_paired_reads.split(',')]
    min_pr_reads_hard, min_pr_reads_soft = (0, pr_vals[0]) if len(pr_vals) == 1 else (pr_vals[0], pr_vals[1])

    sr_vals = [int(x) for x in args.min_split_reads.split(',')]
    min_sr_reads_hard, min_sr_reads_soft = (0, sr_vals[0]) if len(sr_vals) == 1 else (sr_vals[0], sr_vals[1])

    outfile = args.outfile
    if outfile is None:
        outfile = sys.stdout

    # get transcripts
    gene_info_df = get_bed_transcripts(args.bed_file)
    knownGenes = set(gene_info_df['gene'].tolist())

    # tracker for variants to print
    passingRecords = set()

    # tracker for SV insertion hits
    SvInsertionHits = {}

    # Open VCF file for the first time and read all records
    vcf_in = pysam.VariantFile(args.vcf, "r")
    vep = vepHeader(vcf_in.header)

    if "SvInsertion" not in vcf_in.header.filters.keys():
        vcf_in.header.filters.add(
            "SvInsertion", None, None, "BND event overlaps a known insertion event from long-read sequencing"
        )

    if "Imprecise" not in vcf_in.header.filters.keys():
        vcf_in.header.filters.add(
            "Imprecise", None, None, "No contig found so breakend are imprecise"
        )

    if "MinSvReads" not in vcf_in.header.filters.keys():
        vcf_in.header.filters.add(
            "MinSvReads",
            None,
            None,
            f"Fails minimum SR or PR alt-supporting reads ({min_sr_reads_soft} and {min_pr_reads_soft}, respectively)",
        )

    if "MinSvAbundance" not in vcf_in.header.filters.keys():
        vcf_in.header.filters.add(
            "MinSvAbundance",
            None,
            None,
            f"Fails minimum SV abundance (${args.min_abundance})",
        )

    header = vcf_in.header.copy()
    vcf_out = pysam.VariantFile(outfile, "w", header=header)

    for record in vcf_in:
        svtype = record.info.get("SVTYPE")

        svlen = None
        if record.info.get("SVLEN") is not None:
            svlen = abs(record.info.get("SVLEN")[0])
        elif record.stop is not None:
            svlen = record.stop - record.pos
        else:
            svlen = abs(len(record.alts[0])-len(record.ref))

        # this gets genes that overlap SVs from both the standard VEP annotation to get upstream/downstream events and custom BED overlap annotations
        # with gene clusters, like IGH, IGL, TRA, etc.
        genes = []
        consequences = []

        if record.info.get("CSQ") is not None or record.info.get("KnownSvGenes") is not None:
            genes = list(
                set(
                    (parseVepCsq(csq, list(vep.keys()), "SYMBOL", None) if (csq := record.info.get("CSQ")) else []) +
                    ([gene.split("|")[0] for gene in list(ksg)] if (ksg := record.info.get("KnownSvGenes")) else [])
                )
            )
            consequences = "&".join(
                set(parseVepCsq(csq, list(vep.keys()), "Consequence", None) if (csq := record.info.get("CSQ")) else [])
            ).split("&")

        # Record whether a BND event is a known insertion event
        if svtype == "BND" and record.info.get("SvInsertionHits"):
            SvInsertionHits[record.id] = record.info.get("SvInsertionHits")

        PR = (0, 0)
        SR = (0, 0)
        if "PR" in record.format.keys():
            if "DUX" in record.id:
                record.samples[0]["PR"] = (0, record.samples[0]["PR"][1])
            PR = record.samples[0]["PR"]

        if "SR" in record.format.keys():
            if "DUX" in record.id:
                record.samples[0]["SR"] = (0, record.samples[0]["SR"][1])
            SR = record.samples[0]["SR"]

        # Hard filter for SVs with no read support or that fail hard filter rules.
        if (PR[1] == 0 and SR[1] == 0) or PR[1] <= min_pr_reads_hard or SR[1] <= min_sr_reads_hard:
            continue

        # also skip DEL/DUP calls with the SystematicNoise filter set
        # that have either 0 SR or PR reads or that cross the centromere. 
        if (svtype in ["DEL", "DUP"] 
            and (
                (PR[1] == 0 or SR[1] == 0)
                or
                (any(item.find('p') > -1 for item in record.info.get("Cytobands",())) 
                 and any(item.find('q') > -1 for item in record.info.get("Cytobands",())))
            )
        ):
            continue

        geneHits = set(genes) & knownGenes
    
        recordVariant = False
        # get BND records that overlap a known gene or that PASS
        if (svtype == "BND" 
            and ((len(geneHits) > 0 or record.info.get("MATEID")[0] in passingRecords)
                 or ('PASS' in record.filter.keys()
                        and len(set(consequences) & nonSynon) > 0
                        and record.info.get("CONTIG") is not None)
                or 'DUX4' in record.id)
            ):
                recordVariant = True

        # INS/DEL/DUP must pass length filters and affect a known gene
        elif (svtype in ["INS", "DEL", "DUP"] 
              and svlen is not None
              and svlen >= args.min_length 
              and svlen < args.max_length                  
              and len(geneHits) > 0
              and len(set(consequences) & nonSynon) > 0
            ):
                recordVariant = True        

        if recordVariant:
            passingRecords.add(record.id)

            record.filter.clear()

            # apply filters
            if "IMPRECISE" in record.info.keys():
                record.filter.add("Imprecise")            

            if PR[1] < int(min_pr_reads_soft) or SR[1] < int(min_sr_reads_soft):
                record.filter.add("MinSvReads")

            if record.info.get("SvInsertionHits") is not None:
                record.filter.add("SvInsertion")

            if (PR[1] + SR[1]) / (PR[0] + SR[0] + PR[1] + SR[1]) * 100.0 < float(
                args.min_abundance
            ):
                record.filter.add("MinSvAbundance")

            if len(record.filter.keys()) == 0:
                record.filter.add("PASS")

            # Fix BND records
            if record.info.get("SVTYPE") == "BND":
                chr1 = record.contig
                if record.alts[0].startswith("["):
                    chr2, pos2 = record.alts[0].split("[")[1].split(":")
                    if chr1 == chr2:
                        record.info["SVTYPE"] = "INV"

                if record.alts[0].endswith("]"):
                    chr2, pos2 = record.alts[0].split("]")[1].split(":")
                    if chr1 == chr2:
                        record.info["SVTYPE"] = "INV"

                # if DUX4 call, then adjust the VCF records accordingly.
                if 'DUX4' in record.id:
                    # change filter to PASS
                    record.filter.clear()
                    record.filter.add("PASS")
                    if record.info.get("KnownSvGenes") is None:
                        record.info['KnownSvGenes'] = 'DUX4|DUX4|DUX4|transcript|||'

            vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()


if __name__ == "__main__":
    main()
