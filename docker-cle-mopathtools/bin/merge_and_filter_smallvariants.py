#!/usr/bin/env python3

import argparse
import os
import re
import sys

import pysam

__version__ = "1.0.0"


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

def natural_sort_key(s):
    """A sorting key for natural string order."""
    return [
        int(text) if text.isdigit() else text.lower()
        for text in re.split("([0-9]+)", s)
    ]

def convert_symbolic_alt(record, faobj):
    if record.alts[0] == "<DEL>":
        record.ref = faobj.fetch(record.contig, record.pos - 1, record.stop)
        record.alts = (str(faobj.fetch(record.contig, record.pos - 1, record.pos)),)

    elif record.alts[0] == "<DUP>" or record.alts[0] == "<TANDEM:DUP>":
        record.ref = faobj.fetch(record.contig, record.pos - 1, record.pos)
        record.alts = (str(faobj.fetch(record.contig, record.pos - 1, record.stop)),)

    elif (
        record.alts[0] == "<INS>"
        and "CONTIG" in record.info
        and "NN" not in record.info.get("CONTIG")
    ):
        # insertions have a novel sequence, so we must get it from the contig.
        # there are no coordinates for this insertion on the contig, so we must discover them
        # by search for the refseq before and after it.
        record.ref = faobj.fetch(record.contig, record.pos - 1, record.pos)

        contigSeq = record.info.get("CONTIG")
        upstreamSeq = faobj.fetch(record.contig, record.pos - 20, record.pos)
        downstreamSeq = faobj.fetch(record.contig, record.pos, record.pos + 20)
        # create aligner
        aligner = PairwiseAligner()
        aligner.mode = "local"
        # Set the match, mismatch, open gap and extend gap scores
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -1

        # Perform the alignment
        alnupstream = aligner.align(upstreamSeq, contigSeq)
        alndownstream = aligner.align(downstreamSeq, contigSeq)

        # set alt allele if upstream and downstream seqs are found
        if len(alnupstream) > 0 and len(alndownstream) > 0:
            ins_start = alnupstream[0].aligned[1, 0][1] - 1
            ins_end = alndownstream[-1].aligned[1, 0][0]
            if ins_start < ins_end:
                ins_seq = contigSeq[ins_start:ins_end]
                record.alts = (ins_seq,)

    return record

def main():
    parser = argparse.ArgumentParser(
        description="Merge small variants and gene-level SVs and apply filtering criteria"
    )
    parser.add_argument(
        "--small-variant-vcf", "-sm", type=checkfile, help="Small variant VCF file"
    )
    parser.add_argument("--sv-vcf", "-sv", type=checkfile, help="SV VCF file")

    parser.add_argument("-c", "--cramfile", type=checkfile, help="CRAM alignment file")
    parser.add_argument("-r", "--reference", type=checkfile, help="Reference fasta file")
        
    parser.add_argument(
        "--bed-file",
        "-b",
        type=checkfile,
        required=True,
        help="Bed file to filter VCFs",
    )
    parser.add_argument(
        "--max-sv-length",
        "-L",
        type=int,
        default=1000,
        help="Maximum gene-level SV length",
    )
    parser.add_argument(
        "--min-reads",
        "-m",
        type=int,
        default=3,
        help="Minimum alt-supporting reads to pass LowReads filter",
    )
    parser.add_argument(
        "--vaf-hard-filter", 
        type=float, 
        default=0.0, 
        help="Minimum VAF to pass hard cutoff. Variants below this VAF will be filtered out regardless of other criteria. (%)"
    )

    parser.add_argument(
        "--vaf-soft-filter",
        type=float,
        default=2.0,
        help="Minimum VAF to pass MinVAF soft filter. (%)",
    )
    parser.add_argument(
        "--allowed-info-tags", 
        type=str, 
        default=None, 
        help="VCF INFO tags that allow variants to bypass minreads and minvaf filters."
    )
    parser.add_argument("--outfile", "-o", type=fileexists, help="Outfile")
    parser.add_argument(
        "--version", "-v", action="version", version="%(prog)s: " + __version__
    )

    args = parser.parse_args()

    allowed_info_tags = ['FGT'] + (args.allowed_info_tags.split(",") if args.allowed_info_tags is not None else [])

    outfile = args.outfile
    if outfile is None:
        outfile = sys.stdout

    # Open VCF file
    vcf = pysam.VariantFile(args.small_variant_vcf, "r")
    svVcf = pysam.VariantFile(args.sv_vcf, "r")

    header = vcf.header.copy()
    vcf_out = pysam.VariantFile(outfile, "w", header=header)
    vcf_out.header.merge(svVcf.header)

    if "LowReads" not in vcf_out.header.filters.keys():
        vcf_out.header.filters.add(
            "LowReads",
            None,
            None,
            f"Fails minimum alt-supporting reads ({args.min_reads})",
        )

    if "MinVAF" not in vcf_out.header.filters.keys():
        vcf_out.header.filters.add(
            "MinVAF", None, None, f"Fails minimum VAF filter ({args.vaf_soft_filter})"
        )

    if 'ImpreciseSV' not in vcf_out.header.filters.keys():
        vcf_out.header.filters.add("ImpreciseSV",None,None,f'Small SV with imprecise breakends')

    cramfile = None
    if args.cramfile is not None and args.reference is not None:
        cramfile = pysam.AlignmentFile(args.cramfile, "rc", reference_filename=args.reference)

    regions = pysam.TabixFile(args.bed_file)

    variants = []
    variantDict = {}

    for row in regions.fetch(parser=pysam.asTuple()):
        reg = f"{row[0]}:{int(row[1])+1}-{row[2]}" # note that bed is 0-based, region is 1-based: https://pysam.readthedocs.io/en/latest/api.html#pysam.VariantFile.fetch

        for record in vcf.fetch(region=reg):

            # skip duplicates, which can happen if they overlap more than one region
            if f"{record.contig}-{record.pos}-{record.ref}-{record.alts[0]}" in variantDict:
                continue

            variantDict[f"{record.contig}-{record.pos}-{record.ref}-{record.alts[0]}"] = True

            new_record = vcf_out.new_record()

            new_record.chrom = record.chrom
            new_record.pos = record.pos
            new_record.stop = record.stop
            new_record.id = record.id
            new_record.alleles = record.alleles

            for k in record.filter.keys():
                new_record.filter.add(k)

            # change variants filtered only w/ multiallelic to PASS
            if 'multiallelic' in new_record.filter and len(new_record.filter.keys())==1:
                new_record.filter.clear()
                new_record.filter.add("PASS")

            for k in record.info.keys():
                new_record.info[k] = record.info.get(k)

            for k in record.format.keys():
                new_record.samples[0][k] = record.samples[0][k]

            variants.append(new_record)

        # Now get sv variants in this region and add to the list, applying filters for small SVs
        for record in svVcf.fetch(region=reg):
            # skip large SV records and
            if record.info["SVTYPE"] == "BND" or (
                "SVLEN" in record.info.keys()
                and abs(record.info["SVLEN"][0]) > int(args.max_sv_length)
            ):
                continue

            # skip duplicates, which can happen if they overlap more than one region
            if f"{record.contig}-{record.pos}-{record.ref}-{record.alts[0]}" in variantDict:
                continue

            variantDict[f"{record.contig}-{record.pos}-{record.ref}-{record.alts[0]}"] = True
            
            new_record = vcf_out.new_record()

            new_record.chrom = record.chrom
            new_record.pos = record.pos
            new_record.stop = record.stop
            new_record.id = record.id
            new_record.alleles = record.alleles

            # convert symbolic SVs to sequence-resolved small variants if possible, which allows us to apply read-based filters to them
            if new_record.alts[0] in ["<DEL>", "<DUP>", "<TANDEM:DUP>", "<INS>"]:
                new_record = convert_symbolic_alt(new_record, pysam.FastaFile(args.reference))

            for k in record.filter.keys():
                new_record.filter.add(k)

            for k in record.info.keys():
                new_record.info[k] = record.info.get(k)

            for k in record.format.keys():
                new_record.samples[0][k] = record.samples[0][k]

            PR = (0, 0)
            SR = (0, 0)
            if "PR" in record.format.keys():
                PR = record.samples[0]["PR"]

            if "SR" in record.format.keys():
                SR = record.samples[0]["SR"]

            DP = PR[0] + PR[1] + SR[0] + SR[1]
            
            AD = (PR[0] + SR[0], PR[1] + SR[1])

            # skip small SV variants if there's no read support
            if DP == 0 or AD[1] == 0:
                continue

            AF = round((PR[1] + SR[1]) / DP, 4)

            new_record.samples[0]["DP"] = DP
            new_record.samples[0]["AD"] = AD
            new_record.samples[0]["AF"] = AF            

            variants.append(new_record)


    # sort variants by chrom and position
    for record in sorted(variants, key=lambda x: (natural_sort_key(x.chrom), x.pos)):
        
        # If this variant is NOT a prior variant or a force_gt variant (indicated by INFO tags in allowed_info_tags)
        # then do not include variants with < minreads, < hard-vaf cutoff, or if its not a PASS variant and below the soft vaf cutoff
        if not any(tag in record.info for tag in allowed_info_tags):
            if ((not 'AD' in record.format.keys() or 
                record.samples[0]["AD"][1] < args.min_reads or 
                record.samples[0]["AF"][0] * 100 < args.vaf_hard_filter) or 
                (not 'PASS' in record.filter and record.samples[0]["AF"][0] * 100 < args.vaf_soft_filter)):
                continue

        # Variants without AD tags are not called by Dragen.
        # Get the depth for these manually
        if 'AD' not in record.format.keys():
            if args.cramfile is None or args.reference is None:
                raise RuntimeError(f"--cramfile and --reference are required to compute depth for variants without AD (e.g. {record.contig}:{record.pos})")
            depth = int(pysam.depth("--reference",args.reference,"-q",str(1),"-Q",str(15),"-s", "-a", "-r", f"{record.contig}:{record.pos}-{record.pos}",args.cramfile,catch_stdout=True).strip().split('\t')[-1])
            record.samples[0]["DP"] = depth
            record.samples[0]["AD"] = (depth,0)
            record.samples[0]["AF"] = (0.0,)

        # Set filters.
        if record.samples[0]["AD"][1] < args.min_reads:
            record.filter.add("LowReads")

        if record.samples[0]["AF"][0] * 100 < args.vaf_soft_filter:
            record.filter.add("MinVAF")

        vcf_out.write(record)


    vcf_out.close()


if __name__ == "__main__":
    main()
