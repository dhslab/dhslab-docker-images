#!/usr/bin/env python3

import argparse
import itertools
import os
import re
import sys

import pandas as pd
import pyranges as pr
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


# calculates % abundance
def ratio2abundance(cn, nl, l2r):
    t = (nl * (2**l2r - 1)) / (cn - nl)
    return t

def segmentmean2abundance(cn, nl, sm):
    t = (nl * (sm - 1)) / (cn - nl)
    return t

def parse_vcf_header(header_string):
    # Split the header by lines
    header_lines = header_string.split("\n")

    # Filter lines that start with ##
    header_info = {}
    for line in header_lines:
        if line.startswith("##"):
            # Capture key and value in the header line
            match = re.match(r"##([^=]+)=(.*)", line)
            if match:
                key = match.group(1)
                value = match.group(2)
                # Check if the value is a complex structured format (like <ID=...,Description="...">)
                if value.startswith("<") and value.endswith(">"):
                    # Extract key-value pairs inside the angled brackets
                    attributes = re.findall(r'(\w+)=("[^"]*"|[^,>]+)', value[1:-1])
                    # Store as a dictionary
                    header_info[key] = dict(attributes)
                else:
                    # Otherwise, store the value as-is
                    header_info[key] = value
    return header_info

def vcf_to_pyranges(vcf_file):
    vcf = pysam.VariantFile(vcf_file)

    records = pd.DataFrame(
        columns=[
            "Chromosome",
            "Start",
            "End",
            "ID",
            "REF",
            "ALT",
            "CN",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
        ]
    )
    for record in vcf.fetch():

        # skip if its a REF call
        if record.alts is None:
            continue

        chrom = record.contig
        start = record.pos
        end = record.stop
        record_id = record.id
        ref = record.ref
        alt = ",".join(record.alts)
        qual = record.qual
        filter = ";".join(record.filter.keys())

        info_dict = dict(record.info)
        info_dict["set"] = "dragen"

        format_dict = {}
        for field in vcf.header.formats.keys():
            if field in record.samples[0]:
                format_dict[field] = record.samples[0][field]
            else:
                format_dict[field] = None

        # for now, HET records (which means the dragen models failed) will be filtered out.
        # the CN tag will also be set to 1 for DELs and 3 for DUPs
        if "HET" in info_dict or format_dict["CN"] is None:
            if "DEL" in alt:
                format_dict["CN"] = 1
            elif "DUP" in alt:
                format_dict["CN"] = 3
            elif alt == "<DEL>,<DUP>" or alt == "<LOH>":
                format_dict["CN"] = 2

        row_data = [
            chrom,
            start,
            end,
            record_id,
            ref,
            alt,
            format_dict["CN"],
            qual,
            filter,
            info_dict,
            format_dict,
        ]

        records.loc[len(records)] = row_data

    return pr.PyRanges(records)


def merge_records(df, overlap=0, size_filter=2000000):
    if df is None or df.empty:
        return pd.DataFrame()

    if "Cluster" in df.columns:
        df = df.drop(columns="Cluster")

    # unmerged records
    unmergedDf = pd.DataFrame(columns=df.columns)
    # merged records
    mergedDf = pd.DataFrame(columns=df.columns)

    # if there is only one record or they are all non-PASS then dont merge it.
    if (df.End.max() - df.Start.min()) < size_filter or len(
        df[df["FILTER"].isin(["PASS","segmentMean"])]
    ) == 0:
        return pd.DataFrame()

    elif df.shape[0] == 1:
        unmergedDf = df.copy()
    
    else:

        # 251122 (dhs) this omits centromere entries, which are included to merge
        # CNVs of the sample type across centromeres.
        to_merge = df[df["FILTER"].str.contains("PASS",na=False)].copy()
        unmergedDf = df[~df["FILTER"].str.contains("PASS",na=True)].copy()

        # sort to_merge by Start and End
        to_merge = to_merge.sort_values(by=["Start", "End"])

        # if there are records to merge then merge them
        if to_merge.shape[0] > 0:

            merged = {'Chromosome': to_merge['Chromosome'].iloc[0],
                      'Start': int(to_merge['Start'].iloc[0]),
                      'End': int(to_merge['End'].iloc[-1])
                    }

            # make svtype LOSS, GAIN, or CNLOH based on the ALT field
            if to_merge["ALT"].tolist()[0] == "<DEL>":
                merged["ID"] = f"DRAGEN:LOSS:{merged['Chromosome']}:{merged['Start']}-{merged['End']}"
            elif to_merge["ALT"].tolist()[0] == "<DUP>":
                merged["ID"] = f"DRAGEN:GAIN:{merged['Chromosome']}:{merged['Start']}-{merged['End']}"
            elif (to_merge["ALT"].tolist()[0] == "<DEL>,<DUP>" or to_merge["ALT"].tolist()[0] == "<LOH>"):
                if "CNLOH" in to_merge["ID"].tolist()[0]:
                    merged["ID"] = f"DRAGEN:CNLOH:{merged['Chromosome']}:{merged['Start']}-{merged['End']}"
                elif "GAINLOH" in to_merge["ID"].tolist()[0]:
                    merged["ID"] = f"DRAGEN:GAINLOH:{merged['Chromosome']}:{merged['Start']}-{merged['End']}"

            merged["REF"] = "N"
            merged["ALT"] = to_merge["ALT"].tolist()[0]
            merged["QUAL"] = pd.to_numeric(to_merge["QUAL"], errors="coerce").max()

            # if any of the records to merge contain PASS, then it PASSes.
            if to_merge[to_merge["FILTER"]=="PASS"].shape[0] > 0:
                merged["FILTER"] = "PASS"

            # likewise for records that only have segmentMean
            elif to_merge[to_merge["FILTER"]=="segmentMean"].shape[0] > 0:
                merged["FILTER"] = "segmentMean"

            else:
                merged["FILTER"] = ";".join(to_merge["FILTER"].to_list())

            # get format fields and abundance estimate. Abundance will be the mean of the abundance from ICHOR and 2 * the MAF from Dragen.
            infoDf = pd.DataFrame(to_merge["INFO"].tolist())
            formatDf = pd.DataFrame(to_merge["FORMAT"].tolist())

            info_dict = {
                "set": ",".join(set(infoDf["set"].dropna().tolist())),
                "SVTYPE": "CNV",
                "SVLEN": merged["End"] - merged["Start"] + 1,
            }

            if (
                merged["ALT"] == "<DEL>"
                or merged["ALT"] == "<DEL>,<DUP>"
                or merged["ALT"] == "<DUP>,<DEL>"
                or merged["ALT"] == "<LOH>"
            ):
                info_dict["REFLEN"] = info_dict["SVLEN"]
            else:
                info_dict["REFLEN"] = -info_dict["SVLEN"]

            if merged["ALT"] == "<DEL>":
                info_dict["SVLEN"] = -info_dict["SVLEN"]

            # Define your aggregation functions for each column
            aggregations = {
                "GT": lambda x: x.iloc[0],  # Take the first GT (or you could define a more complex logic here)
                "CN": "mean",
                "AS": "sum",  # Summing 'AS'
                "BC": "sum",  # Summing 'BC'
                "SD": "sum",  # Summing 'BC'
                "PE1": "sum",  # Summing the first element of 'PE' tuples
                "PE2": "sum",  # Summing the second element of 'PE' tuples
                "MCN": "mean",
                "CNQ": "mean",
                "MCNQ": "mean",
                "CNF": "mean",
                "MCNF": "mean",
                "SM": "mean",
                "MAF": "mean",
            }

            if "PE" in formatDf.columns:
                formatDf[["PE1", "PE2"]] = pd.DataFrame(
                    formatDf["PE"].dropna().tolist(),
                    index=formatDf["PE"].dropna().index,
                )

            for k in aggregations.keys():
                if k not in formatDf.columns:
                    formatDf[k] = None

            # Perform the aggregation
            formatDf = formatDf.agg(aggregations)
            # take the floor of the CN column
            formatDf["CN"] = int(formatDf["CN"])

            formatDf["PE"] = (formatDf["PE1"], formatDf["PE2"])
            format_dict = formatDf.drop(["PE1", "PE2"], axis=0).to_dict()

            mergedDf = pd.DataFrame([merged])
            mergedDf["INFO"] = [info_dict]
            mergedDf["FORMAT"] = [format_dict]

    # rescale Dragen's MAF to a AF for the unmerged records
    if not unmergedDf.empty:
        for rec in unmergedDf.to_dict(orient="records"):
            dfs_to_merge = [df.dropna(axis=1, how='all') for df in [mergedDf, pd.DataFrame([rec])] if not df.empty]
            mergedDf = pd.concat(dfs_to_merge, axis=0)
        
    if not mergedDf.empty:
        return mergedDf
    else:
        return pd.DataFrame()


def main():
    parser = argparse.ArgumentParser(description="Merge VCF records with overlap")
    parser.add_argument(
        "--vcf", "-vcf", type=checkfile, help="Path to the second input VCF file"
    )
    parser.add_argument(
        "--cytoband", "-c", type=checkfile, required=True, help="Bed file with cytoband intervals"
    )
    parser.add_argument(
        "--sex",
        "-x",
        required=True,
        choices=["male", "female"],
        help="Sex ('male' or 'female')",
    )
    parser.add_argument(
        "--min_size", "-s", type=int, default=5000000, help="Minimum CNA size"
    )
    parser.add_argument(
        "--merge_distance", "-d", type=int, default=5000000, help="Maximum distance between CNAs with same CN to merge"
    )
    parser.add_argument(
        "--min_hard_filter_size",
        "-S",
        type=int,
        default=1000000,
        help="Minimum hard-filtered CNA size. Smallers CNVs will be omitted.",
    )
    parser.add_argument(
        "--min_abund",
        "-f",
        type=float,
        default=10.0,
        help="Minimum CNA cell fraction to pass filter",
    )
    parser.add_argument(
        "--discard_abund", "-m", type=float, default=0.02, help="Discard abundance"
    )
    parser.add_argument(
        "--version", "-v", action="version", version="%(prog)s: " + __version__
    )

    args = parser.parse_args()

    # make new VCF and set headers
    inVcf = pysam.VariantFile(args.vcf)

    # add info tags
    if "set" not in inVcf.header.info.keys():
        inVcf.header.info.add(
            "set", 1, "String", "Callers that identified this variant"
        )

    #    if 'CYTOBANDS' not in inVcf.header.info.keys():
    #        inVcf.header.info.add("CYTOBANDS",1,'String','List of cytobands in order that overlap the variant')

    if "BN" not in inVcf.header.formats.keys():
        inVcf.header.formats.add(
            "BN", 1, "Integer", "IchorCNA bins spanning this variant"
        )

    if "CF" not in inVcf.header.formats.keys():
        inVcf.header.formats.add(
            "CF", 1, "Float", "Estimated tumor cell fraction of variant (%)"
        )

    # add filters
    if "MinCNVSize" not in inVcf.header.filters.keys():
        inVcf.header.filters.add(
            "MinCNVSize",
            None,
            None,
            "CNV does not meet minimum size of " + str(args.min_size) + " bp",
        )

    # add filters
    if "MinCNVAbundance" not in inVcf.header.filters.keys():
        inVcf.header.filters.add(
            "MinCNVAbundance",
            None,
            None,
            "CNV does not meet minimum tumor cell fraction of " + str(args.min_abund),
        )

    if "cnvHetLowQual" not in inVcf.header.filters.keys():
        inVcf.header.filters.add(
            "cnvHetLowQual",
            None,
            None,
            "Dragens HET calling was enabled for this variant and it may be inaccurate.",
        )

    if "centromereFilter" not in inVcf.header.filters.keys():
        inVcf.header.filters.add(
            "centromereFilter",
            None,
            None,
            ">50% of the CNV involves centromeric regions.",
        )

    # read in cytoband file and make a pyranges object
    cytobandPr = pr.read_bed(args.cytoband)
    centromeresPr = cytobandPr[cytobandPr.df.Score.str.contains("cen|var",regex=True)]

    # Create a new VCF writer with added header lines
    outVcfFile = sys.stdout
    outVcf = pysam.VariantFile(outVcfFile, "w", header=inVcf.header)

    vcfPr = vcf_to_pyranges(args.vcf)

    if not vcfPr.df.empty:
        # 251122 (dhs) Add centromere entries for all obserfved variant types.
        # This will make it so CNVs within merge_distance of a centromere get merged.
        unique_alts = vcfPr.df['ALT'].unique().tolist()
        centromeres_df = centromeresPr.df.drop(columns=["Name","Score"]).copy()
        centromeres_df['ALT'] = [unique_alts] * len(centromeres_df)
        centromeres_exploded_df = centromeres_df.explode('ALT')
        centromeresReplicatedPr = pr.PyRanges(centromeres_exploded_df)
        vcfPr = pr.concat([vcfPr, centromeresReplicatedPr])

    mergedCnvDf = vcfPr.cluster(by=['ALT'],slack=args.merge_distance).df

    # merge nearby/overlapping CNVs with similar event
    if not mergedCnvDf.empty:
        mergedCnvDf = mergedCnvDf.groupby('Cluster',group_keys=False).apply(lambda x: merge_records(x,overlap=0,size_filter=args.min_hard_filter_size),include_groups=False)
        if not mergedCnvDf.empty:
            mergedCnvDf = mergedCnvDf[mergedCnvDf['FILTER'].notna()]

    # remove CNVs that just involve the centromeres
    if not mergedCnvDf.empty:
        # remove CNVs that mostly involve centromeres
        mergedCnvPr = pr.PyRanges(mergedCnvDf).coverage(centromeresPr, fraction_col="CentromereOverlap")
        # if more than 50% of the CNV overlaps with a centromere, then add the centromereFilter to FILTER field
        mergedCnvDf = mergedCnvPr.df.copy()
        mergedCnvDf.loc[mergedCnvDf["CentromereOverlap"] >= 0.5, "FILTER"] = mergedCnvDf[mergedCnvDf["CentromereOverlap"] >= 0.5].apply(lambda x: "centromereFilter" if x['FILTER']=="PASS" else x["FILTER"] + ";centromereFilter",axis=1)

        mergedCnvDf = mergedCnvDf.drop(columns=["CentromereOverlap", "NumberOverlaps"])


    for item in mergedCnvDf.to_dict(orient="records"):
        nrec = outVcf.new_record()
        nrec.chrom = item["Chromosome"]
        nrec.pos = item["Start"]
        nrec.stop = item["End"]
        nrec.id = item["ID"]
        nrec.filter.clear()

        # add existing filters
        if "PASS" in item["FILTER"]:
            nrec.filter.add("PASS")

        else:
            filters = item["FILTER"].split(";")
            for v in filters:
                nrec.filter.add(v)

        # add new filters
        if item["End"] - item["Start"] < int(args.min_size):
            if "PASS" not in nrec.filter.keys():
                nrec.filter.add("MinCNVSize")
            else:
                nrec.filter.clear()
                nrec.filter.add("MinCNVSize")

        cf_maf = None
        cf_sm = None
        # Calculate CF from MAF first, if its present and not zero, then the segment mean linear copy number ratio first
        if "MAF" in item["FORMAT"].keys() and not pd.isna(item["FORMAT"]["MAF"]) and item["FORMAT"]["MAF"] > 0:
            cf_maf = round(
                item["FORMAT"]["MAF"] * 2 * 100, 2
            )

        if "SM" in item["FORMAT"].keys() and not pd.isna(item["FORMAT"]["SM"]) and not 'LOH' in item["ID"] and args.sex == "male" and item["FORMAT"]["CN"]!=1 and (item["Chromosome"] == "chrY" or item["Chromosome"] == "chrX"):
            cf_sm = round(
                segmentmean2abundance(item["FORMAT"]["CN"], 1, item["FORMAT"]["SM"]) * 100, 2
            )

        # for autosomes or females
        elif "SM" in item["FORMAT"].keys() and not pd.isna(item["FORMAT"]["SM"]) and not 'LOH' in item["ID"] and item["FORMAT"]["CN"]!=2:
            cf_sm = round(
                segmentmean2abundance(item["FORMAT"]["CN"], 2, item["FORMAT"]["SM"]) * 100, 2
            )

        # if difference between cf_sm and cf_maf is >0.25, then use cf_sm, otherwise take minimum. 
        valid_cfs = [v for v in [cf_maf, cf_sm] if v is not None and v >= 0]
        if len(valid_cfs) > 1 and abs(valid_cfs[0] - valid_cfs[1]) > 0.25:
            item["FORMAT"]["CF"] = cf_sm
        else:        
            item["FORMAT"]["CF"] = min(valid_cfs) if valid_cfs else None

        # sometimes the calculation is >100, so correct it.
        if item["FORMAT"]["CF"] is not None and item["FORMAT"]["CF"] > 100:
            item["FORMAT"]["CF"] = 99.0

        # If CF is not present or below the minimum abundance, then set filter
        if item["FORMAT"]["CF"] is None or item["FORMAT"]["CF"] < float(args.min_abund):
            if "PASS" not in nrec.filter.keys():
                nrec.filter.add("MinCNVAbundance")
            else:
                nrec.filter.clear()
                nrec.filter.add("MinCNVAbundance")

        nrec.ref = item["REF"]
        nrec.alts = (item["ALT"],)

        for k in inVcf.header.info.keys():
            if k == "END":
                continue
            if k in item["INFO"].keys():
                nrec.info[k] = item["INFO"].get(k)
            else:
                nrec.info[k] = None

        for k in inVcf.header.formats.keys():
            if k in item["FORMAT"] and not pd.isna(item["FORMAT"][k]):
                fmt = inVcf.header.formats.get(k)
                val = item["FORMAT"][k]

                if fmt.type == "Integer":
                    if fmt.number == 1:
                        nrec.samples[0][k] = int(val)
                    else:
                        # Cast each element to native Python int
                        nrec.samples[0][k] = tuple(int(v) for v in val)
                else:
                    nrec.samples[0][k] = val
            else:
                nrec.samples[0][k] = tuple(
                    itertools.repeat(None, inVcf.header.formats.get(k).number)
                )

        outVcf.write(nrec)

    outVcf.close()


if __name__ == "__main__":
    main()
