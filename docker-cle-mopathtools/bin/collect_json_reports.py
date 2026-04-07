#!/usr/bin/env python3

import argparse
import os
import sys
from natsort import natsort_keygen
import pandas as pd

__version__ = "1.1.0"

parser = argparse.ArgumentParser(description="Collect ChromoSeq Reports.")
parser.add_argument("files", type=str, nargs="+", help="List of input CSV files.")
parser.add_argument(
    "--outfile", "-o", type=str, help="Output file. If not provided, prints to stdout."
)
parser.add_argument(
    "--tables",
    "-t",
    type=str,
    default=None,
    help="Print individual tables for variants, svs, and qc.",
)
parser.add_argument(
    "--version", "-v", action="version", version="%(prog)s: " + __version__
)

args = parser.parse_args()

variants = pd.DataFrame()
svs = pd.DataFrame()
qc = pd.DataFrame()
outdf = pd.DataFrame()

for file in args.files:
    if not os.path.exists(file):
        print(f"Error: File '{file}' does not exist.")
        continue

    rep = pd.read_json(file)

    caseInfo = pd.DataFrame(data=rep["CASEINFO"].dropna()).to_dict()["CASEINFO"]
    df = pd.DataFrame(
        data=rep["VARIANTS"].dropna().to_dict()["data"],
        columns=rep["VARIANTS"].dropna().to_dict()["columns"],
    )
    df.insert(loc=0, column="name", value=caseInfo["name"])
    
    variants = pd.concat([variants.dropna(axis=1, how='all'), df.dropna(axis=1, how='all')], ignore_index=True)
    # make vaf numeric
    if not variants.empty and 'vaf' in variants.columns:
        variants['vaf'] = pd.to_numeric(variants['vaf'], errors='coerce')

    df = pd.DataFrame(
        data=rep["SV"].dropna().to_dict()["data"],
        columns=rep["SV"].dropna().to_dict()["columns"],
    )
    df.insert(loc=0, column="name", value=caseInfo["name"])

    svs = pd.concat([svs.dropna(axis=1, how='all'), df.dropna(axis=1, how='all')], ignore_index=True, axis=0)

    infoDf = (
        pd.DataFrame(data=rep["CASEINFO"].dropna())
        .transpose()
        .reset_index()
        .drop("index", axis=1)
    )
    qcDf = pd.DataFrame(
        data=rep["QC"]["ASSAY"]["data"], columns=rep["QC"]["ASSAY"]["columns"]
    )
    qcDf = qcDf[qcDf["qcmetric"] == 1][["metric", "value"]]
    
    qcDf = pd.concat([
            qcDf[qcDf['metric'].str.startswith('MAPPING/ALIGNING SUMMARY')],
            qcDf[qcDf['metric'].str.startswith('COVERAGE SUMMARY')]
                    .sort_values(
                    by=['metric'],
                    key=lambda col: natsort_keygen()(col) if col.name == 'metric' else col),
            qcDf[qcDf['metric'].str.startswith('CNV SUMMARY')],
            qcDf[qcDf['metric'].str.startswith('HAPLOTECT')]
        ],axis=0).set_index("metric")

    qc = pd.concat(
        [
            qc,
            pd.concat(
                [
                    pd.DataFrame(data=rep["CASEINFO"].dropna())
                    .replace(r"\s+\(!\)", "", regex=True)
                    .transpose()
                    .reset_index()
                    .drop("index", axis=1),
                    qcDf.transpose().reset_index().drop("index", axis=1),
                ],
                axis=1,
            ),
        ],
        ignore_index=True,
        axis=0,
    )

outdf = pd.merge(
    qc[
        [
            "version",
            "workflow_version",
            "name",
            "mrn",
            "DOB",
            "accession",
            "specimen",
            "exception",
            "date",
            "runid",
        ]
    ],
    variants[variants["category"] == "PASS"]
    .assign(
        passvariants=lambda df: df.apply(lambda row: f"{row['gene']}:{row['psyntax']}[{round(row['vaf'], 2)}%]", axis=1)
    )
    .groupby("name")["passvariants"]
    .agg(lambda x: ",".join(x))
    .reset_index(),
    on="name",
    how="left",
)
outdf = pd.merge(
    outdf,
    svs[svs["category"] == "RECURRENTSV"]
    .groupby("name")["known_genes"]
    .agg(lambda x: ",".join(x.dropna().astype(str)))
    .reset_index()
    .rename(columns={"known_genes": "recurrentsv"}),
    on="name",
    how="left",
)

outdf = pd.merge(
    outdf,
    svs[svs["category"] == "CNV"]
    .assign(
        cnvs=lambda df: df["psyntax"].str.replace("seq[GRCh38] ", "", regex=False)
        + df["known_genes"].fillna("").apply(lambda x: f"[{x}]" if x else "")
    )
    .groupby("name")["cnvs"]
    .agg(lambda x: ",".join(x.dropna().astype(str)))
    .reset_index(),
    on="name",
    how="left",
)

outdf = pd.merge(
    outdf,
    svs[svs["category"] == "OTHERSV"]
    .assign(
        othersvs=lambda df: df["psyntax"].str.replace("seq[GRCh38] ", "", regex=False)
        + df["known_genes"].fillna("").apply(lambda x: f"[{x}]" if x else "")
    )
    .groupby("name")["othersvs"]
    .agg(lambda x: ",".join(x.dropna().astype(str)))
    .reset_index(),
    on="name",
    how="left",
)

if args.tables is not None:
    variants.to_csv(args.tables + "_variants.tsv", sep="\t", index=False)
    svs.to_csv(args.tables + "_svs.tsv", sep="\t", index=False)
    qc.to_csv(args.tables + "_qc.tsv", sep="\t", index=False)
    with pd.ExcelWriter(args.tables + ".xlsx") as writer:
        outdf.to_excel(writer, sheet_name="summary", index=False)
        variants.to_excel(writer, sheet_name="variants", index=False)
        svs.to_excel(writer, sheet_name="svs", index=False)
        qc.to_excel(writer, sheet_name="qc", index=False)

if args.outfile:
    outdf.to_csv(args.outfile, sep="\t", index=False)
else:
    outdf.to_csv(sys.stdout, sep="\t", index=False)
