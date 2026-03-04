#!/usr/bin/env python3

import argparse
import os

import pandas as pd

__version__ = "2.0.0"

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Process spreadsheet file for samples that require demultiplexing, alignment, or summary analysis."
    )
    # add version
    parser.add_argument('-v', '--version', action='version', version='%(prog)s: ' + __version__)
    parser.add_argument(
        "--input_file", "-i", required=True, help="Path to the input CSV file"
    )
    parser.add_argument(
        "--run_directory", "-r", help="Path to run directory"
    )
    parser.add_argument(
        "--output_dir", "-o", help="Directory to save the output CSV files"
    )
    parser.add_argument(
        "--sample_info", "-s", help="Path to sample information CSV file."
    )
    parser.add_argument(
         "--procedure_name",
         "-p",
         default=None,
         help="Sample procedure name (e.g., MyeloSeq, ChromoSeq, etc)."
    )
    parser.add_argument(
        "--procedure_code",
        "-c",
        default=None,
        help="Two character sample procedure code (e.g., CS, MS, etc)."
    )
    return parser.parse_args()


def clean_sample_name(value: str) -> str:
    """Remove whitespace and replace spaces with underscores."""
    if isinstance(value, str):
        return value.strip().replace(" ", "_")
    return value


def save_df(dataframe: pd.DataFrame, output_dir: str, filename: str) -> None:
    """Save dataframe to a CSV file, applying a cleaning function if provided."""
    if not dataframe.empty:
        dataframe.loc[:, "id"] = dataframe["id"].apply(
            clean_sample_name
        )

        dataframe = dataframe.dropna(axis=1, how="all")
        dataframe = dataframe.loc[:, dataframe.ne("").any()]

        dataframe.to_csv(f"{output_dir}/{filename}", index=False)


def process_input(input_file: str, output_dir: str, sample_info: str, procedure_name: str, procedure_code: str, run_directory: str) -> None:
    """Process input file, merge with sample info, and create output files."""
    if procedure_code:
        procedure_code = f"{procedure_code}-"
    else:
        procedure_code = ""

    ext = input_file.split(".")[-1]
    df = pd.read_csv(input_file, sep="\t" if ext == "tsv" else ",",engine="python",on_bad_lines="error")

    id_col = "id" if "id" in df.columns else "Content_Desc"
    if id_col not in df.columns:
        raise ValueError("No 'id' or 'Content_Desc' column found in input file!")

    df.rename(columns={id_col: "id"}, inplace=True)

    # Normalize 'Exceptions' and 'exceptions' columns to 'exception' for consistency
    if "exceptions" in df.columns:
        df.rename(columns={"exceptions": "exception"}, inplace=True)
    elif "Exceptions" in df.columns:
        df.rename(columns={"Exceptions": "exception"}, inplace=True)

    if sample_info:
        sample_df = pd.read_csv(sample_info,dtype={"MRN": str})
        sample_df = sample_df[
            sample_df["Procedure"].str.contains(procedure_name, na=False, regex=True)
        ].loc[:, ["Accession", "gender", "MRN", "Patient DOB", "Part type"]]

        df["Accession"] = df["id"].astype(str).str.extract(r"(^W\d+-\d+)", expand=False)
        df = df.merge(sample_df, on="Accession", how="left")

        for col, new_col in [
            ("gender", "sex"),
            ("MRN", "mrn"),
            ("Patient DOB", "dob"),
            ("Part type", "specimen"),
        ]:
            if col in df.columns:
                if new_col in df.columns:
                    # Prioritize existing 'new_col' values, fill NaNs with 'col' values
                    df[new_col] = df[new_col].combine_first(df[col])
                    df.drop(columns=[col], errors="ignore", inplace=True)
                else:
                    # If 'new_col' doesn't exist, rename 'col'
                    df.rename(columns={col: new_col}, inplace=True)

        df["specimen"] = (
            df["specimen"]
            .fillna("")
            .str.strip()
            .str.replace(r"^[A-Z]\s*-\s*", "", regex=True)
            .replace("", pd.NA)
        )

        specimen_mappings = {"Blood": "PB", "Bone Marrow": "BM", "Bone Marrow-WML": "BMW", "Fresh Tissue": "CP"}
        df["specimen_short_name"] = df["specimen"]
        for key, value in specimen_mappings.items():
            df["specimen_short_name"] = df["specimen_short_name"].str.replace(
                key, value, regex=True
            )

        if df["Accession"].notna().any():
            df.loc[df["Accession"].notna(), "id"] = df.loc[
                df["Accession"].notna()
            ].apply(
                lambda row: (
                    f"{procedure_code}Research-"
                    if "exception" in row and pd.notna(row["exception"])
                    and "research" in str(row["exception"]).lower()
                    else f"{procedure_code}"
                )
                + "-".join(
                    str(row[key]).strip().replace(" ", "_")
                    for key in ["Accession", "mrn", "specimen_short_name"]
                    if pd.notna(row[key])
                ),
                axis=1,
            )

    rename_map = {
        "Content_Desc": "id",
        "Lane": "lanes",
        "Flowcell ID": "flowcell",
        "Index": "index",
        "Accession": "accession"
    }

    # merge fields and keep valid headers
    for src, dest in rename_map.items():
        # Check if the source column exists
        if src in df.columns:
            # If the destination column also exists, merge them
            if dest in df.columns:
                # Prioritize 'src' values, use 'dest' to fill in missing gaps
                df[dest] = df[src].combine_first(df[dest])
                # Drop the original source column
                df.drop(columns=[src], inplace=True)
            # Otherwise, it's a simple rename
            else:
                df.rename(columns={src: dest}, inplace=True)
            
    # modify sex column
    if "sex" in df.columns:
        df["sex"] = df["sex"].str.lower().replace({"m": "male", "f": "female"})
        df.dropna(axis=1, how="all", inplace=True)

    # add run directory
    if run_directory:
        df["run_directory"] = run_directory

    # remove columns if not valid header
    valid_headers = [
        "id",
        "sex",
        "mrn",
        "all_mrns",
        "accession",
        "dob",
        "specimen",
        "exception",
        "flowcell",
        "run_directory",
        "lanes",
        "index",
        "read1",
        "read2",
        "demux_path",
        "fastq_list",
        "fastq_id",
        "cram",
        "bam",
        "dragen_path",
    ]
    df = df.loc[:, [col for col in df.columns if col in valid_headers]]
    
    # Summary analysis samples. This supercedes all other inputs--samples with dragen_path will only be analyzed.
    if "dragen_path" in df.columns:
        dragen_df = df.dropna(subset=["dragen_path"]).dropna(axis=1, how="all")
        save_df(dragen_df, output_dir, "summary_analysis_samples.csv")
        df.drop(columns=["dragen_path"], inplace=True)

    # Alignment samples. Note this includes samples to demux as well as already demuxed fastqs and cram/bam for realignment.
    alignment_cols = ["bam", "cram", "read1", "read2", "demux_path", "fastq_list","index"]
    alignment_df = df.dropna(
        subset=[col for col in alignment_cols if col in df], thresh=1
    )
    save_df(alignment_df, output_dir, "alignment_samples.csv")

    # Demux samples
    demux_columns = ["id", "flowcell", "lanes", "index"]
    demux_df = df.dropna(subset=[col for col in demux_columns if col in df], thresh=3)
    save_df(demux_df, output_dir, "demux_samples.csv")

    # all sample metadata
    meta_columns = ["id","sex","mrn","all_mrns","accession","dob","specimen","exception"]
    meta_df = df[[col for col in meta_columns if col in df]].drop_duplicates()
    save_df(meta_df, output_dir, "sample_metadata.csv") 

def main() -> None:
    args = parse_arguments()

    output_dir = args.output_dir
    if not output_dir or not os.path.exists(output_dir):
        output_dir = os.getcwd()

    sample_info = args.sample_info
    if not args.sample_info or not os.path.exists(sample_info):
        sample_info = None

    process_input(args.input_file, output_dir, sample_info, args.procedure_name, args.procedure_code, args.run_directory)


if __name__ == "__main__":
    main()
