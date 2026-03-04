#!/usr/bin/env python3

import argparse
import glob
import gzip
import os
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter

import pandas as pd

__version__ = "1.0.0"

"""
This script makes a custom-formatted fastq_list.csv file or read1 and read2 FASTQ files. Features include:

-Returns a fastq_list.csv file in this format:

RGID,RGSM,RGLB,Lane,RGPL,Read1File,Read2File

-Generates tags with this format:

RGID: <flowcell>.<i7index>.<i5index>.<lane>
RGSM: <sample id> (passed as argument)
RGLB: <sample id>.<i7index>.<i5index>
Lane: <lane number>
RGL: <runid>.<instrument><side>.<flowcell>.<flowcelltype>.<flowcell_lot>.<reagent_lot>.<runrecipe>
Read1File: <read1 file path>
Read2File: <read2 file path>

-Obtains RGPL info from either a RunParameters.xml file or the reads (info from reads is incomplete)

"""

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Prepare fastq_list.csv file from a passed list file or reads. If a Runparameters.xml file is passed, additional metadata is added."
    )
    parser.add_argument("-i", "--id", type=str, required=True, help="Sample ID")
    parser.add_argument("-1", "--read1", type=check_file, help="Path to read1")
    parser.add_argument("-2", "--read2", type=check_file, help="Path to read2")
    parser.add_argument(
        "-r", "--runinfo", type=str, help="Path to Illumina RunParameters.xml file"
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s: " + __version__
    )

    return parser.parse_args()


def check_file(file_path):
    """Check if a file exists at the given path."""
    if not os.path.exists(file_path):
        raise argparse.ArgumentTypeError(f"The file {file_path} does not exist.")
    return file_path


# NOTE(dhs): get flowcell, lane, read length from reads to make runinfo
def make_runinfo_from_read(readpath):
    # open read 1 fastq and read header and get seq length
    readName = ""
    with gzip.open(readpath, "rt") as gz_file:
        readName = gz_file.readline().strip()

    # if indexes are present in the read header, get the first 100
    # and find the most common one (to account for mismatches/errors in index read)
    parts = readName.split(":")
    index1 = "UNKNOWN"
    index2 = "UNKNOWN"
    index1len = "?"
    index2len = "?"
    indexes = readName.split(" ")
    readlen = "?"
    indexlist = []
    seqlist = []
    with gzip.open(readpath, "rt") as file:
        for i, line in enumerate(file):
            if i % 4 == 0:  # Read names are on every 4th line starting from 0
                read_name = line.strip()                
                indexes = read_name.split(" ")
                if len(indexes) > 1:
                    index = indexes[1].split(":")[-1]
                    if not 'N' in index:
                        indexlist = indexlist + [index]

            if i > 0 and i % 1 == 0:
                seqlist.append(len(line.strip()))

            if i > 40000 or len(indexlist) > 500:  # Stop after reading the first 500 indexes or 10,000 records
                break

    if len(indexlist) > 0:
        counter = Counter(indexlist)
        most_common_index, _ = counter.most_common(1)[0]
        index1, index2 = most_common_index.split("+")

    counter = Counter(seqlist)
    readlen, _ = counter.most_common(1)[0]

    # if indexes are ?, try to recover from the read name
    if index1 == "UNKNOWN":
        pattern = "[ACTG]{8,}"
        matches = re.findall(pattern, os.path.basename(readpath))
        if len(matches) == 2:
            index1, index2 = matches
            index1len = f"{len(index1)}?"
            index2len = f"{len(index2)}?"
        elif len(matches) == 1:
            index1 = matches
            index2 = matches

    # Make run info dict
    runinfo = {
        "RunId": f"RUN_{parts[0][1:]}_{str(int(parts[1])).zfill(4)}_{parts[2]}",
        "Flowcell": parts[2],
        "Lane": parts[3],
        "Instrument": parts[0][1:],
        "Read1Cycles": f"{readlen}",
        "Index1Cycles": index1len,
        "Index1Reverse": "?",
        "Index2Cycles": index2len,
        "Index2Reverse": "?",
        "Read2Cycles": f"{readlen}",
        "FlowcellType": "?",
        "FlowcellLot": "?",
        "ReagentLot": "?",
        "Index1": index1,
        "Index2": index2,
    }

    return runinfo

def make_runinfo_from_run_parameters(filepath):
    runinfo = {}

    tree = ET.parse(filepath)
    root = tree.getroot()
    runinfo = {}
    runinfo["RunId"] = root.find("RunId").text
    for el in root.findall(".//PlannedReads")[0].findall("Read"):
        runinfo[el.attrib["ReadName"] + "Cycles"] = int(el.attrib["Cycles"])

    runinfo["FlowcellType"] = root.find("RecipeName").text.split(" ")[0]
    runinfo["Instrument"] = root.find("InstrumentSerialNumber").text
    runinfo["Side"] = root.find("Side").text
    # fine the serial number of the flowcell:
    for el in root.findall("ConsumableInfo")[0].findall("ConsumableInfo"):
        el2 = el.find("Type")
        if el2 is not None:
            if el2.text == "FlowCell":
                runinfo["Flowcell"] = el.find("SerialNumber").text
                runinfo["FlowcellLot"] = el.find("LotNumber").text
            elif el2.text == "Reagent":
                runinfo["ReagentLot"] = el.find("LotNumber").text

    runinfo["Instrument"] = "".join(
        [runinfo.get("Instrument") or "?", runinfo.get("Side") or "?"]
    )

    return runinfo


def main():
    args = parse_arguments()

    runinfo = {}

    if args.runinfo:
        runparams_path = args.runinfo
        if not os.path.exists(runparams_path):
            raise ValueError(
                "RunParameters.xml file not found."
        )

        runinfo = make_runinfo_from_run_parameters(runparams_path)

    read1 = args.read1
    read2 = args.read2

    if not runinfo or "Lane" not in runinfo or "Index1" not in runinfo:
        runinfo_read1 = make_runinfo_from_read(read1)
        runinfo_read2 = make_runinfo_from_read(read2)

        if (
            runinfo_read1["Lane"] != runinfo_read2["Lane"]
            or runinfo_read1["Flowcell"] != runinfo_read2["Flowcell"]
            or runinfo_read1["Index1"] != runinfo_read2["Index1"]
        ):
            print(f"Read 1 and 2 do not match\n\t{read1}\n\t{read2}")
            sys.exit(1)

        # add runinfo_read1 key/value pairs to runinfo, but dont overwrite
        for key, value in runinfo_read1.items():
            if key not in runinfo:
                runinfo[key] = value

    

    fqlistout = pd.DataFrame(
        columns=["RGID", "RGSM", "RGLB", "Lane", "RGPL", "Read1File", "Read2File"]
    )

    rgid = ".".join(
        [runinfo["Flowcell"], runinfo["Index1"], runinfo["Index2"], runinfo["Lane"]]
    )
    rglb = ".".join([args.id, runinfo["Index1"], runinfo["Index2"]])
    rgpl = (
        f"{runinfo['RunId']}."
        f"{runinfo['Instrument']}."
        f"{runinfo['Flowcell']}."
        f"{runinfo['FlowcellType']}."
        f"{runinfo['FlowcellLot']}."
        f"{runinfo['ReagentLot']}."
        f"{runinfo['Read1Cycles']}x"
        f"{runinfo['Index1Cycles']}x"
        f"{runinfo['Index2Cycles']}x"
        f"{runinfo['Read2Cycles']}"
    )

    fqlistout.loc[len(fqlistout)] = [
        rgid,
        args.id,
        rglb,
        runinfo["Lane"],
        rgpl,
        read1,
        read2,
    ]

    fqlistout.to_csv(f"{args.id}.fastq_list.csv", index=False)



if __name__ == "__main__":
    main()
