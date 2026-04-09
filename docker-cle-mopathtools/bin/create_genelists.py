#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import requests
import json
import math
from natsort import natsort_keygen

def post_to_ensembl(endpoint, ids, expand=0):
    """
    Helper function to perform batched POST requests to Ensembl REST API.
    Handles chunking to respect API limits (max 1000 IDs per request).
    """
    base_url = "https://rest.ensembl.org"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    # Ensembl generic limit is often ~1000, keep it safe at 200
    CHUNK_SIZE = 200
    combined_results = {}
    
    # Remove duplicates and None
    unique_ids = list(set([x for x in ids if x]))
    
    for i in range(0, len(unique_ids), CHUNK_SIZE):
        chunk = unique_ids[i:i + CHUNK_SIZE]
        
        # Determine payload key based on endpoint
        if "symbol" in endpoint:
            payload = {"symbols": chunk, "expand": expand}
        else:
            payload = {"ids": chunk, "expand": expand}
            
        try:
            response = requests.post(base_url + endpoint, headers=headers, data=json.dumps(payload))
            if response.status_code == 200:
                combined_results.update(response.json())
            else:
                print(f"Warning: Batch request to {endpoint} failed (Code {response.status_code}).", file=sys.stderr)
        except Exception as e:
            print(f"Warning: Request error in batch fetch: {e}", file=sys.stderr)
            
    return combined_results

def parse_region_string(region_string):
    """
    Parses a genomic region string in the format 'chrom:start-end'.
    """
    if not region_string:
        raise ValueError("Region string is empty.")

    clean_region = region_string.replace(',', '')
    
    if ':' not in clean_region:
        raise ValueError(f"Invalid format: '{region_string}'. Expected 'chrom:start-end'.")
    
    chrom, coords = clean_region.rsplit(':', 1)
    
    if not chrom:
        raise ValueError(f"Invalid format: '{region_string}'. Chromosome name is empty.")
    
    if '-' not in coords:
        raise ValueError(f"Invalid format: '{region_string}'. Expected 'start-end'.")
    
    try:
        start_str, end_str = coords.split('-')
        start = int(start_str)
        end = int(end_str)
    except ValueError:
        raise ValueError(f"Invalid coordinates in '{region_string}'.")
        
    if start > end:
        raise ValueError(f"Invalid coordinates: Start ({start}) cannot be greater than End ({end}).")
    
    return chrom, start, end

def generate_transcript_bed(gene_id, gene_name, start, end, transcript_data, region=None):
    """
    Generates a single BED-like format record for the entire transcript.
    """
    chrom = ""
    strand = ""
    transcript_id = ""
    trx_start = ""
    trx_end = ""
    type = ""
    
    if not transcript_data and region is not None:
        chrom, start, end = parse_region_string(region)
        transcript_id = gene_id
        type = "transcript_region"
        trx_start = start
        trx_end = end
        strand = "0"

    else:
        # Get transcript-level information
        chrom = "chr" + transcript_data["seq_region_name"]
        strand_val = transcript_data["strand"]
        strand = "1" if strand_val == 1 else "-1"
        transcript_id = transcript_data['id']
        type = "transcript"
        if region:
            region_chrom, start, end = parse_region_string(region)
            type = "transcript_region"

        # Safely access translation data
        translation_data = transcript_data.get("Translation", {})
        trx_start = translation_data.get("start", start)
        trx_end = translation_data.get("end", end)

    info = f"sv|{type}|{gene_name}|{gene_id}|{transcript_id}|{trx_start}|{trx_end}|{strand}"

    bed_line_parts = [chrom, start-1, end, gene_name, info]
    return ["\t".join(map(str, bed_line_parts))]

def generate_exon_region_bed(gene_id, gene_name, transcript_data, region=None):
    """
    Generates BED format records for a region relative to a transcript.
    """
    bed_records = []

    if not region:
        return bed_records
        
    chrom, start, end = parse_region_string(region)

    translation_data = transcript_data.get("Translation", {})
    coding_start = translation_data.get("start", start)
    coding_end = translation_data.get("end", end)
    strand = transcript_data.get("strand", "1")

    info = f"gene|exon_region|{gene_name}|{gene_id}|{transcript_data['id']}|{coding_start}|{coding_end}|{strand}"

    bed_line = ["chr" + transcript_data["seq_region_name"], start-1, end, gene_name, info]
    bed_records.append("\t".join(map(str, bed_line)))

    return bed_records

def generate_exon_bed(gene_id, gene_name, transcript_data, exon_list=None):
    """
    Generates BED format records for exons with 2 bp flanking sequence.
    """
    bed_records = []
    exons = transcript_data.get("Exon", [])
    if not exons:
        return bed_records
        
    for exon_number, exon in enumerate(exons, start=1):

        if exon_list and exon_number not in exon_list:
            continue

        start = int(exon["start"]) - 3
        end = int(exon["end"]) + 2
        strand = "1" if exon["strand"] == 1 else "-1"

        translation_data = transcript_data.get("Translation", {})
        coding_start = translation_data.get("start", start + 3)
        coding_end = translation_data.get("end", end - 2)

        if ((strand == "1" and start > coding_end) or
            (strand == "-1" and end < coding_start)): 
            continue

        if strand == "1" and end > coding_end: 
            end = coding_end + 1
        
        if strand == "-1" and start < coding_start:
            start = coding_start - 1

        info = f"gene|exon_{exon_number}|{gene_name}|{gene_id}|{transcript_data['id']}|{coding_start}|{coding_end}|{strand}"

        bed_line = ["chr" + exon["seq_region_name"], start, end, gene_name, info]
        bed_records.append("\t".join(map(str, bed_line)))

    return bed_records

def get_codon_range_for_region(exons, cds_start, cds_end, region_start, region_end, strand):
    # Sort exons in transcript order
    exons_sorted = sorted(exons, key=lambda e: e['start'], reverse=(strand == -1))
    cds_bases = []
    for exon in exons_sorted:
        # Get overlap of exon with CDS
        exon_start = max(exon['start'], cds_start)
        exon_end = min(exon['end'], cds_end)
        if exon_start > exon_end:
            continue
        # For -1 strand, walk from exon_end to exon_start
        if strand == -1:
            exon_bases = list(range(exon_end, exon_start - 1, -1))
        else:
            exon_bases = list(range(exon_start, exon_end + 1))
        cds_bases.extend(exon_bases)
    # Map region bases to CDS positions
    region_cds_positions = []
    for i, base in enumerate(cds_bases):
        if region_start <= base <= region_end or region_end <= base <= region_start:
            region_cds_positions.append(i + 1)  # 1-based
    if not region_cds_positions:
        return None
    first_codon = math.ceil(min(region_cds_positions) / 3)
    last_codon = math.ceil(max(region_cds_positions) / 3)
    return first_codon, last_codon

def generate_codon_region_bed(gene_id, gene_name, transcript_data, region=None):
    """
    Generates a BED record for the codon range overlapping a specified region of a gene.
    Returns a BED record with a string like 'codonX[-Y]' in the info field.
    """
    bed_records = []
    if not region or not transcript_data:
        return bed_records

    chrom, region_start, region_end = parse_region_string(region)
    translation_data = transcript_data.get("Translation", {})
    coding_start = translation_data.get("start")
    coding_end = translation_data.get("end")
    strand = transcript_data.get("strand", "1")
    exons = transcript_data.get("Exon", [])
    transcript_id = transcript_data.get("id", "NA")

    codon_range = get_codon_range_for_region(exons, coding_start, coding_end, region_start, region_end, strand)
    if not codon_range:
        return bed_records

    first_codon, last_codon = codon_range

    if first_codon == last_codon:
        codon_str = f"codon{first_codon}"
    else:
        codon_str = f"codons{first_codon}-{last_codon}"

    info = f"hotspot|{codon_str}|{gene_name}|{gene_id}|{transcript_id}|{coding_start}|{coding_end}|{strand}"
    bed_line = ["chr" + transcript_data["seq_region_name"], region_start-1, region_end, gene_name, info]
    bed_records.append("\t".join(map(str, bed_line)))
    return bed_records

def generate_other_region_bed(gene_id, region, target_type):
    """
    Generates BED format records for a region relative to a transcript.
    """
    bed_records = []

    if not region:
        return bed_records
        
    chrom, start, end = parse_region_string(region)

    info = f"{target_type}|{gene_id}|.|.|.|.|.|."

    bed_line = [chrom, start-1, end, gene_id, info]
    bed_records.append("\t".join(map(str, bed_line)))

    return bed_records


def parse_targets_from_args(args):
    """
    Normalize inputs into a list of target dictionaries and performs
    BATCH FETCHING from Ensembl to populate 'gene_data' and 'transcript_data'.
    """
    targets = []

    # --- 1. Construct Initial Target List ---
    
    # Process CLI Gene Identifiers
    if args.gene_identifier:
        if (args.region or args.exons) and len(args.gene_identifier) > 1:
            print(f"Error: Region or exons specified via CLI for multiple genes. This isn't allowed.", file=sys.stderr)
            sys.exit(1)

        exon_list = []
        if args.exons:
            try:
                exon_list = [int(x) for x in args.exons.split(",")]
            except ValueError:
                print(f"Error: Invalid CLI exon format '{args.exons}'", file=sys.stderr)
                sys.exit(1)

        for gid in args.gene_identifier:
            # If region provided with generic identifiers, treat as region lookup
            if args.region:                
                targets.append({
                    'lookup_id': gid,
                    'lookup_type': 'region',
                    'exons': exon_list,
                    'region': args.region,
                    'target_type': args.target_type
                })
            else:
                targets.append({
                    'lookup_id': gid,
                    'lookup_type': 'transcript' if args.transcript else 'gene',
                    'exons': exon_list,
                    'region': args.region,
                    'target_type': args.target_type
                })

    # Process Input File
    if args.file:
        try:
            input_df = pd.read_csv(args.file, sep='\t')
            input_df.columns = [c.lower() for c in input_df.columns]
            
            required_cols = {'gene_name', 'gene_id', 'transcript_id', 'target'}
            if not any(col in input_df.columns for col in required_cols) or not 'target' in input_df.columns:
                print(f"Error: Input file must contain at least one of: {required_cols} and a target column.", file=sys.stderr)
                sys.exit(1)

            input_df = input_df.astype(object).where(pd.notnull(input_df), None)

            for _, row in input_df.iterrows():
                t_id = row.get('transcript_id')
                g_id = row.get('gene_id')
                g_name = row.get('gene_name')
                target_type = row.get('target')
                                    
                if t_id and g_name and t_id == g_name:
                    lookup_id = g_name
                    lookup_type = 'region'
                elif g_id and g_name and g_id == g_name:
                    lookup_id = g_name
                    lookup_type = 'region'
                elif target_type == 'sv' and g_id:
                    lookup_id = g_id
                    lookup_type = 'gene'
                elif t_id:
                    lookup_id = t_id
                    lookup_type = 'transcript'
                elif g_id:
                    lookup_id = g_id
                    lookup_type = 'gene'
                elif g_name:
                    lookup_id = g_name
                    lookup_type = 'gene'
                else:
                    lookup_id = None
                    lookup_type = None

                ex_str = row.get('exons', None)
                row_exons = []
                if ex_str:
                    try:
                        row_exons = [int(float(x)) for x in str(ex_str).split(",")]
                    except ValueError:
                        print(f"Warning: malformed exon list '{ex_str}' for {lookup_id}, ignoring exons.", file=sys.stderr)

                targets.append({
                    'lookup_id': lookup_id,
                    'lookup_type': lookup_type,
                    'exons': row_exons,
                    'region': row.get('region'),
                    'hotspot': row.get('hotspot'),
                    'target_type': target_type
                })

        except FileNotFoundError:
            print(f"Error: Input file '{args.file}' not found.", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"Error parsing input file: {e}", file=sys.stderr)
            sys.exit(1)

    # --- 2. Batch Fetch Data from Ensembl ---
    
    # Bucketize identifiers
    transcript_ids = set()
    gene_ids = set()
    gene_symbols = set()
    
    for t in targets:
        # Initialize data placeholders
        t['gene_data'] = None
        t['transcript_data'] = None
        
        lid = t['lookup_id']
        ltype = t['lookup_type']
        
        if ltype == 'transcript':
            transcript_ids.add(lid)
        elif ltype == 'gene':
            # Heuristic: Ensembl IDs start with ENS...
            if lid.startswith('ENS'):
                gene_ids.add(lid)
            else:
                gene_symbols.add(lid)
    
    # Perform Batch Requests
    fetched_transcripts = {}
    fetched_genes = {}
    
    # 1. Fetch Transcripts (Direct) - Need expand=1 for exons
    if transcript_ids:
        print(f"Batch fetching {len(transcript_ids)} transcript IDs...", file=sys.stderr)
        fetched_transcripts.update(post_to_ensembl("/lookup/id", list(transcript_ids), expand=1))

    # 2. Fetch Genes (ID) - No expand needed yet (just need canonical id)
    if gene_ids:
        print(f"Batch fetching {len(gene_ids)} gene IDs...", file=sys.stderr)
        fetched_genes.update(post_to_ensembl("/lookup/id", list(gene_ids), expand=0))

    # 3. Fetch Genes (Symbol)
    if gene_symbols:
        print(f"Batch fetching {len(gene_symbols)} gene symbols...", file=sys.stderr)
        fetched_genes.update(post_to_ensembl("/lookup/symbol/homo_sapiens", list(gene_symbols), expand=0))

    # --- 3. Determine Secondary Lookups (Canonical Transcripts) ---
    canonical_transcript_ids = set()
    
    # Map back to targets to find which canonical transcripts we need
    for t in targets:
        lid = t['lookup_id']
        ltype = t['lookup_type']
        
        if ltype == 'gene':
            # Retrieve the gene object from our batch results
            g_obj = fetched_genes.get(lid)
            if g_obj:
                t['gene_data'] = g_obj
                # Identify canonical transcript
                if 'canonical_transcript' in g_obj:
                    # Canonical ID might look like ENST00000.1, split version
                    canon_id = g_obj['canonical_transcript'].split('.')[0]
                    canonical_transcript_ids.add(canon_id)
    
    # 4. Fetch Canonical Transcripts
    if canonical_transcript_ids:
        # Filter out ones we already fetched in step 1 to save bandwidth
        needed = [x for x in canonical_transcript_ids if x not in fetched_transcripts]
        if needed:
            print(f"Batch fetching {len(needed)} canonical transcripts...", file=sys.stderr)
            fetched_transcripts.update(post_to_ensembl("/lookup/id", needed, expand=1))

    # --- 4. Final Assignment to Targets ---
    for t in targets:
        lid = t['lookup_id']
        ltype = t['lookup_type']
        
        if ltype == 'transcript':
            t['transcript_data'] = fetched_transcripts.get(lid)
            # Try to populate gene_data from parent if missing? 
            # (Optional, but main loop handles display_name from transcript data)
            
        elif ltype == 'gene':
            g_obj = t['gene_data'] # populated earlier
            if g_obj and 'canonical_transcript' in g_obj:
                canon_id = g_obj['canonical_transcript'].split('.')[0]
                t['transcript_data'] = fetched_transcripts.get(canon_id)

    return targets

def main():
    parser = argparse.ArgumentParser(
        description="Generate BED file for coding exons or transcripts with flanking sequence."
    )
    parser.add_argument("gene_identifier", nargs="*", help="One or more gene names or Ensembl gene/transcript IDs")
    parser.add_argument("-t", "--transcript", action="store_true", default=False, help="Treat CLI input identifiers as transcript IDs")
    parser.add_argument("-e", "--exons", help="Comma-separated list of exon numbers")
    parser.add_argument("-r", "--region", help="Interval to target")
    parser.add_argument("-T", "--target-type", default="gene", help="Type of target for custom region handling (e.g. sv, gene, hotspot)")
    parser.add_argument("-f", "--file", help="Tab-delimited file with headers: gene_name, gene_id, transcript_id, exons, region, target_type")
    parser.add_argument("-a", "--add-to-file", help="Add entries to an existing BED file")
    parser.add_argument("-o", "--outfile", help="Output BED file name")

    args = parser.parse_args()

    if not args.gene_identifier and not args.file:
        parser.print_help()
        print("\nError: You must provide either gene identifiers via CLI or an input file via -f.", file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame(columns=["chrom", "start", "end", "name", "info"])

    if args.add_to_file:
        try:
            df = pd.read_csv(args.add_to_file, sep="\t", header=None, names=["chrom", "start", "end", "name", "info"])
        except FileNotFoundError:
            print(f"Error: File '{args.add_to_file}' not found", file=sys.stderr)
            sys.exit(1)

    targets = parse_targets_from_args(args)
    all_bed_records = []

    for target in targets:
        identifier = target['lookup_id']
        lookup_type = target['lookup_type']
        target_type = target['target_type']
        gene_start = None
        gene_end = None
        
        # Retrieve pre-fetched data
        gene_data = target.get('gene_data')
        transcript_data = target.get('transcript_data')

        try:
            gene_name = "NA"
            gene_id = "NA"
            
            # --- Logic to extract Gene Name/ID/Coords based on what data we have ---
            if lookup_type == "region":
                gene_id = identifier
                gene_name = identifier
                transcript_data = None
                region_chrom, gene_start, gene_end = parse_region_string(target['region'])

            elif lookup_type == 'transcript':
                if not transcript_data:
                    raise ValueError(f"Transcript data not found for {identifier}")
                
                gene_id = transcript_data.get("Parent", "NA")
                gene_name = transcript_data.get("display_name", "NA").split("-")[0]
                
                if target_type == 'sv':
                    # SV mode needs gene coordinates, but we only looked up transcript.
                    # We might lack gene_start/end if we didn't fetch the gene object separately.
                    # Fallback: use transcript start/end or fetch gene specifically if needed.
                    # (Current batch logic didn't fetch gene obj for transcript inputs to save time, 
                    # assuming transcript data is sufficient for most BED needs).
                    pass 

            elif lookup_type == 'gene':
                if not gene_data:
                    raise ValueError(f"Gene data not found for {identifier}")
                
                gene_id = gene_data["id"]
                gene_name = gene_data["display_name"]
                gene_start = gene_data.get("start")
                gene_end = gene_data.get("end")

                if not transcript_data:
                     # Check if it was because of missing canonical
                     if "canonical_transcript" not in gene_data:
                         raise ValueError(f"No canonical transcript found for gene '{gene_name}'")
                     else:
                         raise ValueError(f"Canonical transcript data failed to fetch for '{gene_name}'")
            
            # --- Generate BED ---
            bed_records = []
            if target_type == 'sv':
                # For SVs, we often want the gene start/end, which we have from gene_data
                bed_records = bed_records + generate_transcript_bed(gene_id, gene_name, gene_start, gene_end, transcript_data, target['region'])

            elif target_type == 'gene':
                if target['region']:
                    for sub_region in target['region'].split(","):
                        bed_records = bed_records + generate_exon_region_bed(gene_id, gene_name, transcript_data, sub_region)

                if target['exons']:
                    bed_records = bed_records + generate_exon_bed(gene_id, gene_name, transcript_data, target['exons'])

                # If neither region nor exons specified for 'gene' mode, usually we default to all exons?
                if not target['region'] and not target['exons']:
                    bed_records = bed_records + generate_exon_bed(gene_id, gene_name, transcript_data, None)

                if target['hotspot']:
                    for sub_region in target['hotspot'].split(","):
                        bed_records = bed_records + generate_codon_region_bed(gene_id, gene_name, transcript_data, sub_region)
            
            else:
                for sub_region in target['region'].split(","):
                    bed_records = bed_records + generate_other_region_bed(gene_id, sub_region, target_type)

            for record in bed_records:
                all_bed_records.append(record.split("\t"))

        except ValueError as e:
            print(f"Error processing '{identifier}': {e}", file=sys.stderr)
        except Exception as e:
            print(f"Unexpected error processing '{identifier}': {e}", file=sys.stderr)

    if all_bed_records:
        new_df = pd.DataFrame(all_bed_records, columns=["chrom", "start", "end", "name", "info"])
        df = pd.concat([df, new_df], ignore_index=True)

    if not df.empty:
        natsort_key = natsort_keygen()
        df['start'] = pd.to_numeric(df['start'])
        df = df.sort_values(
            by=["chrom", "start"], 
            key=lambda col: col.map(natsort_key) if col.name == 'chrom' else col
        ).reset_index(drop=True)

        if args.outfile:
            df.to_csv(args.outfile, sep="\t", index=False, header=False)
        else:
            df.to_csv(sys.stdout, sep="\t", index=False, header=False)

if __name__ == "__main__":
    main()