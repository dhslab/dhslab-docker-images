import csv
import sys

def tsv_to_vcf(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\\t')
        
        # Write VCF header
        outfile.write("##fileformat=VCFv4.2\\n")
        outfile.write("#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n")
        
        for row in reader:
            chrom = row[0]
            pos = row[1]
            # VCF is 1-based, so we assume the input is 1-based start
            
            # Create a simple representation for structural variants
            ref_allele = 'N'  # Placeholder
            alt_allele = '<INS>' # Placeholder for insertion
            
            # ID can be a dot or can be constructed if needed
            vcf_id = '.'
            
            # QUAL and FILTER are often placeholders
            qual = '.'
            filter_status = 'PASS'
            
            # Info can be constructed from the rest of the TSV
            info = f"END={row[2]};SVTYPE=INS" # Example
            
            outfile.write(f"{chrom}\\t{pos}\\t{vcf_id}\\t{ref_allele}\\t{alt_allele}\\t{qual}\\t{filter_status}\\t{info}\\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python tsv_to_vcf.py <input.tsv> <output.vcf>")
        sys.exit(1)
    
    input_tsv = sys.argv[1]
    output_vcf = sys.argv[2]
    tsv_to_vcf(input_tsv, output_vcf)
