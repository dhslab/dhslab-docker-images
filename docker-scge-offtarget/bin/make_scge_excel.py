#!/usr/bin/env python3

import argparse
import json
import pandas as pd
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Generate Excel report for SCGE pipeline')
    parser.add_argument('--report_json', required=True, help='Path to report JSON file')
    parser.add_argument('--circos_plot', help='Path to Circos plot image')
    parser.add_argument('--cna_plot', help='Path to CNA plot image')
    parser.add_argument('--baf_plot', help='Path to BAF plot image')
    parser.add_argument('--indel_freq_plot', help='Path to Indel Frequency plot image')
    parser.add_argument('--off_targets_plot', help='Path to Off-targets plot image')
    parser.add_argument('--output', required=True, help='Output Excel filename')
    return parser.parse_args()

def add_image_to_sheet(worksheet, image_path, cell, x_scale=0.5, y_scale=0.5):
    """Helper to safely add an image if it exists."""
    if image_path and os.path.exists(image_path) and os.path.getsize(image_path) > 0:
        try:
            worksheet.insert_image(cell, image_path, {'x_scale': x_scale, 'y_scale': y_scale})
            return True
        except Exception as e:
            print(f"Warning: Could not insert image {image_path}: {e}")
    return False

def main():
    args = parse_args()

    # Load JSON data
    try:
        with open(args.report_json, 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error loading JSON: {e}")
        sys.exit(1)

    # Create Excel Writer
    # Check for available engines
    try:
        import xlsxwriter
        engine = 'xlsxwriter'
    except ImportError:
        try:
            import openpyxl
            engine = 'openpyxl'
        except ImportError:
            print("Error: Neither xlsxwriter nor openpyxl is installed.")
            sys.exit(1)

    writer = pd.ExcelWriter(args.output, engine=engine)
    workbook = writer.book

    # --- Sheet 1: Summary & Plots ---
    # Only xlsxwriter supports insert_image easily with this API. 
    # If openpyxl, image insertion is different. We will focus on xlsxwriter logic 
    # as it's standard for formatting. If openpyxl, we skip images or adapt.
    
    summary_ws = workbook.add_worksheet('Summary')
    
    # Formats (only works with xlsxwriter)
    if engine == 'xlsxwriter':
        bold = workbook.add_format({'bold': True})
        title_format = workbook.add_format({'bold': True, 'font_size': 14})
    else:
        bold = None
        title_format = None

    # Write Metadata
    if engine == 'xlsxwriter':
        summary_ws.write('A1', 'Somatic Editing Genome Report', title_format)
    else:
        summary_ws.append(['Somatic Editing Genome Report'])

    metadata = [
        ('Sample ID', data.get('sample_id', 'N/A')),
        ('Control Sample', data.get('metadata', {}).get('control_sample', 'N/A')),
        ('Transgene', data.get('transgene_description', 'N/A')),
        ('Tumor Mean Coverage', f"{data.get('metadata', {}).get('mean_coverage', {}).get('tumor', 'N/A')}x"),
        ('Normal Mean Coverage', f"{data.get('metadata', {}).get('mean_coverage', {}).get('normal', 'N/A')}x")
    ]

    row = 2
    for key, val in metadata:
        if engine == 'xlsxwriter':
            summary_ws.write(row, 0, key, bold)
            summary_ws.write(row, 1, val)
        else:
            summary_ws.append([key, val])
        row += 1

    # Place Images (xlsxwriter only for now)
    if engine == 'xlsxwriter':
        current_row = row + 2
        summary_ws.write(current_row, 0, 'On-Target Indel Frequency', bold)
        add_image_to_sheet(summary_ws, args.indel_freq_plot, f'A{current_row+2}')
        
        summary_ws.write(current_row, 8, 'Transgene Integrations (Circos)', bold)
        add_image_to_sheet(summary_ws, args.circos_plot, f'I{current_row+2}')

        current_row += 25 
        
        summary_ws.write(current_row, 0, 'Off-Target Sites', bold)
        add_image_to_sheet(summary_ws, args.off_targets_plot, f'A{current_row+2}')

        current_row += 25
        summary_ws.write(current_row, 0, 'Copy Number Analysis (CNA)', bold)
        add_image_to_sheet(summary_ws, args.cna_plot, f'A{current_row+2}')
        
        summary_ws.write(current_row, 8, 'B-Allele Frequency (BAF)', bold)
        add_image_to_sheet(summary_ws, args.baf_plot, f'I{current_row+2}')


    # --- Sheet 2: On-Target Variants ---
    transgene_data = data.get('tables', {}).get('on_target_sv_transgene', [])
    if transgene_data:
        if isinstance(transgene_data, list):
            df_on_target = pd.DataFrame(transgene_data)
        else:
            df_on_target = pd.DataFrame([transgene_data])
            
        cols_to_keep = ['Location', 'SYMBOL', 'Consequence', 'EXON', 'INTRON', 'HGVSc', 'HGVSp']
        existing_cols = [c for c in cols_to_keep if c in df_on_target.columns]
        if existing_cols:
            df_on_target = df_on_target[existing_cols]
        
        df_on_target.to_excel(writer, sheet_name='On-Target Variants', index=False)
    else:
        # Empty sheet
        pd.DataFrame({'Message': ['No on-target variants found']}).to_excel(writer, sheet_name='On-Target Variants', index=False)


    # --- Sheet 3: Off-Target Indels ---
    off_target_data = data.get('tables', {}).get('off_target_indels', [])
    if off_target_data:
        df_off = pd.DataFrame(off_target_data)
        df_off.to_excel(writer, sheet_name='Off-Target Indels', index=False)
    else:
        pd.DataFrame({'Message': ['No off-target indels found']}).to_excel(writer, sheet_name='Off-Target Indels', index=False)

    # --- Sheet 4: Targeted Mutations ---
    targeted_muts = data.get('tables', {}).get('targeted_gene_mutations', [])
    if targeted_muts:
        df_muts = pd.DataFrame(targeted_muts)
        df_muts.to_excel(writer, sheet_name='Targeted Gene Mutations', index=False)
    else:
        pd.DataFrame({'Message': ['No targeted gene mutation data available']}).to_excel(writer, sheet_name='Targeted Gene Mutations', index=False)

    writer.close()
    print(f"Successfully generated {args.output}")

if __name__ == "__main__":
    main()
