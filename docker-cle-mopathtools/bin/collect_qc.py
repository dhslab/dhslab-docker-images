#!/usr/bin/env python

import sys, os, re, csv, pysam, json, binascii, math, argparse
import sqlite3
import pandas as pd
import numpy as np
from time import gmtime, strftime
from natsort import natsort_keygen
from pathlib import Path

__version__ = '1.0.0'

def sanitize_for_json(obj):
    """Recursively replace inf and nan with None."""
    if isinstance(obj, dict):
        return {k: sanitize_for_json(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [sanitize_for_json(v) for v in obj]
    elif isinstance(obj, float):
        if math.isinf(obj) or math.isnan(obj):
            return None  # Will become null in JSON
    return obj

def parse_dragen_metrics_file(metrics_file):
    df = pd.read_csv(metrics_file,sep=',',names=['category','readgroup','metric','value','percent'])
    # add sequential index as column for sorting later
    df.insert(0,'index',range(len(df)))

    df = df[df['readgroup'].isna()].drop(columns='readgroup')
    df['metric'] = df['category'] + ': ' + df['metric']
    dfpct = df[df['percent'].notna()].copy()
    df = df.drop(columns='percent')
    dfpct['metric'] = dfpct['metric'].apply(lambda x: x + ' (%)')
    dfpct['value'] = dfpct['percent']
    dfpct = dfpct.drop(columns='percent')
    df = pd.concat([df,dfpct])

    # change metric with format [\d+x: inf) to >\d+x
    df.loc[df['metric'].str.contains(r'\[\s*\d+x: inf\)'), 'metric'] = df[df['metric'].str.contains(r'\[\s*\d+x: inf\)')]['metric'].str.replace(r'\[\s*(\d+)x: inf\)', r'>\1x', regex=True)

    # change 'PCT of' to "Percent of"
    df.loc[df['metric'].str.contains('PCT of'), 'metric'] = df[df['metric'].str.contains('PCT of')]['metric'].str.replace('PCT of', 'Percent of', regex=False)

    # if metric contains "Number of" and contains "(%)" replace "Number of" with "Percent"
    df.loc[df['metric'].str.contains('Number of') & df['metric'].str.contains('\(%\)'), 'metric'] = df[df['metric'].str.contains('Number of') & df['metric'].str.contains('\(%\)')]['metric'].str.replace('Number of', 'Percent', regex=False)

    # convert value to numeric where possible
    df['value'] = pd.to_numeric(df['value'], errors='coerce')

    # Duplicate metrics for large values with new units, keeping originals
    df_g = df[df['value'] > 5e9].copy()
    df_g['metric'] = df_g['metric'] + ' (G)'
    df_g['value'] = (df_g['value'] / 1e9).round(1)
    df_m = df[(df['value'] > 5e6) & (df['value'] <= 5e9)].copy()
    df_m['metric'] = df_m['metric'] + ' (M)'
    df_m['value'] = (df_m['value'] / 1e6).round(1)
    df = pd.concat([df, df_g, df_m], ignore_index=True)

    # order by index column
    df = df.sort_values(by=['index','metric']).drop(columns='index').reset_index(drop=True)

    return df[['category','metric','value']]

def parse_coverage_report(coverage_report_file,aggregate_keys=['gene','biomarker']):
    df = pd.read_csv(coverage_report_file, header=None, sep="\t")
    df.columns = df.iloc[0]
    df.columns = df.columns.str.replace('^#', '', regex=True)
    df = df.drop(df.index[0])
    df.fillna(0, inplace=True)
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    # Convert columns 5 onward to float in place to avoid dtype conflict
    for col in df.columns[5:]:
        df[col] = pd.to_numeric(df[col], errors='coerce')
        
    df.insert(4,'region',df['info'].str.split('|',expand=True).loc[:,1]) # parse region from info column
    df.insert(5,'region_type',df['info'].str.split('|',expand=True).loc[:,0]) # parse region from info 
    df.reset_index(drop=True, inplace=True)

    # if df['region'] contains any of the aggregate keys
    if aggregate_keys is not None and len(df['region_type'].isin(aggregate_keys)) > 0:
        for key in aggregate_keys:
            sdf = df[df['region_type']==key].copy()
            sdf['length'] = sdf['end']-sdf['start']+1
            pct_cov_columns = [x for x in sdf.columns if "pct_above" in x]
            sdf[pct_cov_columns] = sdf[pct_cov_columns].astype(float).div(100).mul(sdf['length'], axis=0)
            aggregate_funcs = {key: 'sum' for key in pct_cov_columns + ['total_cvg','length']}
            aggregate_funcs['min_cvg'] = 'min'
            aggregate_funcs['max_cvg'] = 'max'
            aggregate_funcs['start'] = 'min'
            aggregate_funcs['end'] = 'max'
            sdf = sdf.groupby(['gene']).agg(aggregate_funcs).reset_index()
            sdf[pct_cov_columns] = sdf[pct_cov_columns].div(sdf['length'], axis=0).mul(100)
            sdf[pct_cov_columns] = sdf[pct_cov_columns].round(1)
            sdf['mean_cvg'] = sdf['total_cvg'] / sdf['length']

            sdf['region'] = key
            sdf['region_type'] = key

            sdf = sdf.drop(columns='length')
            sdf = sdf.reindex(columns=df.columns)

            df = pd.concat([df,sdf],axis=0)

    df[['mean_cvg','min_cvg','max_cvg']] = df[['mean_cvg','min_cvg','max_cvg']].astype(int)

    return df

def parse_wgs_histogram(wgs_hist_file, coverage_depths=[1,5,10,20,30,50,100,200,500,1000]):
    dtype = {'Depth': str, 'Overall': str}
    histDf = (pd.read_csv(
                wgs_hist_file,
                sep=",",
                dtype=dtype
            )
            .iloc[:-1]  # Omits the last row (replaces the "2000+" string filter)
            .assign(
                Depth=lambda x: pd.to_numeric(x['Depth'], errors='coerce') + 1,
                Overall=lambda x: pd.to_numeric(x['Overall'], errors='coerce')
            )
            .assign(Fraction=lambda x: 100 - (x['Overall'].cumsum() / x['Overall'].sum()) * 100)
            .query('Depth in @coverage_depths')
            .sort_values('Depth', ascending=False)
    )

    histDf['category'] = "COVERAGE SUMMARY"
    histDf['metric'] = histDf.apply(lambda x: f"COVERAGE SUMMARY: Percent of genome with coverage >{x['Depth']}x",axis=1)
    histDf['value'] = histDf.apply(lambda x: round(x['Fraction'], 1), axis=1)

    return histDf[['category','metric','value']]

# parse haplotect loci file or dataframe into per-site genotypes
def pack_haplotect(haplotectlocidf=None,haplotect_file=None):
    df = pd.DataFrame()
    if haplotectlocidf is not None:
        df = haplotectlocidf[['chr','snp1','snp2','total_reads','haplotype_counts']].copy()
    elif haplotect_file is not None and Path(haplotect_file).is_file():
        haplotectlocidf = pd.read_csv(haplotect_file,sep='\t')
        haplotectlocidf.columns = haplotectlocidf.columns.str.replace('#', '')
        df = haplotectlocidf.iloc[:-2]
    else:
        return None
    
    # parse haplotype strings
    df['sites'] = df['haplotype_counts'].str.strip(';')
    df['sites'] = df['sites'].str.split(';')
    df = df.explode('sites')
    df[['hap','counts']] = df['sites'].str.split(':',expand=True)
    df = df.dropna(subset=['counts'])
    df['counts'] = df['counts'].astype(int)
    # filter to get only homozygous and major heterozygous haplotypes (avoid contaminating haps)
    df['fraction'] = df['counts']/df['total_reads'].astype(int)
    df = df[df['fraction']>.40] # 40% is arbitrary but should be very specific
    # parse haplotypes into sites and make genotypes
    df[['allele1','allele2']] = df['hap'].apply(lambda x: pd.Series(list(x)))
    df1 = df[['chr','snp1','allele1']]
    df1.columns = ['chr','pos','gt']
    df2 = df[['chr','snp2','allele2']]
    df2.columns = ['chr','pos','gt']
    df = pd.concat([df1,df2],ignore_index=True).drop_duplicates()
    df = df.groupby(['chr','pos'])[['gt']].agg(list).reset_index()
    df['gt'] = df['gt'].apply(lambda x: "".join(sorted(x)) if len(x)==2 else "".join(x+x))
    df = df.sort_values(
        by=['chr','pos'],
        key=natsort_keygen())
    
    # pack into string chr1,pos1,gt1|chr2,pos2,gt2|...
    genotypes = '|'.join(df.apply(lambda row: f"{row['chr']},{row['pos']},{row['gt']}", axis=1))

    return genotypes

#
# Script
#

def main():

    parser = argparse.ArgumentParser(description='Collect Dragen QC and coverage metrics.')
    parser.add_argument('-c','--coverage-report',help='Mopath coverage report file')
    parser.add_argument('-m','--mapping-metrics',help='Dragen mapping metrics file')
    parser.add_argument('-n','--cnv-metrics',help='Dragen CNV metrics file')
    parser.add_argument('-u','--umi-metrics',help='Dragen UMI metrics file')
    parser.add_argument('-w','--wgs-coverage-metrics',help='Dragen WGS metrics file')
    parser.add_argument('-f','--wgs-fine-hist',help='Dragen WGS fine histogram file')
    parser.add_argument('-t','--haplotect',help='Haplotect file')

    parser.add_argument('-l','--coverage-qc-levels',default='10,20,40,60,80,100,150,200',help='Genome-wide coverage levels to collect for WGS assays.')
    parser.add_argument('-s','--coverage-summary-keys',default='gene,biomarker',help='Strings to use for target coverage summarization. [default: gene,biomarker]')

    # outfile name
    parser.add_argument('-o','--outfile',help='Output file name.',default=None)

    parser.add_argument('-v', '--version', action='version', version='%(prog)s: ' + __version__)

    args = parser.parse_args()

    # Exit and print usage if no arguments provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    qcDf = pd.DataFrame(columns=['category','metric','value'])
    covDf = pd.DataFrame(columns=['gene','region','start','end','region_type','mean_cvg','min_cvg','max_cvg'])

    if args.mapping_metrics is not None and Path(args.mapping_metrics).is_file():
        mapping_df = parse_dragen_metrics_file(args.mapping_metrics)
        if not qcDf.empty and not mapping_df.empty:
            qcDf = pd.concat([qcDf, mapping_df])
        elif not mapping_df.empty:
            qcDf = mapping_df

    if args.cnv_metrics is not None and Path(args.cnv_metrics).is_file():
        cnv_df = parse_dragen_metrics_file(args.cnv_metrics)
        #
        # Chromosome number and ploidy
        #
        ploidy = float(cnv_df[cnv_df["metric"].str.contains("Overall ploidy")]["value"].tolist()[0])
        dragen_chromosome_number = int(round(ploidy * 23, 0))
        chromosomes = dragen_chromosome_number if abs(dragen_chromosome_number - 46) > 2 else max(46, dragen_chromosome_number)
        cnv_df = pd.concat([cnv_df, pd.DataFrame([{'category': 'CNV SUMMARY', 'metric': 'CNV SUMMARY: Chromosome number', 'value': chromosomes}])], ignore_index=True)

        if not qcDf.empty and not cnv_df.empty:
            qcDf = pd.concat([qcDf, cnv_df])
        elif not cnv_df.empty:
            qcDf = cnv_df

    if args.umi_metrics is not None and Path(args.umi_metrics).is_file():
        umiDf = parse_dragen_metrics_file(args.umi_metrics)
        
        consensusReads = round(umiDf.loc[umiDf.metric=='UMI STATISTICS: Consensus pairs emitted','value'].astype(int).tolist()[0] * 2 / umiDf.loc[umiDf.metric=='UMI STATISTICS: Number of reads','value'].astype(int).tolist()[0] * 100,1)
        duplicateReads = round(100-consensusReads,1)

        umiDf = pd.concat([umiDf,pd.DataFrame.from_dict({0:['UMI STATISTICS','UMI STATISTICS: Consensus reads (%)',consensusReads,0],
                                                        1:['UMI STATISTICS','UMI STATISTICS: Duplicate reads (%)',duplicateReads,1]},orient='index',columns=['category','metric','value','qcmetric'])],axis=0)

        # Custom on-target read calculation
        # Enrichment Rate we use for UMI collapsed data is from the “On target number of reads ” from umi_metrics.csv divided by “Mapped reads” from mapping_metrics.csv.  As a first step to troubleshooting on-target rate, it would great to confirm if you are reporting this value.

        ontargetReads = umiDf.loc[umiDf.metric=='UMI STATISTICS: On target number of reads','value'].astype(int).tolist()[0]
        mappedReads = qcDf.loc[qcDf.metric=='MAPPING/ALIGNING SUMMARY: Mapped reads','value'].astype(int).tolist()[0]
        umiOnTargetRate = round(ontargetReads / mappedReads * 100,1)

        umiDf = pd.concat([umiDf,pd.DataFrame.from_dict({0:['UMI STATISTICS','UMI STATISTICS: On-target rate (%)',umiOnTargetRate,0]},orient='index',columns=['category','metric','value'])],axis=0)

        qcDf = pd.concat([qcDf,umiDf])

    if args.wgs_coverage_metrics is not None and Path(args.wgs_coverage_metrics).is_file():
        qcDf = pd.concat([qcDf,parse_dragen_metrics_file(args.wgs_coverage_metrics)])

    if args.wgs_fine_hist is not None and Path(args.wgs_fine_hist).is_file():
        qcDf = pd.concat([qcDf,parse_wgs_histogram(args.wgs_fine_hist, coverage_depths=[int(x) for x in args.coverage_qc_levels.split(',')])])

    if args.coverage_report is not None and Path(args.coverage_report).is_file():
        covDf = parse_coverage_report(args.coverage_report, aggregate_keys=args.coverage_summary_keys.split(',') if args.coverage_summary_keys else None)

        # Calculate coverage summary for assay targets
        coverageLevelLabels = covDf.columns[13:].tolist()
        covDf['bases'] = covDf['end'] - covDf['start'] + 1
        assayCov = round(covDf.apply(lambda x: x[coverageLevelLabels] / 100 * x['bases'],axis=1).sum() / covDf['bases'].sum() * 100,1)
        assayCov.index = assayCov.index.str.replace("pct_above_","COVERAGE SUMMARY: Percent of assay with coverage >") + 'x'
        assayCovDf = assayCov.to_frame(name='value').reset_index().rename(columns={0:'metric'})
        assayCovDf['category'] = 'COVERAGE SUMMARY'
        qcDf = pd.concat([qcDf,assayCovDf.reindex(columns=qcDf.columns)])

    if args.haplotect is not None and Path(args.haplotect).is_file():
        haplotectlocidf = pd.read_csv(args.haplotect,sep='\t')
        haplotectlocidf.columns = haplotectlocidf.columns.str.replace('#', '')
        haplotectdf = pd.DataFrame([haplotectlocidf.iloc[-1,:-2].tolist()],columns=haplotectlocidf.iloc[-2,:-2].tolist())
        haplotectdf.columns = haplotectdf.columns.str.replace('#', '')
        haplotectdf['informative_snppairs'] = haplotectdf['informative_snppairs'].astype(int)
        haplotectdf['mle_estimate'] = haplotectdf['mle_estimate'].astype(float)
        haplotectdf['contamination_fraction'] = haplotectdf['contamination_fraction'].fillna(0)
        haplotectdf['contamination_fraction'] = haplotectdf['contamination_fraction'].astype(float)

        haplotectlocidf = haplotectlocidf.iloc[:-2]
        # make distance, total_reads, haplotype_counts, contamination_fraction columns numeric
        cols_to_numeric = ['distance','total_reads','contamination_fraction']
        for col in cols_to_numeric:
            haplotectlocidf[col] = pd.to_numeric(haplotectlocidf[col], errors='coerce')

        haplotectdf = haplotectdf.transpose().reset_index().drop(index=0)
        haplotectdf.columns = ['metric','value']
        haplotectdf['category'] = 'HAPLOTECT'
        haplotectdf['metric'] = haplotectdf['metric'].apply(lambda v: f'HAPLOTECT: {v}')
        # add haplotect genotypes
        haplotectdf = pd.concat([
            haplotectdf,
            pd.DataFrame.from_dict(
                {0: ['HAPLOTECT', 'HAPLOTECT: Genotypes', pack_haplotect(haplotectlocidf=haplotectlocidf)]},
                orient='index',
                columns=['category', 'metric', 'value']
            )
        ], axis=0)

        qcDf = pd.concat([qcDf,haplotectdf.reindex(columns=qcDf.columns)])

    qcDf['qcmetric'] = 0
    qcDf['qcflag'] = 0

    ########################
    #
    # Start report
    #
    ########################

    if qcDf.empty and covDf.empty:
        print("No QC or coverage data found to report.",file=sys.stderr)
        sys.exit(0)

    # make dict for report and redirect output for text report
    jsonout = {}

    jsonout['ASSAY'] = qcDf.to_dict('split')
    jsonout['ASSAY'].pop('index', None)

    coverageLevelLabels = [x for x in covDf.columns if "pct_above" in x]

    xdf = covDf[(covDf['region_type']=='hotspot')][['gene','region','info','mean_cvg','min_cvg','max_cvg'] + coverageLevelLabels]
    jsonout['HOTSPOT_COVERAGE'] = xdf.to_dict('split')
    jsonout['HOTSPOT_COVERAGE'].pop('index', None)

    xdf = covDf[(covDf['region_type']=='gene') & (covDf['region']!='gene')][['gene','region','info','mean_cvg','min_cvg','max_cvg'] + coverageLevelLabels]
    jsonout['EXON_COVERAGE'] = xdf.to_dict('split')
    jsonout['EXON_COVERAGE'].pop('index', None)

    xdf = covDf[(covDf['region_type']=='gene') & (covDf['region']=='gene')][['gene','region','info','mean_cvg','min_cvg','max_cvg'] + coverageLevelLabels]
    jsonout['GENE_COVERAGE'] = xdf.to_dict('split')
    jsonout['GENE_COVERAGE'].pop('index', None)

    xdf = covDf[(covDf['region_type']=='sv')][['gene','region','info','mean_cvg','min_cvg','max_cvg'] + coverageLevelLabels]
    jsonout['TRANSCRIPT_COVERAGE'] = xdf.to_dict('split')
    jsonout['TRANSCRIPT_COVERAGE'].pop('index', None)

    # biomarker coverage
    xdf = covDf[(covDf['region_type']=='biomarker')][['gene','region','info','mean_cvg','min_cvg','max_cvg'] + coverageLevelLabels]
    jsonout['BIOMARKER_COVERAGE'] = xdf.to_dict('split')
    jsonout['BIOMARKER_COVERAGE'].pop('index', None)

    jsonout['HAPLOTECT_LOCI'] = haplotectlocidf.to_dict('split')
    jsonout['HAPLOTECT_LOCI'].pop('index', None)

    # dump json to outfile or stdout
    if args.outfile is None:
        json.dump(sanitize_for_json(jsonout),sys.stdout,indent=" ",allow_nan=False)
    
    else:
        j = open(args.outfile, "w")
        json.dump(sanitize_for_json(jsonout),j,indent=" ",allow_nan=False)
        j.close()


if __name__ == "__main__":
    main()
