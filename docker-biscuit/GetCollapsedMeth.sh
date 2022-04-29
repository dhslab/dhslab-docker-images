#!/bin/bash

VCF=$1
BED=$2

dir=$(dirname $VCF)
base=$(basename $VCF .meth.vcf.gz)

set -euf -o pipefail && /usr/bin/biscuit vcf2bed -e -k 1 $VCF | ( [[ "${BED}" ]] && intersectBed -a stdin -b $BED || tee /dev/null ) | \
	awk -v OFS="\t" '{ if ($4=="C") { $3=$3+1; } else { $2=$2-1; } print $1,$2,$3,sprintf("%d",$8*$9),$9; }' | \
	bedtools groupby -g 1,2,3 -c 4,5 -o sum,sum | awk -v OFS="\t" '{ print $1,$2,$3,sprintf("%.3f",$4/$5),$5; }' | \
	sortBed -i stdin | tee $dir/$base".meth.bed" | cut -f 1-4 | bgzip -c > $dir/$base".meth.bg.gz" && \
    bgzip -c $dir/$base".meth.bed" > $dir/$base".meth.bed.gz" && rm -f $dir/$base".meth.bed" && \
    Rscript --default-packages=data.table,bsseq -e 'args <- commandArgs(TRUE); stopifnot(all(length(args)==2,file.exists(args[1]))); x <- data.table(read.table(args[1],colClasses=c("character","integer","integer","numeric","integer"),col.names=c("chr","start","end","ratio","coverage"))); x$numCs <- round(x$ratio * x$coverage,0); bs <- BSseq(M = as.matrix(x$numCs),Cov = as.matrix(x$coverage),gr=GRanges(x$chr,IRanges(x$start+1,x$start+1)),sampleNames=args[2]); save(bs,file=paste0(dirname(args[1]),"/",args[2],".bs.rda"))' $dir/$base".meth.bed.gz" $base
