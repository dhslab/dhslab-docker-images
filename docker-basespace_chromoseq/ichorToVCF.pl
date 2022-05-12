#!/usr/bin/perl

use File::Basename;
use Getopt::Long;
use strict;

my $refseq = "~/refdata/hg38/all_sequences.fa";

GetOptions("r=s" => \$refseq);

my $segsfile = $ARGV[0];

my $name = basename($segsfile,".seg.txt");

print <<EOF;
##fileformat=VCFv4.2
##source=ichorCNA
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=LOG2RATIO,Number=.,Type=Float,Description="Log2 ratio from CNV analysis">
##INFO=<ID=CNBINS,Number=1,Type=Integer,Description="Number of CN bins in CNA call">
##INFO=<ID=POS,Number=1,Type=Integer,Description="Position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=set,Number=1,Type=String,Description="Source VCF for the merged record in CombineVariants">
##INFO=<ID=CN,Number=.,Type=Integer,Description="Integer copy number estimate from CNV analysis">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$name
EOF

# now open CNV file and print as VCF records
open(CNV,"$segsfile") || die "cant open CNV file";
<CNV>; # remove header
while(<CNV>){
    chomp;
    next if /NEUT/;

    my @F = split("\t",$_);
    my $svtype = "DEL";
    $svtype = "DUP" if $F[5] > 0;
    my $pos = $F[2]-1;
    my $refnt = `samtools faidx $refseq $F[1]:$pos-$pos | tail -n 1`;
    chomp $refnt;
    my $svlen = ($F[3]-$F[2]+1);
    $svlen = -$svlen if $svtype eq 'DEL';
    print join("\t",$F[1],$F[2]-1,"ICHOR:$F[1]_$F[2]_$F[3]",$refnt,"<$svtype>",".","PASS",
		 join(";","SVTYPE=$svtype","LOG2RATIO=$F[5]","CN=$F[6]","CNBINS=$F[4]","POS=".($F[2]-1),"END=$F[3]","SVLEN=$svlen","IMPRECISE"),
		 "GT","./."),"\n";
}
