#!/usr/bin/perl

use strict;
use JSON;
use Getopt::Long;

my $usage = "Usage: $0 -p|project <project-name> -o|out <outdir> (required)\nAny one of:\n\t<LIMS dir 1> [ <LIMS dir 2> ...] (positional LIMS directories)\n\t-b|biosample <biosamplename> (launch entire process with uploaded data)\n\t-d|dragen <dragen session id> (launch chromoseq on aligned data)\n\t-c|chromoseq <chromoseq session id> (download chromoseq session data)\n";

my $inputs = '/gscmnt/gc2555/spencer/dhs/git/docker-basespace_chromoseq/inputs.v17.local.json';
my $wdl = '/gscmnt/gc2555/spencer/dhs/git/docker-basespace_chromoseq/Chromoseq.v17.wdl';
my $conf = '/gscmnt/gc2555/spencer/dhs/projects/wdltest/application.new.conf';
my $fa = '/gscmnt/gc2555/spencer/refdata/hg38/all_sequences.fa';

my $group = '/dspencer/chromoseq';

my $dir = '';
my $name = '';
my $out = '';

my $testing = 0;
my $debug = 0;
my $force = 0;

my $gender = '';

GetOptions("debug" => \$debug,
	   "gender=s" => \$gender,
	   "i|inputs=s" => \$inputs,
	   "w|wdl=s" => \$wdl,
	   "d|dir=s" => \$dir,
	   "n|name=s" => \$name,
	   "o|out=s" => \$out,
	   "g|group=s" => \$group,
	   "t" => \$testing,
	   "f" => \$force);

die "$0 -d dir -n name -o out" if !-d $dir and !$name and !-d $out;

$wdl = `readlink -f $wdl`;
chomp $wdl;

my $inputs = from_json(`cat $inputs`);

if ($gender =~ /^male$|^female$/){
    $inputs->{'ChromoSeq.Gender'} = $gender;

} elsif (`grep -c estimation $dir/*.ploidy_estimation_metrics.csv` =~ /1/) {
    `grep estimation $dir/*.ploidy_estimation_metrics.csv | cut -d ',' -f 4` =~ /(XX|XY)/;
    if ($1 eq 'XX'){
        $inputs->{'ChromoSeq.Gender'} = "female";
    } elsif ($1 eq 'XY') {
        $inputs->{'ChromoSeq.Gender'} = "male";
    }
} elsif (`grep -c XX $dir/*.wgs_ploidy.csv` =~ /1/){
    $inputs->{'ChromoSeq.Gender'} = "female";

} elsif (`grep -c XY $dir/*.wgs_ploidy.csv` =~ /1/){
    $inputs->{'ChromoSeq.Gender'} = "male";
    
} elsif (`grep -c XAvgCov/YAvgCov $dir/*.wgs_coverage_metrics.csv` =~ /1/) {
    `grep XAvgCov/YAvgCov $dir/*.wgs_coverage_metrics.csv | cut -d ',' -f 4` =~ /([0-9.]+)/;
    if ($1 > 1.75){
	$inputs->{'ChromoSeq.Gender'} = "female";
    } elsif ($1 < 1.75) {
	$inputs->{'ChromoSeq.Gender'} = "male";
    }
} elsif (`grep -c XAvgCov/YAvgCov $dir/*.mapping_metrics.csv` =~ /1/) {
    `grep XAvgCov/YAvgCov $dir/*.mapping_matrics.csv | cut -d ',' -f 4` =~ /([0-9.]+)/;
    if ($1 > 1.75){
	$inputs->{'ChromoSeq.Gender'} = "female";
    } elsif ($1 < 1.75) {
	$inputs->{'ChromoSeq.Gender'} = "male";
    }
}

die "Cant find gender\n" if $inputs->{'ChromoSeq.Gender'} eq '';

my $mapsum = `readlink -f $dir/*.mapping_metrics.csv`;
chomp $mapsum;
die if ! -e $mapsum;
$inputs->{'ChromoSeq.MappingSummary'} = $mapsum;

my $covsum = `readlink -f $dir/*.wgs_coverage_metrics.csv `;
chomp $covsum;
$inputs->{'ChromoSeq.CoverageSummary'} = $covsum if -e $covsum;

my $cram = `readlink -f $dir/*.cram`;
chomp $cram;
die if ! -e $cram;
$inputs->{'ChromoSeq.Cram'} = $cram;

my $crai = `readlink -f $cram.crai`;
chomp $crai;
die if ! -e $crai;
$inputs->{'ChromoSeq.CramIndex'} = $crai;

my $counts = `readlink -f $dir/*.target.counts`;
chomp $counts;
die "$name" if ! -e $counts;
$inputs->{'ChromoSeq.TumorCounts'} = $counts;

$inputs->{'ChromoSeq.Name'} = $name;

$out = `readlink -f $out`;
chomp $out;
mkdir $out if ! -e $out;

$inputs->{'ChromoSeq.OutputDir'} = $out;

open(JS,">$out/chromoseq.$name.json") || die;
print JS to_json($inputs);
close JS;

my $cmd = "bsub -g $group -oo $out/out.$name.log -eo $out/err.$name.log -q research-hpc -a \"docker(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:7)\" /usr/bin/java -Dconfig.file=$conf -jar /opt/cromwell.jar run -t wdl -i $out/chromoseq.$name.json $wdl";

print STDERR $cmd,"\n";

`$cmd` if !$testing;
