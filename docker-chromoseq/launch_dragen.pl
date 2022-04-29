#!/usr/bin/perl

use strict;
use JSON;
use YAML::Tiny;

my $BS = '/usr/local/bin/bs';

my $DRAGENAPP = 6495489;
my $CHROMOSEQAPP = 6984978;

my $RefId = 14321367413;

my $debug = "TESTING6-";

my $ProjectName = "TestProject";
#my $ProjectId = 126915789;

# launch dragen and wait

$align_appsession = `$BS launch application -i $DRAGENAPP -o app-session-name:"$ts" -o project-id:$ProjectId -o ht-ref:custom.v7 -o ht-id:14321367413 -o input_list.tumor-sample:$Name -o dupmark_checkbox:1 -o bai_checkbox:1 -o output_format:CRAM

$align_output = ./bs await appsession $align_appsession

# get completion of alignment and check to make sure its complete
./bs appsession get -i $align_appsession

# get cramfiles
./bs contents dataset -i $align_appsession{Id} -f tsv

# parse to get $cram_id
$chromoseq_appsession = ./bs launch application -i $CHROMOSEQAPP -o app-session-name:"$run_chromoseq_name" -o project-id:$Project -o file-id:$cram_id

$chromoseq_output = ./bs await appsession $chromoseq_appsession

# get status and if is good, download files
./bs appsession get -i $chromoseq_appsession

# if status complete
./bs download dataset -i $chromoseq_output{'Id'} -o ./
