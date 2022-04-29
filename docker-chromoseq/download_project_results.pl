#!/usr/bin/perl

use strict;
use JSON;
use YAML::Tiny;

my $BS = '/usr/local/bin/bs';

my $DRAGENAPP = 6495489;
my $CHROMOSEQAPP = 6984978;

my $RefId = 14321367413;

my $debug = "";

my $ProjectName = shift @ARGV;

mkdir "$ProjectName";

my $batch_json = from_json(`$BS appsession list --project-name=$ProjectName --filter-term=DRAGEN --exec-status=Complete -f json`);

foreach my $case (@{$batch_json}){
  my $dataset_json = from_json(`bs appsession property get -i $case->{Id} --property-name=Output.Datasets -f json`);
  `$BS download dataset -i $dataset_json->{Id} --extension=csv -o \"$ProjectName/$dataset_json->{Name}\"`;
}

#my $batch_json = from_json(`$BS appsession list --project-name=$ProjectName --filter-term=Chromoseq --exec-status=Complete -f json`);

#foreach my $case (@{$batch_json}){
#  my $dataset_json = from_json(`bs appsession property get -i $case->{Id} --property-name=Output.Datasets -f json`);
#  `$BS download dataset -i $dataset_json->{Id} -o \"$ProjectName/$dataset_json->{Name}\"`;
#}
