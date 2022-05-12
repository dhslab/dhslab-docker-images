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
#my $ProjectId = 126915789; # a dummy project

my $outdir = shift @ARGV;

mkdir $outdir if !-e $outdir;

my @dirs = @ARGV;

map { die "$_ doesnt exist" if ! -e $_ } @dirs;

# go through datasets and make sure they're all from the same library
my %datasets = ();

foreach my $dir (@dirs){ 
    # parse manifest.
    my $yaml = `ls $dir/*.yaml`;
    chomp $yaml;
    my $y = new YAML::Tiny;
    my %dat = %{($y->read($yaml))->[0]};

    my $name = $dat{library_summary}{full_name};
    $datasets{$name} = \%dat;
}

die "Multiple directories with different sample/library information detected!" if scalar keys %datasets > 1;

# get dataset info, make biosample, and upload

my %manifest = %{$datasets{(keys %datasets)[0]}};
my $biosamplename = $debug . $manifest{sample}{full_name};

# first check to see if biosample is present and create it if necessary
my $biosample_json = '';
my @biosample_check = `$BS list biosample --project-name $ProjectName -f csv`;
chomp @biosample_check;
my %biosamples = map { my @l = split(",",$_); $l[0] => $l[1] } @biosample_check[1..$#biosample_check];

if (!defined($biosamples{$biosamplename})){
    
    print STDERR "No samples called $biosamplename found. Creating it...\n";
    
    my $metadata = join(" ",(map { "--metadata Sample.$_:$manifest{sample}{$_}" } keys %{$manifest{sample}}),
			(map { "--metadata LibrarySummary.$_:$manifest{library_summary}{$_}" } keys %{$manifest{library_summary}}));
    
    `$BS create biosample -n $biosamplename $manifest{library}{full_name} -p $ProjectName $metadata | tee $biosamplename.$manifest{index_illumina}{analysis_id}.create_biosample.json`;
    
}

sleep 10;

my $biosample_json = from_json(`$BS biosample get -n $biosamplename -f json`);

# get project id
my $ProjectId = `bs project list --filter-term \"^$ProjectName\$\" --terse`;
chomp $ProjectId;

my $dirs = join(" ", @dirs);

my $label = "Upload $biosamplename " . localtime();
print STDERR "Uploading $biosamplename...\n";
my $upload_json = `$BS upload dataset -p $ProjectId --recursive --biosample-name=$biosamplename --library-name=$manifest{library_summary}{full_name} -l \"$label\" $dirs`;

sleep 30;

$label = "Dragen $biosamplename " . localtime();
my $align_json = from_json(`$BS launch application -i $DRAGENAPP -o app-session-name:\"$label\" -o project-id:$ProjectId -o ht-ref:custom.v7 -o ht-id:$RefId -o input_list.tumor-sample:$biosample_json->{Id} -o dupmark_checkbox:1 -o bai_checkbox:1 -o output_format:CRAM -f json | tee $biosamplename.$manifest{index_illumina}{analysis_id}.dragen.json`);

print STDERR "Launched Dragen: $label. Waiting...\n";

# now wait
#my $align_result = from_json(`$BS await appsession $align_json->{Id} -f json`);

#die "Dragen failed! Check logs for AppSession: $align_result->{AppSession}{Name}" if ($align_result->{AppSession}{ExecutionStatus} !~ /Complete/);

#print STDERR "Dragen finished.\n";

# get cram files
#my @cram = `$BS contents dataset -i $align_result->{Id} -f csv | grep cram | cut -d ',' -f 1`;
#chomp @cram;
#my $files = join(",",@cram);

# launch chromoseq
#$label = "Chromoseq $biosamplename " . localtime();
#my $chromoseq_session = from_json(`$BS launch application -i $CHROMOSEQAPP -o app-session-name:\"$label\" -o project-id:$ProjectId -o file-id:$files -f json`);

#$chromoseq_result = from_json(`$BS await appsession $chromoseq_session->{Id} -f json`);

#print STDERR "Launched Chromoseq: $label. Waiting...\n";

#die "Chromoseq failed! Check logs for AppSession: $chromoseq_result->{AppSession}{Name}" if ($chromoseq_result->{AppSession}{ExecutionStatus} !~ /Complete/);

#print STDERR "Chromoseq finished.\nDownloading files to $outdir.\n";

# download files

#`$BS download dataset -i $align_json->{Id} --extension=csv -o \"$ProjectName/$dataset_json->{Name}\"`;

#`$BS download dataset -i $chromoseq_result->{Id} -o $outdir`;

#print STDERR "Done.\n";

