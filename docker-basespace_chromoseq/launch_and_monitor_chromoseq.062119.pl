#!/usr/bin/perl

use strict;
use JSON;
use YAML::Tiny;
use Getopt::Long;

my $usage = "Usage: $0 -p|project <project-name> -o|out <outdir> (required)\nAny one of:\n\t<LIMS dir 1> [ <LIMS dir 2> ...] (positional LIMS directories)\n\t-b|biosample <biosamplename> (launch entire process with uploaded data)\n\t-d|dragen <dragen session id> (launch chromoseq on aligned data)\n\t-c|chromoseq <chromoseq session id> (download chromoseq session data)\n";

my $BS = '/usr/local/bin/bs';

my $DRAGENSOMATICAPP = 6495489;
my $DRAGENGERMLINEAPP = 6840834;
my $CHROMOSEQAPP = 6984978;

my $RefId = 14321367413;
my $ref_fasta = 14477433053;

my $debug = "";

my $ProjectName = "TestProject";
my $outdir = "./";

my $dragensession = '';
my $biosamplename = '';
my $chromoseqsession = '';

GetOptions("debug" => \$debug,
	   "b|biosample=s" => \$biosamplename,
	   "p|project=s" => \$ProjectName,
	   "d|dragen=s" => \$dragensession,
	   "c|chromoseq=s" => \$chromoseqsession,
	   "o|out=s" => \$outdir);

die $usage if (!$ProjectName || !$outdir || (scalar @ARGV < 1 && !$biosamplename && !$dragensession && !$chromoseqsession)) || ($biosamplename && $dragensession) || ($biosamplename && $chromoseqsession) || ($chromoseqsession && $dragensession);

# get project id
my $ProjectId = `bs project list --filter-term \"^$ProjectName\$\" --terse`;
chomp $ProjectId;

die "failed to get project id for $ProjectName" if !$ProjectId;

mkdir $outdir if !-e $outdir;
$outdir = `readlink -f $outdir`;
chomp $outdir;

# passed positional arguments. assume upload and launch workflow
if (scalar @ARGV > 0){
  my @dirs = @ARGV;

  map { die "$_ doesnt exist" if ! -e $_ } @dirs;
  
  # go through datasets and make sure they're all from the same library
  my %datasets = ();

  foreach my $dir (@dirs){ 
    # parse manifest.
    my $yaml = `ls $dir/*.yaml`;
    chomp $yaml;
    die "No metadata file found in $dir" if !$yaml;
    
    my $y = new YAML::Tiny;
    my %dat = %{($y->read($yaml))->[0]};
    
    my $name = $dat{library_summary}{full_name};
    $datasets{$name} = \%dat;
  }
  
  die "Multiple directories with different sample/library information detected!" if scalar keys %datasets > 1;
  
  ## get dataset info, make biosample, and upload
  my %manifest = %{$datasets{(keys %datasets)[0]}};
  $biosamplename = $manifest{library_summary}{full_name};

  die "Sample name not found" if !$biosamplename;
  
  $biosamplename = $biosamplename . "-DEBUG" if $debug;
  
  $outdir = $outdir . '/' . $biosamplename;
  mkdir $outdir if !-e $outdir;
  
  # first check to see if biosample is present and create it if necessary
  my $biosample_json = '';
  my @biosample_check = `$BS list biosample -f csv`;
  chomp @biosample_check;
  my %biosamples = map { my @l = split(",",$_); $l[0] => $l[1] } @biosample_check[1..$#biosample_check];
  
  die "Biosample $biosamplename is already on basespace. Exiting" if (defined($biosamples{$biosamplename}));
  
  print STDERR "No samples called $biosamplename found. Creating it...\n";
  
  my $metadata = join(" ",(map { "--metadata Sample.$_:$manifest{sample}{$_}" } keys %{$manifest{sample}}),
		      (map { "--metadata LibrarySummary.$_:$manifest{library_summary}{$_}" } keys %{$manifest{library_summary}}));
  
  `$BS create biosample -n $biosamplename $manifest{library}{full_name} -p $ProjectName $metadata | tee $outdir/$biosamplename.$manifest{index_illumina}{analysis_id}.create_biosample.json`;
  
  # for some reason it takes a few seconds for the new biosample to register
  sleep 10;
  
  my $dirs = join(" ", @dirs);
  
  my $label = '';
  
  my $label = "Upload $biosamplename " . localtime();
  print STDERR "Uploading $biosamplename...\n";
  `$BS upload dataset -p $ProjectId --recursive --biosample-name=$biosamplename --library-name=$manifest{library_summary}{full_name} -l \"$label\" $dirs`;
  
  # added another sleep just to be safe
  sleep 10;
}

# if just uploaded data to biosample or passed one, run dragen
if ($biosamplename){
  my $biosample_json = from_json(`$BS biosample get -n $biosamplename -f json`);
  
  my $label = "Dragen $biosamplename " . localtime();
  my $align_json = from_json(`$BS launch application -i $DRAGENGERMLINEAPP -o app-session-name:\"$label\" -o project-id:$ProjectId -o pipeline-mode:0 -o ht-ref:custom.v7 -o ht-id:$RefId -o input_list.sample-id:$biosample_json->{Id} -o dupmark_checkbox:1 -o bai_checkbox:1 -o output_format:CRAM -f json | tee $outdir/$biosamplename.dragen.json`);
  
  print STDERR "Launched Dragen: $label. Waiting...\n";

  # now wait
  my $align_result = from_json(`$BS await appsession $align_json->{Id} -f json`);
  
  die "Dragen failed! Check logs for AppSession: $align_result->{AppSession}{Name}" if ($align_result->{AppSession}{ExecutionStatus} !~ /Complete/);

  print STDERR "Dragen finished. Downloading data\n";
  
  `$BS dataset download -i $align_result->{Id} -o $outdir`;

  $dragensession = $align_result->{AppSession}{Id};
}

# if passed a dragen session id or just did alignment, launch chromoseq
if ($dragensession){
  my $dragenresult = `$BS appsession property get -i $dragensession --property-name=Output.datasets --terse`;
  chomp $dragenresult;

  die "failed to get dataset for dragen session $dragensession" if !$dragenresult;

  # get cram files, first by getting the appresult for the dragen process, then the file id from that process
  my $appresultId = `$BS dataset get -i $dragenresult -F V1pre3Id -f csv`;
  chomp $appresultId;
  $appresultId =~ s/[^0-9]+//g;

  die "failed to get appresult for dragen session $dragensession" if !$appresultId;
  
  my @cram = `$BS appresult content -i $appresultId --extension=cram,crai --terse`;
  chomp @cram;
  
  die "Error getting cram and cram index dragen session $dragensession" if scalar @cram != 2;
  
  my $files = join(",",@cram);
  
  # launch chromoseq
  my $label = "Chromoseq $biosamplename " . localtime();
  my $chromoseq_session = from_json(`$BS launch application -i $CHROMOSEQAPP -o app-session-name:\"$label\" -o project-id:$ProjectId -o ref-fa-id:$ref_fasta -o file-id:$files -f json | tee $outdir/$biosamplename.chromoseq.json`);

  print STDERR "Launching Chromoseq: $label. Waiting...\n";
  
  my $chromoseq_result = from_json(`$BS await appsession $chromoseq_session->{Id} -f json`);
  
  die "Chromoseq failed! Check logs for AppSession: $chromoseq_result->{AppSession}{Name}" if ($chromoseq_result->{AppSession}{ExecutionStatus} !~ /Complete/);
  
  print STDERR "Chromoseq finished.\n";
  
  $chromoseqsession = $chromoseq_result->{AppSession}{Id};

}

# if passed chromoseq session, download files
if ($chromoseqsession){
  
  my $chromoseqresult = `$BS appsession property get -i $chromoseqsession --property-name=Output.datasets --terse`;
  chomp $chromoseqresult;

  die "failed to get dataset for chromoseq session $chromoseqsession" if !$chromoseqresult;

  # download files
  print STDERR "Downloading Chromoseq results.\n";
  `$BS download dataset -i $chromoseqresult -o $outdir`;

}

print STDERR "Done.\n";
  
