#!/usr/bin/perl

use strict;
use JSON qw(from_json to_json);
use YAML::Tiny;
use Getopt::Long;
use File::Basename qw(dirname basename);
use File::Spec;
use Cwd qw(abs_path realpath);

BEGIN {
    package CS;

    use File::Basename qw(dirname basename);
    use JSON qw(from_json to_json);
    
    sub new {
	my ($caller, %in) = @_;
	my $class = ref($caller) || $caller;
	
	my $usage = "$0 <command> [options]\n\nCommands:\n\tprocess -d|dir <flowcell dir or fastq dir\n\tdemux -d|dir <flowcell directory> -o|out <outdir>\n\tmap -o|out <outdir> -d|dir <fastq dir> or -f <fastq samplesheet> [-n|name <sample>\n\tanalyze -c|cram <cram> -n|name <sample name> -o|out <outdir>";	 
	
	my $obj = { _inputs => '/gscmnt/gc2555/spencer/dhs/git/docker-basespace_chromoseq/inputs.v17.local.json',
		    _wdl => '/gscmnt/gc2555/spencer/dhs/git/docker-basespace_chromoseq/Chromoseq.v17.wdl',
		    _conf => '/gscmnt/gc2555/spencer/dhs/projects/wdltest/application.new.conf',
		    _ref => '/gscmnt/gc2555/spencer/refdata/hg38/all_sequences.fa',
		    _dragenref => "/staging/garza_testing/reference/",
		    _fastqdir => "/staging/cle/fastqs",
		    _cramout => "/staging/cle/outputs",
		    _instrument_data => "/gscmnt/gc2648/wgs/dragen_testing/instrument_data/",
		    _drtmp => "/mnt/gc2648/wgs/dragen_user/tmp",
		    _drtmponlocal => "/gscmnt/gc2648/wgs/dragen_user/tmp",
		    _tmp => "/gscmnt/gc2648/wgs/dragen_testing/tmp",
		    _tmpondr => "/mnt/gc2648/wgs/dragen_testing/tmp",
		    _fastqstotransfer => [],
		    _toclean => [],
		    _run => join("",map { sprintf q|%X|, rand(16) } 1 .. 10),
		    _outdir => '',
		    _samples => {},
		    _genders => [],
		    _group => '/dspencer/chromoseq',
		    _force => 0,
		    _testing => 0
	};
    
	foreach my $k (keys %in){
	    $obj->{"_" . $k} = $in{$k} if $in{$k};
	}
	
	bless $obj, $class;
	$obj;
    }
    
    sub demux {
	my $obj = shift;
	my @dirs = @_;
	
	# make final fastqfile file on dragen where fastqs will be listed for subsequent mapping
	my $drfastqfile = $obj->{_drtmp} . "/cs.fastqsheet." . $obj->{_run} . ".csv";
	my $rc = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "echo RGID,RGSM,RGLB,Lane,Read1File,Read2File > $drfastqfile" || echo 1`;
	die "Failed making fastq sheet $drfastqfile\n" if $rc =~ /1/;
	
	# query run folder on dragen until sequencing is complete and then demux and return the fastqfile for mapping
	my $done = 0;
	while($done < scalar @dirs){
	    foreach my $d (@dirs){
		
		$d =~ s/gscmnt/mnt/; s/production\///;
		
		my $rc = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "test -f $d/*/SequenceComplete.txt && echo 1"`;
		
		if ($rc =~ /1/){
		    my $dir = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "ls $d/*/SequenceComplete.txt"`;
		    chomp $dir;
		    $dir = dirname($dir);
		    
		    my $fcname = basename($dir);
		    
		    # find samplesheet
		    my $ss = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "test -f $d/*/*.csv || echo 1"`;
		    die "Cant find samplesheet for $dir\n" if $ss =~ /1/;
		    
		    $ss = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "ls $d/*/*.csv"`;
		    chomp $ss;
		    
		    my $fastqdir = $obj->{_fastqdir} . "/$fcname";
		    my $fclog = $obj->{_fastqdir} . "/$fcname.log";
		    
		    # see if demux has already been run                                                                                       
		    my $rv = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "test -f $fclog && test -d $fastqdir && echo 1"`;
		    if ($rv =~ /1/){
			$rv = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "grep -c normally $fclog"`;
			if ($rv =~ /1/){
			    print STDERR "Demux finished normally in $fastqdir. Will proceed.\n\n";
			    
			    # add to fastqfile
			    
			    $rc = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "tail -n +2 $fastqdir/Reports/fastq_list.csv >> $drfastqfile" || echo 1`;
			    die "Cant add $fcname to fastq sample sheet\n" if $rc =~ /1/;
			    
			    push @{$obj->{_fastqstotransfer}}, $fastqdir;
			    
			    $done++;
			    next;
			    
			} elsif ($obj->{_force}){
			    $rv = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "rm -Rf $fastqdir || echo 1"`;
			    die "Failed removing $fastqdir prior to demux" if $rv =~ /1/;
			}
		    }
		    
		    print STDERR "Starting demux with the dragen for $d\n\n";
		    
		    my $cmd = "dragen --bcl-conversion-only true --sample-sheet $ss --bcl-input-directory $dir --output-directory $fastqdir &> $fclog";
		    $cmd = "dragen --first-tile-only true --bcl-conversion-only true --sample-sheet $ss --bcl-input-directory $dir --output-directory $fastqdir &> $fclog" if $obj->{_testing};
		    
		    print STDERR "$cmd\n\n";
		    
		    $rc = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "$cmd || echo 1"`;
		    die "Demux for $d didnt work!\n" if $rc =~ /1/;
		    
		    # add to master fastqfile
		    $rc = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "tail -n +2 $fastqdir/Reports/fastq_list.csv >> $drfastqfile" || echo 1`;
		    die "Cant add $fcname to fastq sample sheet\n" if $rc =~ /1/;
		    
		    push @{$obj->{_fastqstotransfer}}, $fastqdir;
		    
		    $done++;
		    
		} else {
		    print STDERR "\nStill waiting...\n";
		    sleep 300;
		    
		}
	    }
	}
	return $drfastqfile;
    }

    sub process_instrument_data {
	my $obj = shift;
	my @dirs = map { Cwd::realpath($_) } @_;

	my %samples = ();
	foreach my $i (0..$#dirs){
	    my $d = $dirs[$i];

	    my $rc = '';
	    
	    my $yaml = `ls $d/*.yaml`;
	    chomp $yaml;

	    if ($yaml and -e $yaml){
		
		my $y = new YAML::Tiny;
		my %dat = %{($y->read($yaml))->[0]};
		
		if (-e "$obj->{_outdir}/$dat{read_group_header}{LB}" and $rc = `grep -c normally $obj->{_outdir}/$dat{read_group_header}{LB}/run_dragen.log` and $rc =~ /1/ and !$obj->{_force}){
		    print STDERR "Looks like valid Dragen output exists for $dat{read_group_header}{LB}, not including this sample.\n\n";
		    next;
		    
		} elsif (-e "$obj->{_outdir}/$dat{read_group_header}{LB}" and $obj->{_force}){
		    `rm -Rf "$obj->{_outdir}/$dat{read_group_header}{LB}"`;
		}
		$dat{data_dir} = $d;
		$samples{$dat{read_group_header}{LB}} = 1;
		push @{$obj->{_samples}{$dat{read_group_header}{LB}}}, \%dat;
	    }
	}
	return (keys %samples);
    }
    
    sub prepare_fastqs {
	my $obj = shift;
	my $s = shift;       
	
	# Make the fastq dir and samplesheet containing all fastqs/samples to process
	my $fqdir = '';
	my $fastqfile = $obj->{_tmp} . "/" . $obj->{_run} . "." . $s . ".fastq_list.csv";

	if (-e $fastqfile){
	    unlink($fastqfile) or die "fastq samplesheet exists and cant remove!";
	}
	my $rc = `echo RGID,RGSM,RGLB,Lane,Read1File,Read2File > $fastqfile || echo 1`;
	die "Failed making fastq sheet $fastqfile\n" if $rc =~ /1/;
	
	# make dragen-accessible version
	my $drfastqfile = $fastqfile;
	$drfastqfile =~ s/gscmnt/mnt/;
	
	foreach my $i (0..$#{$obj->{_samples}{$s}}){	    
	    my %dat = %{$obj->{_samples}{$s}[$i]};
	    my $d = $dat{data_dir};
	    
	    my $fqdirname = basename($d);
	    
	    print STDERR "copying fastqs from $d to dragen...";
	    
	    $rc = `scp -r $d gtac\@compute1-dragen-1.ris.wustl.edu:$obj->{_fastqdir}/ || echo 1`;
	    die "Failed to transfer fastqs\n" if $rc =~ /1/;	     
	    
	    my $R1 = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "readlink -f $obj->{_fastqdir}/$fqdirname/*R1*.fastq.gz"`;
	    chomp $R1;
	    my $R2 = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "readlink -f $obj->{_fastqdir}/$fqdirname/*R2*.fastq.gz"`;
	    chomp $R2;
	    
	    my $s = join(',',$dat{read_group_header}{ID},$dat{read_group_header}{SM},$dat{read_group_header}{LB},$dat{index_illumina}{lane},$R1,$R2);
	    `echo $s >> $fastqfile`;

	    push @{$obj->{_toclean}}, "$obj->{_fastqdir}/$fqdirname";
	    
	    print STDERR "done.\n\n";
	    
	}
	return $drfastqfile;
    }
    
    sub align {
	my $obj = shift;
	my $fastqfile = shift @_;	

	my @cramdirs = ();
	
	# chech to be sure the fastqfile is there
	my $rc = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "test -e $fastqfile || echo 1"`;
	die "fastqfile $fastqfile does not exist" if $rc =~ /1/;

	my %samples = ();
	my @ss = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "cat $fastqfile"`;
	chomp @ss;
	shift @ss;
	foreach my $l (@ss){
	    my @x = split(',',$l);
	    $samples{$x[0]} = $x[2];
	}

	foreach my $samplename (keys %samples){
	    my $sample = $samples{$samplename};
	    my $samplename = $samples{$samplename};

	    my $out = $obj->{_outdir};
	    
	    # check to see if the cram has been generated already
	    if (-e "$out/$sample" and $rc = `grep -c normally $out/$sample/run_dragen.log` and $rc =~ /1/){
		print STDERR "Looks like valid Dragen output exists for $sample here: $out, continuing.\n\n";
		
	    } else {
		
		my $cramout = $obj->{_cramout} . "/$sample";
		my $cramtmp = $obj->{_drtmp} . "/$sample";
		
		# check to see if the dragen ran normally on the dragen
		if ($rc = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "test -d $cramout && grep -c normally $cramout/run_dragen.log"` and $rc =~ /1/){
		    print STDERR "Looks like valid Dragen output exists for $sample here: $cramout, continuing.\n\n";
		    
		    # what about in the temp dragen location? 
		} elsif ($rc = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "test -d $cramtmp && grep -c normally $cramtmp/run_dragen.log"` and $rc =~ /1/){
		    print STDERR "Looks like valid Dragen output exists for $sample here: $cramtmp, continuing.\n\n";
		    
		} else {
		    if ($rc = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "test -d $cramout && echo 1"` and $rc =~ /1/ and $obj->{_force}){
			
			$rc = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "rm -Rf $cramout || echo 1"`;
			die "Cant make Output folder $cramout" if $rc =~ /1/;
			$rc = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "mkdir $cramout || echo 1"`;
			die "Cant make Output folder $cramout" if $rc =~ /1/;
			
		    } else {
			$rc = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "mkdir $cramout || echo 1"`;
			die "Cant make Output folder $cramout" if $rc =~ /1/;
		    }
		    
		    # map reads with the dragen
		    my $cmd = "dragen -r $obj->{_dragenref} --fastq-list $fastqfile --fastq-list-sample-id $samplename --enable-cnv true --cnv-target-bed /staging/garza_testing/reference/all_sequences.fa.bed --cnv-interval-width 500000 --output-directory $cramout --output-file-prefix $sample --output-format CRAM --enable-bam-indexing true --enable-duplicate-marking true &> $cramout/run_dragen.log";
		    
		    print STDERR "Running the dragen on $sample.\n\n";
		    
		    print STDERR "$cmd\n\n";
		    
		    if (!$obj->{_testing}){
			$rc = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "$cmd || echo 1"`;
			die "Failed running the dragen" if $rc =~ /1/;
			
			# move to mounted disk
			print STDERR "Moving files off dragen\n\n";
			$rc = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "mv $cramout $obj->{_drtmp}/ || echo 1"`;
			die "Cant move to mounted disk" if $rc =~ /1/;
			
		    }
		}
		
		print STDERR "Copying files to final directory\n\n";
		
		# copy files to new folder
		if (-e "$out/$sample" and $obj->{_force}){
		    `rm -Rf $out && echo 1` or die "Cant overwrite existing directory: $out";
		    
		} elsif (-e "$out/$sample" and !$obj->{_force}) {
		    die "Cant overwrite existing directory without force: $out";
		}
		
		if (!$obj->{_testing}){
		    `cp -r $obj->{_drtmponlocal}/$sample $obj->{_outdir}/ && echo 1` or die "Cant copy $sample to $obj->{_outdir}";
		    
		    # delete dragen copy
		    $rc = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "rm -Rf $obj->{_drtmp}/$sample || echo 1"`;
		    die "Cant remove $obj->{_drtmp}/$sample on dragen" if $rc =~ /1/;
		}
	    }
	    push @cramdirs, "$obj->{_outdir}/$sample";
	}
	$fastqfile =~ s!/mnt!/gscmnt!;
	`mv $fastqfile $obj->{_outdir}/`;    
	return @cramdirs;	
    }
    
    sub launch {
	my $obj = shift;
	my @dirs = @_;
	
	foreach my $i (0..$#dirs){
	    my $dir = $dirs[$i];
	    
	    die "$dir doesnt exist" if ! -d $dir;
	    
	    # get input json
	    my $inputs = from_json(`cat $obj->{_inputs}`);
	    
	    # get gender
	    if ($#{$obj->{_gender}} == $#dirs){
		$inputs->{'ChromoSeq.Gender'} = $obj->{_gender}[$i];
		
	    } elsif (`grep -c estimation $dir/*.ploidy_estimation_metrics.csv` =~ /1/) {
		`grep estimation $dir/*.ploidy_estimation_metrics.csv | cut -d ',' -f 4` =~ /(XX|XY|X0)/;
		if ($1 eq 'XX' or $1 eq 'X0'){
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
	    die "no counts file" if ! -e $counts;
	    $inputs->{'ChromoSeq.TumorCounts'} = $counts;
	    
	    my $name = basename($cram,".cram");
	    
	    $inputs->{'ChromoSeq.Name'} = $name;
	    
	    my $out = $obj->{_outdir} . '/' . $name;
	    mkdir $out if ! -e $out;
	    
	    $inputs->{'ChromoSeq.OutputDir'} = $out;
	    
	    open(JS,">$out/chromoseq.$name.json") || die;
	    print JS to_json($inputs);
	    close JS;
	    
	    my $cmd = "bsub -g $obj->{_group} -oo $out/out.$name.log -eo $out/err.$name.log -q research-hpc -a \"docker(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:7)\" /usr/bin/java -Dconfig.file=$obj->{_conf} -jar /opt/cromwell.jar run -t wdl -i $out/chromoseq.$name.json $obj->{_wdl}";
	    
	    print STDERR $cmd,"\n";
	    
	    `$cmd` if !$obj->{_testing};
	    
	}     
    }
    
    sub cleanup {
	my $obj = shift;
	
	# remove temp fastqs
	foreach my $d (@{$obj->{_toclean}}){
	    my $rc = `ssh gtac\@compute1-dragen-1.ris.wustl.edu "rm -Rf $d || echo 1"`;
	    die "Couldnt remove $d" if $rc =~ /1/;	
	}
    }
    
}
    
my $usage = "$0 <command> [options]\n\nCommands:\n\tprocess -d|dir <flowcell dir or fastq dir\n\tdemux -d|dir <flowcell directory> -o|out <outdir>\n\tmap -o|out <outdir> -d|dir <fastq dir> or -f <fastq samplesheet> [-n|name <sample>\n\tanalyze -c|cram <cram> -n|name <sample name> -o|out <outdir>";

my $debug = 0;
my $force = 0;
my $testing = 0;
my $out = '';
my $run = '';
my $group = '';
my $genders = '';
my $file = '';
my $inputs = ''; # for custom chromoseq run

my $cmd = shift @ARGV;

die "$usage" if $cmd !~ /flowcells|fastqs|alignments/;

GetOptions("debug" => \$debug,
           "o|out=s" => \$out,
	   "i|inputs=s" => \$inputs,
           "j|jobname=s" => \$group,
	   "g|genders=s" => \$genders,
	   "r|run=s" => \$run,
	   "f|file=s" => \$file,
           "t" => \$testing,
           "f" => \$force);

die $usage if (!$out or !-e $out or ($cmd ne 'alignments' and !$run));

$out = Cwd::realpath($out);
$out = $out . '/' . $run;

mkdir($out) if !-e $out;

my @in = @ARGV;

if ($file and -e $file){
    @in = `cat $file`;
    chomp @in;
}

my @genders = split(',',$genders);
if (scalar @genders > 0){
    map { die "gender must be male/female" if $_ !~/^male$|^female$/ } @genders;
}

die "number of genders must equal number of inputs" if scalar @genders > 0 and scalar @genders != scalar @in;

my $cs = CS->new(run=>$run,testing=>$testing,force=>$force,outdir=>$out,group=>$group,gender=>\@genders,inputs=>$inputs);

# do demux, if requested
if ($cmd eq 'flowcells'){

    my $fastqfile = $cs->demux(@in);
    my @cramdirs = $cs->align($fastqfile);
    $cs->launch(@cramdirs);
    $cs->cleanup;
}

# prepare and map fastqs that are on the dragen, or not.
elsif ($cmd eq 'fastqs'){

    my @samples = $cs->process_instrument_data(@in);

    foreach my $s (@samples){
	my $fastqfile = $cs->prepare_fastqs($s);
	my @cramdirs = $cs->align($fastqfile);
	$cs->cleanup;
	$cs->launch(@cramdirs);
    }
    
} elsif ($cmd eq 'alignments'){
    $cs->launch(@in);
}

print STDERR "done.\n\n";
