#!/usr/bin/perl

use File::Basename;

my ($reg,$R1,$R2) = @ARGV;

die "Cant find file $R1" if !-s $R1;
die "Cant find file $R2" if !-s $R2;
die "Region must be in format chr:start-end!" if $reg !~!-s $R2;

my $N = basename($R1,"_R1_001.fastq.gz");

# collapse reads
`/usr/local/bin/pear -f $R1 -r $R2 -o $N`;

# trim reads
`/usr/local/bin/trim_galore --path_to_cutadapt /opt/cutadapt-1.8.1/bin/cutadapt $N.assembled.fastq`;

# map reads
`minimap2 -a ~dspencer/refdata/GRCh37/all_sequences.fa $N.assembled_trimmed.fq | samtools view -b - | samtools sort - > $N.bam && samtools index $N.bam`;

# collect indels

my %lengths = ();
my %muts = ();

open(B,"samtools view $N.bam $reg |") || die "Cant run samtools to get indels";
 
while(<B>){
    chomp;
    my @l = split("\t",$_);

    my $pos = $l[3];
    my $len = 0;
    while($l[5] =~ /(\d+)[DIM]/g){ $len+=$1; }
    $lengths{$len}++;
    $muts{$pos . ":" . $l[5]}++;
}
close B;

my %sum = ();
my $modelen = (sort { $lengths{$b} <=> $lengths{$a} } keys %lengths)[0];
foreach my $mut (keys %muts){
    my $len = 0;
    while($mut =~ /(\d+)[DIM]/g){ $len+=$1; }
    next if $len < $modelen *.9;

    if ($mut !~ /[DI]/){
	$sum{wt}+=$muts{$mut};
    } else {
	$sum{$mut}+=$muts{$mut};
    }
}

open(O,">$N.indels.txt") || die "cant make indel output file";
my $total = 0;
map { $total+=$sum{$_} } keys %sum;
foreach my $k (sort {$sum{$b} <=> $sum{$a}} keys %sum){
    print O join("\t",$k,$sum{$k},sprintf("%.2f",$sum{$k}/$total*100)),"\n";
}
close O;

