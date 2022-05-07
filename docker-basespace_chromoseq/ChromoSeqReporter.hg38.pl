#!/usr/bin/perl

use strict;
use JSON;

my %aa3to1 = qw(Ala A Arg R Asn N Asp D Asx B Cys C Glu E Gln Q Glx Z Gly G His H Ile I Leu L Lys K Met M Phe F Pro P Ser S Thr T Trp W Tyr Y Val V Xxx X Ter *);

sub lookup_bed {
    my ($c,$p,$h) = @_;
    my @ret = ();
    foreach my $k (keys %{$h}){
	if ($c eq $h->{$k}{chr} and $p > $h->{$k}{start} and $p < $h->{$k}{end}){
	    push @ret, $h->{$k};
	}
    }
    @ret;
}

sub lookup_translocation {
    my ($c1,$p1,$c2,$p2,$slop,$h) = @_;

    my @ret = ();
    my @oneside = ();
    foreach my $k (keys %{$h}){
	if ($c1 eq $h->{$k}{chr1} and $p1 > ($h->{$k}{start1} - $slop) and $p1 < ($h->{$k}{end1} + $slop)
	    and $c2 eq $h->{$k}{chr2} and $p2 > ($h->{$k}{start2} - $slop) and $p2 < ($h->{$k}{end2} + $slop)){
	    push @ret, $h->{$k};
	} elsif ($c1 eq $h->{$k}{chr1} and $p1 > ($h->{$k}{start1} - $slop) and $p1 < ($h->{$k}{end1} + $slop) and $h->{$k}{chr2} eq '0'){
	    push @oneside, $h->{$k};
	}	    
    }
    @ret = @oneside if scalar @ret == 0; # only allow the 'onesided' comparison if there are no twosided ones.
    @ret;      
}
    
my $CNVBED = "/opt/files/ChromoSeq.hg38.bed";
my $TRANSBED = "/opt/files/ChromoSeq.translocations.fixed.v3.sorted.hg38.bedpe";

my $slop = 0;
my $frac_for_whole_del = 0.75;
my $lowabundance = 10;
my $minSVlen = 100000;
my $maxSVlen = 5000000;
my $maxBlacklistfreq = 0.02;

my %genelist = ();

my ($name,$variants,$cnvs,$translocations) =  @ARGV;

print "ChromoSeq Report for $name ---- Generated on: " . localtime() . "\n\n";

my %arms = ();
my %chroms = ();
open(BED,$CNVBED) || die;
while(<BED>){
    chomp;
    my @l = split("\t",$_);
    $arms{$l[4]}++;
    $chroms{$l[0]}{$l[3]} = 1;
    
    $genelist{$l[5]}++ if $l[5] ne '.';
    
}
close BED;

print "Copy number alterations:\n\n";

my $anycnvs = 0;

open(F,$cnvs) || die "cant open file: $cnvs";

while(<F>){
    chomp;
    my ($chr,$start,$end,$segs,$l2r,$cn,$change,$bands,$band_count,$arms,$arm_count,$genes) = split("\t",$_);

    my $c = $chr;
    $c =~ s/chr//g;
    
    $genes =~ s/^\.,|\.//g;
    my $abund = '';
    
    if ($l2r > 0){
	$abund = ((2**$l2r - 1) / (($cn/2 - 1)))*100;
	# whole chromosome
	if ($band_count > (scalar(keys %{$chroms{$chr}}) * $frac_for_whole_del)){	    
	    print "seq[GRCh38] gain($c)\n\t[ " . "ploidy: $cn, " . "est. abundance: " . sprintf("%.1f\%",((2**$l2r - 1) / (($cn/2 - 1)))*100) . ", Genes affected: " . $genes . " ]";
	   	    
	} else {
	    my @bands = split(",",$bands);
	    print "seq[GRCh38] gain($c)($bands[0]$bands[$#bands])\nchr$c:g." . $start . "_" . $end . "gain\n\t[ " . "ploidy: $cn, " . "est. abundance: " . sprintf("%.1f\%",((2**$l2r - 1) / (($cn/2 - 1)))*100) . ", Genes affected: " . $genes . " ]";

	}
    } elsif ($l2r < 0){
	$cn = 1 if $cn > 2;
	$abund = ((2**$l2r - 1) / (($cn/2 - 1)))*100;
	# whole chromosome
	if ($band_count > (scalar(keys %{$chroms{$chr}}) * $frac_for_whole_del)){

#	    seq[GRCh38] type(#)(bkptband1bkptband2)
#            chr#:g.basepairbkpt1_basepairbkpt2type
# For example:  seq[GRCh38] del(X)(q21.31q22.1)
#	   chrX:g.89555676_100352080del
	    
	    print "seq[GRCh38] del($c)\n\t[ " . "ploidy: $cn, " . "est. abundance: " . sprintf("%.1f\%",((2**$l2r - 1) / (($cn/2 - 1)))*100) . ", Genes affected: " . $genes . " ]";

	} else {
	    my @bands = split(",",$bands);
	    print "seq[GRCh38] del($c)($bands[0]$bands[$#bands])\n" . "chr$c:g." . $start . "_" . $end . "del\n\t[ " . "ploidy: $cn, " . "est. abundance: " . sprintf("%.1f\%",((2**$l2r - 1) / (($cn/2 - 1)))*100) . ", Genes affected: " . $genes . " ]";
	    
	}
    }
    if ($abund < $lowabundance){
	print "\t*** LOW ABUNDANCE FINDING ***\n";
    } else {
	print "\n";
    }
    $anycnvs++;
}
close F;

if ($anycnvs == 0){
    print "***NO COPY NUMBER CHANGES IDENTIFIED***\n\n";
} else {
    print "\n\n";
}

print "Known Translocations:\n\n";

my %bands = ();
open(BED,$CNVBED) || die;
while(<BED>){
    chomp;
    my @l = split("\t",$_);
    $bands{"$l[0]:$l[1]-$l[2]"} = { chr => $l[0],
		      start => $l[1],
		      end => $l[2],
		      band => $l[3],
		      arm => $l[4] };
}
close BED;

my %trans = ();
open(BED,$TRANSBED) || die;
while(<BED>){
    chomp;
    my @l = split("\t",$_);
    $trans{$l[6]} = { chr1 => $l[0],
		      start1 => $l[1],
		      end1 => $l[2],
		      chr2 => $l[3],
		      start2 => $l[4],
		      end2 => $l[5],
		      band1 => $l[7],
		      band2 => $l[8],
		      genes => $l[6]};
    my ($g1,$g2) = split("_",$l[6]);

    $genelist{$g1}++;
    $genelist{$g2}++;
    
}
close BED;

my $anytrans = 0;

my %t2 = ();

my %types = ("BND" => "t",
	     "INV" => "inv",
	     "DEL"=> "del",
	     "DUP" => "dup");

my %isknown = ();
    
open(T,"gunzip -c $translocations |") || die;
while(<T>){
  next if /^#/;
  chomp;
  my @l = split("\t",$_);

  my $id = $l[2];
  
  my $foundtrans = 0;
  
  my ($chr1,$pos1,$chr2,$pos2,$svtype);
  
  if (/SVTYPE=BND/){
      $svtype = "BND";
      ($chr1,$pos1) = @l[0..1];
      $l[4] =~ /[\[\]](\S+):(\d+)[\[\]]/;
      $chr2 = $1;
      $pos2 = $2;
  } elsif (/SVTYPE=(INV|DEL|DUP|INS)/){
      $svtype = $1;
      ($chr1,$pos1) = @l[0..1];
      $chr2 = $chr1;
      $l[7] =~ /END=(\d+)/;
      $pos2 = $1;
  }

  my $filter = $l[6];

  my $len = 0;
  my $popfreq = 0.0;
  my $csblfreq = 0.0;
  
  if ($l[7] =~ /SVLEN=-*(\d+);/){
    $len = $1;
  }
  if ($l[7] =~ /POPFREQ_AF=(\S+?);/){
    $popfreq = $1;
  }
  
  if ($l[7] =~ /BLACKLIST_AF=(\S+?);/){
    $csblfreq = $1;
  }
  
  my $contig = 'None';
  if ($l[7] =~ /CONTIG=([ACTGNactgn]+);/){
      $contig = $1;
  }

  # get support
  my $paired_support = 0;
  my $paired_fraction = 0.0;
  my $split_support = 0;
  my $split_fraction = 0.0;  
  
  if ($l[9] =~ /(\d+),(\d+):(\d+),(\d+)$/){
      $paired_support = "$2/" . ($1+$2);
      $paired_fraction = ($2+$1) > 0 ? $2/($2+$1) : 0.0;
      $split_support = "$3/" . ($3+$4);
      $split_fraction = ($3+$4) > 0 ? $4/($3+$4) : 0.0;
  } elsif ($l[9] =~ /(\d+),(\d+)$/){
      $paired_support = "$2/" . ($1+$2);
      $paired_fraction = ($2+$1) > 0 ? $2/($2+$1) : 0.0;
  }  
  
  my @t = lookup_translocation($chr1,$pos1,$chr2,$pos2,$slop,\%trans);
  if (scalar @t > 0){
      
      my @p1 = lookup_bed($chr1,$pos1,\%bands);
      my @p2 = lookup_bed($chr2,$pos2,\%bands);
      
    foreach my $t (@t){
	my $c1 = $chr1;
	$c1 =~ s/chr//;
	my $c2 = $chr2;
	$c2 =~ s/chr//;
	my ($gene1,$gene2) = split("_",$t->{genes});
	print join("\t","seq[GRCh38] $types{$svtype}($c1;$c2)($p1[0]->{band};$p2[0]->{band})\n",
		   "$gene1--$gene2","$chr1:$pos1;$chr2:$pos2",
		   "PAIRED_READS: $paired_support (" . sprintf("%.1f\%",$paired_fraction*100) . ")",
		   "SPLIT_READS: $split_support (" . sprintf("%.1f\%",$split_fraction*100) . ")",
		   "Flags: " . $filter,
		   "Population frequency: " . ($popfreq * 100) . "%",
		   "ChromoSeq frequency: " . ($csblfreq  * 100) . "%",
		   "Contig: " . $contig),"\n\n";
	
	$isknown{$id} = 1;
    }
      $anytrans++;
  }


  if ($filter eq 'PASS' && $popfreq == 0 && $csblfreq < $maxBlacklistfreq && ($svtype eq 'BND' || ($svtype =~ /DEL|DUP|INS|INV/ && $len > $minSVlen))){   
    
    my $chr1 = $l[0];
    my $pos1 = $l[1];
    my $chr2 = '';
    my $pos2 = '';
    
    my $type = '';
    
    if ($l[7] =~ /SVTYPE=BND/){
      $type = "BND";
      $l[4] =~ /[\[\]](\S+):(\d+)[\[\]]/;
      $chr2 = $1;
      $pos2 = $2;
    } elsif ($l[7] =~ /SVTYPE=(DEL|DUP|INS|INV)/){
      $type = $1;
      $l[7]=~/END=(\d+)/;
      $chr2 = $chr1;
      $pos2 = $1;
    }
    
    my $pref = 0;
    my $palt = 0;
    my $sref = 0;
    my $salt = 0;
    if ($l[9] =~ /(\d+),(\d+):(\d+),(\d+)$/){
	$pref = $1;
	$palt = $2;
	$sref = $3;
	$salt = $4;
    }
    
    my @p1 = lookup_bed($chr1,$pos1,\%bands);
    $l[7] =~ /CSQ=(\S+)$/;
    my @csq = split(",",$1);
    my @fields = split /\|/, $csq[0];
    my $g = $fields[3];
    my $consequence = $fields[1];
    my $exon = $fields[8];
    my $intron = $fields[9];

    if ($type eq "BND"){
	$l[7] =~ /MATEID=([^;]+)/;
	my $mateid = $1;	
	$t2{$id} = [ $chr1, $pos1, $p1[0]->{band}, $type, $g, $consequence, $exon, $intron, $pref, $palt, $sref, $salt, $id, $mateid, $isknown{$id}, $contig, $csblfreq ];
	
    } elsif ($type =~ /DEL|DUP|INV|INS/){
	
	next if $type =~ /DEL|DUP/ and $pos2 - $pos1 > $maxSVlen; # skip if its a huge del/dup since we have CNV analysis for this
	
	my @p2 = lookup_bed($chr2,$pos2,\%bands);
	my %genes = ();
	foreach my $c (@csq){
	    my @f = split /\|/, $c;
	    $genes{$f[3]} = 1 if $f[7] eq "protein_coding";
	}
	
	if (scalar keys %genes > 10){
	    my %knowngenes = ();
	    map { $knowngenes{$_} = 1 if defined($genelist{$_}) } keys %genes;
	    $g = int(($pos2 - $pos1) / 1000) . " kbp, ". (scalar keys %genes) . " genes";
	    $g .= " (including: " . join(",",sort keys %knowngenes) . ")" if scalar keys %knowngenes > 0;
	} elsif (scalar keys %genes > 0) {
	    $g = int(($pos2 - $pos1) / 1000) . " kbp, genes: ". join(", ",sort keys %genes);
	} else {
	    $g = int(($pos2 - $pos1) / 1000) . " kbp";
	}
        $t2{$id} = [ $chr1, $pos1, $p1[0]->{band}, $type, $g, '', '', '', $pref, $palt, $sref, $salt, $id, $id, $isknown{$id}, $contig, $csblfreq ];
#	$t2{$id . "-END"} = [ $chr2, $pos2, $p2[0]->{band}, $type, $g, '', '', '', $pref, $palt, $sref, $salt, $id . "-END", $id, $isknown{$id}, $contig, $csblfreq ];
    }
    
  }
}

if ($anytrans == 0){
  print "***NO TRANSLOCATIONS IDENTIFIED***\n\n" 
} else {
  print "\n\n";
}

print "Gene-level hotspot analysis\n\n";

my @out = ();
open(F,$variants) || die;
<F>;
while(<F>){
    chomp;
    my @F = split("\t",$_);

#    splice(@F,6,4);

    next if $F[10]=~/synonymous|UTR|stream/ || $F[5] eq 'FilteredInAll';
    
    $F[14]=~s/\S+:(c\.\S+?)/\1/g; 
    $F[15]=~s/\S+:(p\.\S+?)/\1/g;
    my $HGVSpShort = $F[15];
    while( $HGVSpShort and my ( $find, $replace ) = each %aa3to1 ) {
        eval "\$HGVSpShort =~ s{$find}{$replace}g";
    }
 
    $F[10]=~/^(\S+?)_/; 
    my $var=uc($1); 
    $var .= " INSERTION" if length($F[4])>length($F[3]); 
    $var .= " DELETION" if length($F[4])<length($F[3]);
    $var .= " SITE VARIANT" if $var =~ /SPLICE/;
    #    $F[15]=~s/\/\d+//g;
    push @out, join("\t",$F[11],uc($var),$HGVSpShort,$F[14],$F[9],$F[0],$F[1],$F[3],$F[4]); 
}
close F;

if (scalar @out > 0){
    print join("\t",qw(GENE VARIANT HGVSp HGVSc VAF CHROM POS REF ALT)),"\n", join("\n",@out),"\n";
} else {
    print "***NO GENE-LEVEL VARIANTS IDENTIFIED***\n\n";
}


#
# unknown SVs that pass muster
#

my @list1 = ();
my @list2 = ();

map { 
  (($genelist{$t2{$_}[4]}>1)) ? 
    push @list1, $_ : push @list2, $_ } (sort { $t2{$a}[0] cmp $t2{$b}[0] } keys %t2);

if (scalar (@list1) > 0){
  print "\n\Previously unreported high-confidence structural variants that involve known genes\n\n";
  
} elsif (scalar (@list2) > 0) {
  print "\n\nPreviously unreported high-confidence structural variants\n\n";

} else {
  print "No other findings identified\n\n";
  exit;
}

# sort by SV type
@list2 = sort { $t2{$a}[3] cmp $t2{$b}[3] } @list2;

my @list = ();
if (scalar @list1 > 0){
  @list = (@list1,"space",@list2);
} elsif (scalar @list2 > 0){
  @list = @list2;
}

my $svt = '';

foreach my $v (@list){

  next if defined($t2{$v}[14]) and $t2{$v}[14] ne ''; # if this one was already reported as a hotspot event then skip

  next if !defined($t2{$t2{$v}[13]}); # if only one end passed all the filters then skip
  
  if ($v eq 'space' and scalar @list2 > 0){
    print "\n\nPreviously unreported high-confidence structural variants\n\n";
    next;
  } elsif ($svt ne '' and $t2{$v}[3] ne $svt){
      print "\n";
  }

  $svt = $t2{$v}[3];

  my @x = @{$t2{$v}};
  my @y = @{$t2{$t2{$v}[13]}};
  if (($x[0]=~/chr(\S+)/)[0] > ($y[0]=~/chr(\S+)/)[0]){
    my @tmp = @y;
    @y = @x;
    @x = @tmp;
  }
  my ($chr1,$pos1,$b1,$type1,$g1,$consequence1,$exon1,$intron1,$pref1,$palt1,$sref1,$salt1,$id11,$id12,$id13,$contig1,$freq1) = @x; #@{$t2{$v}};
  my ($chr2,$pos2,$b2,$type2,$g2,$consequence2,$exon2,$intron2,$pref2,$palt2,$sref2,$salt2,$id21,$id22,$id23,$contig2,$freq2) = @y; #@{$t2{$t2{$v}[13]}};

  if ($type1 eq 'BND'){
    $pref1 = $pref1 + $pref2;
    $palt1 = $palt1 + $palt2;
    $sref1 = $sref1 + $sref2;
    $salt1 = $salt1 + $salt2;
  }
  
  my $paired_support = $palt1 . "/" . ($palt1 + $pref1);
  my $paired_fraction = $palt1 + $pref1 > 0 ? $palt1 / ($palt1 + $pref1) : 0.0;
  my $split_support = $salt1 . "/" . ($salt1 + $sref1);
  my $split_fraction = $salt1 + $sref1 > 0 ? $salt1 / ($salt1 + $sref1) :  0.0;
  
  $exon1 = "exon $1" if ($exon1 =~ /(\d+)\/\d+/);
  $exon2 = "exon $1" if ($exon2 =~ /(\d+)\/\d+/);    
  $intron1 = "intron $1" if ($intron1 =~ /(\d+)\/\d+/);
  $intron2 = "intron $1" if ($intron2 =~ /(\d+)\/\d+/);
  
  if ($consequence1 =~ /intergenic/){
    $consequence1 = "INTERGENIC";
  } elsif ($consequence1 =~ /stream/) {
    $consequence1 = "$g1($consequence1)";
    } else {
      $consequence1 = "$g1($exon1$intron1)";
    }
  
  if ($consequence2 =~ /intergenic/){
    $consequence2 = "INTERGENIC";
  } elsif ($consequence2 =~ /stream/) {
    $consequence2 = "$g2($consequence2)";
  } else {
    $consequence2 = "$g2($exon2$intron2)";
  }
    
  if ($type1 eq 'BND'){

    print join("\t","seq[GRCh38] t(" . ($chr1=~/chr(\S+)/)[0] .";" . ($chr2=~/chr(\S+)/)[0] . ")($b1;$b2)\n",
	       "$consequence1--$consequence2","$chr1:$pos1;$chr2:$pos2",
	       "PAIRED_READS: $paired_support (" . sprintf("%.1f\%",$paired_fraction*100) . ")",
	       "SPLIT_READS: $split_support (" . sprintf("%.1f\%",$split_fraction*100) . ")",
	       "ChromoSeq frequency: " . ($freq1 * 100) . "%",
	       "Contig: " . $contig1 . "\n"),"\n";
  } else {
    
    print join("\t","seq[GRCh38] " . lc($type1) . "(" . ($chr1=~/chr(\S+)/)[0] . ")" . "(" . $b1 . $b2 . ")\n" .
	       sprintf("%s:g.%d_%d%s",$chr1,$pos1,$pos2,lc($type1)),"\n",$g1,
	       "PAIRED_READS: $paired_support (" . sprintf("%.1f\%",$paired_fraction*100) . ")",
	       "SPLIT_READS: $split_support (" . sprintf("%.1f\%",$split_fraction*100) . ")",
               "ChromoSeq frequency: " . ($freq1 * 100) . "%",
               "Contig: " . $contig1 . "\n"),"\n";
  }

}


#           seq[GRCh38] type(#)(bkptband1bkptband2)
#            chr#:g.basepairbkpt1_basepairbkpt2type
# For example:  seq[GRCh38] del(X)(q21.31q22.1)
#          chrX:g.89555676_100352080del
