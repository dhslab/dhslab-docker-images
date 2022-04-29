#!/usr/bin/perl

use strict;

sub sum {
  my $s = 0;
  map {$s+=$_} @_;
  return $s;  
}

my ($hic,$NAME,$c,$res) = @ARGV;

#my @chroms=qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21);

#foreach my $c (@chroms){

  print STDERR "running $c\n";
 
  open(I,"/usr/bin/java -Xmx6g -jar /usr/local/bin/juicer_tools.jar dump observed NONE $hic $c $c BP $res |") or die;
  open(IO,"| gzip > $NAME.$c.int.txt.gz") or die;
  my %h = ();
  while(<I>){
    chomp;
    my @F = split("\t",$_);
    $h{$F[0]} += $F[2];
    $h{$F[1]} += $F[2] unless $F[0] == $F[1];
    print IO join("\t",$c,$F[0],$c,$F[1],$F[2]),"\n";
  }
  close I;
  close IO;

  open(OF,"| gzip -c > $NAME.$c.frag.txt.gz");
  foreach my $i (sort {$a<=>$b} keys %h){
    print OF join("\t",$c,"0",$i,$h{$i},"0"),"\n";
  }
  close OF;
  
  # run bias generation
#  `python HiCKRy.py -f $NAME.$c.frag.txt.gz -i $NAME.$c.int.txt.gz -o $NAME.$c.bias.txt` or die;
#  `fithic -f $NAME.$c.frag.txt.gz -i $NAME.$c.int.txt.gz -r $res -t $NAME.$c.bias.txt -U 10000000 -L 20000 -v -l $NAME.$c -o ./fithic` or die;

#}
