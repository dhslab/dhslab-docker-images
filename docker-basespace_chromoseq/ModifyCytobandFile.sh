#!/bin/bash

in=$1
out=$(basename $in .bed)".mod.bed"
			
awk 'BEGIN { last=""; } { if(/chrM|alt|HLA|Un|andom/){ print; } else if (NR==1){ print $1,$2,$3,"pter",$5; print; last=$0; lastchr=$1; } else if ($1 != lastchr) { gsub("q[0-9.]+","qter",last); print last; print $1,$2,$3,"pter",$5; print; lastchr=$1; last=$0; } else { print; lastchr=$1; last=$0; } }' $in > $out

bgzip $out
tabix -p bed $out.gz

