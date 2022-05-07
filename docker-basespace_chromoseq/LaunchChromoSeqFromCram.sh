#!/bin/bash

FILE=$1
OUTDIR=$2
NAME=$3

if [[ ! -e $OUTDIR/ ]];
then
    mkdir $OUTDIR/
fi

cat <<EOF > $OUTDIR/chromoseq.$NAME.json
{
    "ChromoSeq.CramIndex": "$FILE.crai",
    "ChromoSeq.Cram": "$FILE",
    "ChromoSeq.Name": "$NAME",
    "ChromoSeq.OutputDir": "$OUTDIR"
}
EOF

#LSF_PRESERVE_ENVIRONMENT=false bsub -g /dspencer/chromoseq -oo $OUTDIR/out.$NAME.log -eo $OUTDIR/err.$NAME.log -q research-hpc -a "docker(sleongmgi/cromwell:develop-with-mysql)" /usr/bin/java -Dconfig.file=/gscuser/dspencer/projects/wdltest/application.conf -jar /cromwell/cromwell.jar run /gscuser/dspencer/projects/chromoseq/Chromoseq.v7.hg38.wdl $OUTDIR/chromoseq.$NAME.json

LSF_PRESERVE_ENVIRONMENT=false bsub -g /dspencer/chromoseq -oo $OUTDIR/out.$NAME.log -eo $OUTDIR/err.$NAME.log -q research-hpc -a "docker(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:7)" /usr/bin/java -Dconfig.file=/gscmnt/gc2555/spencer/dhs/projects/wdltest/application.new.conf -jar /opt/cromwell.jar run -t wdl -i $OUTDIR/chromoseq.$NAME.json /gscuser/dspencer/projects/chromoseq/Chromoseq.v8.cromwell34.hg38.wdl
