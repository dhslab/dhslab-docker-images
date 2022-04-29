FROM dhspence/docker-genomic-analysis2
MAINTAINER David H. Spencer <dspencer@wustl.edu>

LABEL \
  description="image for small rnaseq analysis"


COPY subread-1.6.2-Linux-x86_64.tar.gz /opt/
COPY bwa-0.7.17.tar.bz2 /opt/

WORKDIR /opt

RUN tar -zxvf subread-1.6.2-Linux-x86_64.tar.gz && \
    mv subread-1.6.2-Linux-x86_64/bin/featureCounts /usr/local/bin/

RUN tar -jxvf bwa-0.7.17.tar.bz2 && \
    cd bwa-0.7.17 && \
    make && \
    mv bwa /usr/local/bin/


RUN mkdir /opt/cutadapt

RUN pip install --install-option="--prefix=/opt/cutadapt" --upgrade cutadapt && \
    cp /opt/cutadapt/bin/cutadapt /usr/local/bin/cutadapt

RUN export PYTHONPATH=/opt/cutadapt/lib/python3.6/site-packages/