FROM dhspence/docker-baseimage:031822
MAINTAINER David Spencer <dspencer@wustl.edu>

LABEL Image with freebayes

RUN conda install -c bioconda freebayes
RUN export PATH=$PATH:/opt/conda/bin/
