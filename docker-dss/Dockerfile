FROM ghcr.io/dhslab/docker-baseimage:231101
MAINTAINER David Spencer <dspencer@wustl.edu>

LABEL Image with DSS

RUN /usr/local/bin/Rscript -e 'install.packages("BiocManager"); BiocManager::install("DSS")'
