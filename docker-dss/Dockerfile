FROM ghcr.io/dhslab/docker-baseimage:231101
MAINTAINER David Spencer <dspencer@wustl.edu>

LABEL Base Spencer lab image with DSS added test5

RUN /usr/local/bin/Rscript -e 'install.packages("BiocManager"); BiocManager::install("DSS")'
