# docker-dss

FROM dhspence/docker-genomic-analysis
MAINTAINER David Spencer <dspencer@wustl.edu>

LABEL Image with DSS

RUN /usr/local/bin/Rscript -e 'install.packages("BiocManager"); BiocManager::install("DSS")'

