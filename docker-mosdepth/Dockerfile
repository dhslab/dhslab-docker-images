FROM ghcr.io/dhslab/docker-htslib:latest

LABEL Basic image with htslib, samtools and mosdepth for bioinformatics

ENV DEBIAN_FRONTEND=noninteractive

################
# mosdepth 0.3.6
################

RUN wget --no-check-certificate https://github.com/brentp/mosdepth/releases/download/v0.3.6/mosdepth && \
    mv mosdepth /usr/local/bin/ && \
    chmod a+x /usr/local/bin/mosdepth
    
