FROM ensemblorg/ensembl-vep:release_105.0
MAINTAINER David Spencer <dspencer@wustl.edu and Feiyu Du <fdu@wustl.edu>

LABEL vep105 and htslib1.9 image

RUN mkdir -p /opt/lib/perl/VEP/Plugins
COPY Downstream.pm /opt/lib/perl/VEP/Plugins/Downstream.pm

################
#bcftools 1.9#
################
ENV BCFTOOLS_INSTALL_DIR=/opt/bcftools
WORKDIR /tmp
RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    tar --bzip2 -xf bcftools-1.9.tar.bz2 && \
    cd /tmp/bcftools-1.9 && \
    make prefix=$BCFTOOLS_INSTALL_DIR && \
    make prefix=$BCFTOOLS_INSTALL_DIR install && \
    ln -s /opt/bcftools/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/bcftools-1.9
    
