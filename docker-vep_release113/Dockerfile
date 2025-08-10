FROM ensemblorg/ensembl-vep:release_113.0

LABEL vep113 and htslib1.21 image

USER root
RUN mkdir -p /opt/lib/perl/VEP/Plugins
COPY Downstream.pm /opt/lib/perl/VEP/Plugins/Downstream.pm
COPY TranscriptCoordinates.pm /opt/lib/perl/VEP/Plugins/TranscriptCoordinates.pm
RUN apt-get update && apt-get install -y wget libncurses5-dev libncursesw5-dev libcurl4-openssl-dev

##############
#HTSlib 1.21#
##############
ENV HTSLIB_INSTALL_DIR=/opt/htslib

WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.21.tar.bz2 && \
    cd /tmp/htslib-1.21 && \
    ./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/ && \
    rm -Rf /tmp/htslib-1.21

################
#Samtools 1.21#
################
ENV SAMTOOLS_INSTALL_DIR=/opt/samtools

WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar --bzip2 -xf samtools-1.21.tar.bz2 && \
    cd /tmp/samtools-1.21 && \
    ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$SAMTOOLS_INSTALL_DIR && \
    make && \
    make install && \
    ln -s /opt/samtools/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/samtools-1.21

################
#bcftools 1.21#
################
ENV BCFTOOLS_INSTALL_DIR=/opt/bcftools
WORKDIR /tmp
RUN wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 && \
    tar --bzip2 -xf bcftools-1.21.tar.bz2 && \
    cd /tmp/bcftools-1.21 && \
    make prefix=$BCFTOOLS_INSTALL_DIR && \
    make prefix=$BCFTOOLS_INSTALL_DIR install && \
    ln -s /opt/bcftools/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/bcftools-1.21
