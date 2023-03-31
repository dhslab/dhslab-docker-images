FROM ubuntu:20.04
MAINTAINER David Spencer <dspencer@wustl.edu>

LABEL Image for micro-c/hi-c tools

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get install -y --no-install-recommends python3 \
    python3-dev \
    python3-pip \
    build-essential \
    bedtools \
    zlib1g-dev \
    libcurl4 \
    liblz4-tool \
    make \
    less \
    curl \
    libcurl4-gnutls-dev \
    git \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    bedtools \
    wget \ 
    libhdf5-dev \
    autotools-dev \
    automake \
    pypy \
    unzip

###############
# HTSlib 1.16 #
###############

ENV HTSLIB_INSTALL_DIR=/usr/local

WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.16.tar.bz2 && \
    cd /tmp/htslib-1.16 && \
    ./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/ && \
    rm -Rf /tmp/htslib-1.16

#################
# Samtools 1.16.1 #
#################

ENV SAMTOOLS_INSTALL_DIR=/usr/local

WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
    tar --bzip2 -xf samtools-1.16.1.tar.bz2 && \
    cd /tmp/samtools-1.16.1 && \
    ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$SAMTOOLS_INSTALL_DIR && \
    make && \
    make install && \
    cd / && \
    rm -rf /tmp/samtools-1.16.1

RUN pip install pysam tabulate numpy scipy py2bit matplotlib pyBigWig pandas pairtools cooler h5py PyGLtools

RUN git clone https://github.com/billgreenwald/pgltools.git && mv pgltools /usr/local/ && \
    ln -s /usr/local/pgltools/sh/pgltools /usr/local/bin/pgltools && \
    sed -i -e 's!$script_path!/usr/local/pgltools/sh!g' /usr/local/pgltools/sh/pgltools


