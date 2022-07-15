FROM ubuntu:20.04
MAINTAINER Alex Miller <miller.alexander@wustl.edu>

LABEL Image for micro-c pipeline

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
    curl \
    git \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    bwa \
    bedtools \
    wget \ 
    libhdf5-dev \
    autotools-dev \
    automake

RUN pip install pysam && \
    pip install tabulate && \
    pip install numpy && \
    pip install scipy && \
    pip install py2bit && \
    pip install matplotlib && \
    pip install pyBigWig && \
    pip install deeptools && \
    pip install pandas && \
    pip install pairtools && \
    pip install cooler && \
    pip install h5py

#####################
# Micro-C get_qc.py #
#####################

WORKDIR /tmp
RUN git clone https://github.com/dovetail-genomics/Micro-C.git && \
cp Micro-C/get_qc.py /usr/bin && \
rm -r /tmp/Micro-C

###############
# HTSlib 1.13 #
###############

ENV HTSLIB_INSTALL_DIR=/opt/htslib

WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.13.tar.bz2 && \
    cd /tmp/htslib-1.13 && \
    ./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/ && \
    ln -s $HTSLIB_INSTALL_DIR/bin/tabix /usr/bin/tabix && \
    ln -s $HTSLIB_INSTALL_DIR/bin/bgzip /usr/bin/bgzip && \
    rm -Rf /tmp/htslib-1.13

#################
# Samtools 1.13 #
#################

ENV SAMTOOLS_INSTALL_DIR=/opt/samtools

WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && \
    tar --bzip2 -xf samtools-1.13.tar.bz2 && \
    cd /tmp/samtools-1.13 && \
    ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$SAMTOOLS_INSTALL_DIR && \
    make && \
    make install && \
    ln -s /opt/samtools/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/samtools-1.13

################
# Preseq 3.1.2 #
################

ENV PRESEQ_INSTALL_DIR=/opt/preseq

WORKDIR /tmp
RUN wget https://github.com/smithlabcode/preseq/releases/download/v3.1.2/preseq-3.1.2.tar.gz && \
    tar -xvf preseq-3.1.2.tar.gz && \
    cd preseq-3.1.2/ && \
    ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$PRESEQ_INSTALL_DIR && \
    make && \
    make install && \
    ln -s /opt/preseq/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/preseq-3*

##########
# PBGZIP #
##########

ENV PBGZIP_INSTALL_DIR=/opt/pbgzip

WORKDIR /tmp
RUN git clone https://github.com/nh13/pbgzip.git && \
    cd ./pbgzip/ && \
    sh autogen.sh && \
    ./configure --prefix=$PBGZIP_INSTALL_DIR && \
    make && \
    make install && \
    ln -s /opt/pbgzip/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/pbgzip*
