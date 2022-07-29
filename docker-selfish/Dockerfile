FROM ubuntu:hirsute
MAINTAINER David Spencer <dhspence@gmail.com>

LABEL \
    description="Image for tools used in the CLE"

#LDAP and NSS
RUN echo "America/Chicago" > /etc/timezone
RUN apt-get update -y && apt-get install -y libnss-sss tzdata
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime

#LSF: Java bug that need to change the /etc/timezone.
#     The above /etc/localtime is not enough.
RUN dpkg-reconfigure --frontend noninteractive tzdata


RUN apt-get update -y && apt-get install -y \
    ant \
    apt-utils \
    git \
    libncurses5-dev \
    ncurses-dev \
    unzip \
    wget \
    zlib1g-dev

RUN apt-get update -y && apt-get install -y python3.6 python3-pip python3-dev build-essential nodejs libbz2-dev liblzma-dev
RUN pip3 install --upgrade pip
RUN pip3 install biopython
RUN pip3 install pysam
RUN pip3 install scipy
RUN pip3 install statsmodels
RUN ln -s $(which python3) /usr/local/bin/python

##############
#HTSlib 1.10.2#
##############
ENV HTSLIB_INSTALL_DIR=/opt/htslib

WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.10.2.tar.bz2 && \
    cd /tmp/htslib-1.10.2 && \
    ./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/ && \
    cp $HTSLIB_INSTALL_DIR/bin/tabix /usr/bin/ && \
    cp $HTSLIB_INSTALL_DIR/bin/bgzip /usr/bin/

################
#Samtools 1.10#
################
ENV SAMTOOLS_INSTALL_DIR=/opt/samtools

WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
    tar --bzip2 -xf samtools-1.10.tar.bz2 && \
    cd /tmp/samtools-1.10 && \
    ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$SAMTOOLS_INSTALL_DIR && \
    make && \
    make install && \
    cp /opt/samtools/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/samtools-1.10

WORKDIR /tmp
RUN git clone https://github.com/arq5x/bedtools2.git && \
    cd bedtools2 && \
    git checkout v2.27.0 && \
    make && \
    cp bin/* /usr/local/bin/ && \
    cd / && rm -rf /tmp/bedtools2
    

RUN pip3 install selfish-hic

RUN pip3 install serpentine

RUN pip3 install chess-hic

RUN apt-get update && apt-get install -y libnss-sss && apt-get clean all

