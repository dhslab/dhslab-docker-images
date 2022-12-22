FROM ubuntu:20.04
MAINTAINER "Dave Spencer" <dhspence@gmail.com>

LABEL description="merqury"	

ENV DEBIAN_FRONTEND=noninteractive

#some basic tools
RUN apt-get update && apt-get install -y --no-install-recommends locales && \
    echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    locale-gen en_US.UTF-8 && \
    LC_ALL=en_US.UTF-8 && \
    LANG=en_US.UTF-8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8 && \
    TERM=xterm
    
RUN apt-get update -y && apt-get install -y \
    wget \
    git \
    unzip \
    bzip2 \
    tar \
    build-essential \
    less \
    libcurl4-openssl-dev \
    libpng-dev \
    libssl-dev \
    libxml2-dev \
    make \
    ncurses-dev \
    nodejs \
    pkg-config \
    python3.9 \
    python3-pip \
    python3-dev \
    rsync \
    unzip \
    zip \
    libbz2-dev \
    bash-completion \
    ca-certificates \
    file \
    fonts-texgyre \
    gfortran \
    gsfonts \
    libbz2-1.0 \
    libjpeg-turbo8 \
    libopenblas-dev \
    libpangocairo-1.0-0 \
    libpcre3 \
    libtiff5 \
    liblzma5 \
    locales \
    zlib1g \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libpango1.0-dev \
    libjpeg-dev \
    libicu-dev \
    libpcre3-dev \
    libpng-dev \
    libreadline-dev \
    libtiff5-dev \
    liblzma-dev \
    libx11-dev \
    libxt-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    ncurses-dev \
    libcurl4-openssl-dev \
    bedtools \
    seqtk \
    cmake \
    autoconf

WORKDIR /tmp

ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

ENV PATH=$CONDA_DIR/bin:$PATH

RUN conda install -c bioconda graphaligner

##############
#HTSlib 1.13#
##############
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

################
#Samtools 1.13#
################
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


RUN cd /opt/ && \
    git clone https://github.com/lh3/minimap2 && \
    cd minimap2 && make && \
    cp minimap2 /usr/local/bin/


# svim-asm

WORKDIR /tmp

RUN pip install --upgrade pip && \
    git clone https://github.com/eldariont/svim-asm.git && \
    cd svim-asm && \
    pip install .

# gfatools

WORKDIR /tmp

RUN git clone https://github.com/lh3/gfatools && \
    cd gfatools && make && \
    mv gfatools /usr/local/bin/

# vg

WORKDIR /tmp

RUN wget https://github.com/vgteam/vg/releases/download/v1.43.0/vg && \
    chmod a+x vg && mv vg /usr/local/bin/


# minigraph

WORKDIR /tmp

RUN git clone https://github.com/lh3/minigraph && \
    cd minigraph && make && mv minigraph /usr/local/bin

# whatshap

RUN pip3 install whatshap

ENV PATH="/usr/local/bin:${PATH}"

RUN pip install statsmodels && \
    pip install modbamtools && \
    pip install Liftoff
