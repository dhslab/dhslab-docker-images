# Use a standard Ubuntu 20.04 base image
FROM ubuntu:focal-20221130
MAINTAINER David Spencer <dspencer@wustl.edu>

LABEL cle_utilities="bioinformatics tools"

# Set non-interactive frontend for package installations
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Chicago

# Install basic dependencies and add the deadsnakes PPA for newer Python versions
RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    bzip2 \
    curl \
    wget \
    g++ \
    less \
    libcurl4-openssl-dev \
    libpng-dev \
    libssl-dev \
    libxml2-dev \
    libnss-sss \
    make \
    ncurses-dev \
    pkg-config \
    unzip \
    zip \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    bc \
    python3.11 \
    python3.11-dev \
    python3.11-distutils \
    python3-pip \
    gawk \
    openssh-client \
    bedtools && \
    # Clean up apt cache to reduce image size
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Update alternatives to make python3.11 the default python3
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1 && \
    update-alternatives --set python3 /usr/bin/python3.11

# Ensure pip is up-to-date
RUN python3 -m pip install --upgrade pip

############################################################
# BUILD HTSLIB, SAMTOOLS, AND BCFTOOLS FROM SOURCE (v1.21)
############################################################

ENV HTSLIB_INSTALL_DIR=/opt/htslib
ENV SAMTOOLS_INSTALL_DIR=/opt/samtools
ENV BCFTOOLS_INSTALL_DIR=/opt/bcftools
ENV HTS_VERSION=1.21

# HTSlib ${HTS_VERSION}
WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/${HTS_VERSION}/htslib-${HTS_VERSION}.tar.bz2 && \
    tar -xjf htslib-${HTS_VERSION}.tar.bz2 && \
    cd htslib-${HTS_VERSION} && \
    ./configure --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    rm -rf /tmp/htslib-${HTS_VERSION}*

# Samtools ${HTS_VERSION}
RUN wget https://github.com/samtools/samtools/releases/download/${HTS_VERSION}/samtools-${HTS_VERSION}.tar.bz2 && \
    tar -xjf samtools-${HTS_VERSION}.tar.bz2 && \
    cd samtools-${HTS_VERSION} && \
    ./configure --prefix=$SAMTOOLS_INSTALL_DIR --with-htslib=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    rm -rf /tmp/samtools-${HTS_VERSION}*

# bcftools ${HTS_VERSION}
RUN wget https://github.com/samtools/bcftools/releases/download/${HTS_VERSION}/bcftools-${HTS_VERSION}.tar.bz2 && \
    tar -xjf bcftools-${HTS_VERSION}.tar.bz2 && \
    cd bcftools-${HTS_VERSION} && \
    ./configure --prefix=$BCFTOOLS_INSTALL_DIR --with-htslib=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    rm -rf /tmp/bcftools-${HTS_VERSION}*

# Add the new tools to the system PATH and library path
ENV PATH=$HTSLIB_INSTALL_DIR/bin:$SAMTOOLS_INSTALL_DIR/bin:$BCFTOOLS_INSTALL_DIR/bin:$PATH
ENV LD_LIBRARY_PATH=$HTSLIB_INSTALL_DIR/lib:$LD_LIBRARY_PATH

############################################################
# INSTALL PYTHON PACKAGES
############################################################

RUN python3 -m pip install --no-cache-dir \
    numpy \
    cython \
    pandas \
    pysam \
    cyvcf2 \
    pyranges \
    scipy \
    openpyxl \
    biopython \
    natsort

# Set a default working directory for when the container starts
WORKDIR /data