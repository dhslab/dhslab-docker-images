FROM ubuntu:focal-20221130

LABEL maintainer="Andrew Bonney abonney@wustl.edu"
LABEL description="Docker image for HiCorr micro-C/Hi-C analysis"

ENV DEBIAN_FRONTEND=noninteractive

# Install essential packages
RUN apt-get update -y && apt-get install -y --no-install-recommends \
    build-essential \
    bzip2 \
    curl \
    g++ \
    git \
    less \
    libcurl4-openssl-dev \
    libpng-dev \
    libssl-dev \
    make \
    perl \
    wget \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && apt-get clean all \
    && rm -rf /var/lib/apt/lists/*

# Install conda (mamba version)
ENV CONDA_DIR /opt/conda
RUN cd /tmp && mkdir -p $CONDA_DIR && \
    wget "https://github.com/conda-forge/miniforge/releases/download/24.3.0-0/Mambaforge-24.3.0-0-Linux-x86_64.sh" && \
    bash Mambaforge-24.3.0-0-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm -f Mambaforge-24.3.0-0-Linux-x86_64.sh

# Install required packages via conda
RUN $CONDA_DIR/bin/mamba install -y --channel conda-forge --channel bioconda \
    bedtools \
    bowtie \
    samtools \
    numpy \
    pandas \
    scipy \
    perl-statistics-basic \
    perl-json \
    perl-yaml-tiny \
    && $CONDA_DIR/bin/conda clean -y --all

# Link conda binaries
RUN ln -s /opt/conda/bin/* /usr/local/bin/

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# Set locale
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
