# Start from Ubuntu 22.04 LTS for a modern, stable base
FROM ubuntu:22.04

# Add labels for metadata
LABEL maintainer="David Spencer <dspencer@wustl.edu>"
LABEL description="Bioinformatics tools with Python 3.11 and Mambaforge"

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Chicago

# --- 1. Install Basic System Dependencies ---
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
    bzip2 \
    ca-certificates \
    curl \
    wget \
    gawk \
    tzdata \
    openssh-client \
    grep && \
    # Set timezone
    ln -sf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone && \
    dpkg-reconfigure --frontend noninteractive tzdata && \
    # Clean up to reduce image size
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# --- 2. Install Mambaforge and Create Conda Environment ---
ENV CONDA_DIR /opt/conda
WORKDIR /tmp
RUN mkdir ${CONDA_DIR} && \
        wget --quiet "https://github.com/conda-forge/miniforge/releases/download/25.3.1-0/Miniforge3-25.3.1-0-Linux-x86_64.sh" -O /tmp/mamba.sh

RUN /bin/bash /tmp/mamba.sh -f -b -p $CONDA_DIR && \
    rm /tmp/mamba.sh

# Add /opt/conda/bin to PATH
ENV PATH=/opt/conda/bin:$PATH

# htslib version
ARG HTSLIB_VERSION=1.21

# python version
ARG PYTHON_VERSION=3.11

# --- 3. Install All Software with a Single Mamba Command ---
RUN mamba install -y -c conda-forge -c bioconda \
    "python=${PYTHON_VERSION}" \
    # Command-line tools
    "htslib=${HTSLIB_VERSION}" \
    "samtools=${HTSLIB_VERSION}" \
    "bcftools=${HTSLIB_VERSION}" \
    bedtools \
    # Python libraries
    numpy \
    cython \
    pandas \
    pysam \
    cyvcf2 \
    pyranges \
    scipy \
    openpyxl \
    biopython && \
    # Clean up conda cache to save space
    mamba clean -y --all && \
    ln -s /opt/conda/bin/* /usr/local/bin/

