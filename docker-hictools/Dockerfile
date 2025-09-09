FROM ubuntu:focal-20221130

LABEL Image for micro-c/hi-c tools

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get install -y --no-install-recommends ca-certificates \
    build-essential \
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
    unzip \
    tzdata \
    pkg-config \
    libssl-dev \
    libpng-dev

#set timezone to CDT
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime
RUN echo "America/Chicago" > /etc/timezone
RUN dpkg-reconfigure --frontend noninteractive tzdata

## Conda installation
# Configure environment
ENV CONDA_DIR /opt/conda
RUN cd /tmp && mkdir -p $CONDA_DIR && \
    wget "https://github.com/conda-forge/miniforge/releases/download/24.3.0-0/Mambaforge-24.3.0-0-Linux-x86_64.sh" && \
    bash Mambaforge-24.3.0-0-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm -f Mambaforge-24.3.0-0-Linux-x86_64.sh && \
    $CONDA_DIR/bin/mamba install -y --channel conda-forge --channel bioconda python=3.11 htslib==1.21 samtools==1.21 bedtools pysam tabulate numpy scipy py2bit matplotlib pandas h5py pyBigWig pairtools cooler cooltools hic2cool && \
    $CONDA_DIR/bin/conda clean -y --all && \
    ln -s /opt/conda/bin/* /usr/local/bin/
    


