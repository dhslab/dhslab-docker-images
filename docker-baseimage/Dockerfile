FROM ubuntu:focal-20221130

LABEL version="1.0"
LABEL description="Image with basic bioinformatics tools"

ENV DEBIAN_FRONTEND=noninteractive

# Install essential packages
#some basic tools
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
    libxml2-dev \
    make \
    ncurses-dev \
    nodejs \
    pkg-config \
    rsync \
    unzip \
    wget \
    zip \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    bc \
    curl \
    tzdata \
    emacs \
    gawk \
    openssh-client \
    grep \
    evince \
    libnss-sss \
    default-jdk \
    default-jre \
    git-lfs \
    smbclient \
    && apt-get clean all

#set timezone to CDT
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime
RUN echo "America/Chicago" > /etc/timezone
RUN dpkg-reconfigure --frontend noninteractive tzdata

# Install conda (mamba version) with Jupyter
ENV CONDA_DIR /opt/conda
RUN cd /tmp && mkdir -p $CONDA_DIR && \
    wget "https://github.com/conda-forge/miniforge/releases/download/23.1.0-1/Mambaforge-23.1.0-1-Linux-x86_64.sh" && \
    bash Mambaforge-23.1.0-1-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm -f Mambaforge-23.1.0-1-Linux-x86_64.sh && \
    $CONDA_DIR/bin/mamba install -y --channel conda-forge --channel bioconda bcftools==1.18 htslib==1.18 bedtools=2.31 pandas==2.0.2 pybedtools==0.9.0 openpyxl==3.1.2 samtools==1.18 scipy==1.11.4 cython cyvcf2 pysam statsmodels pyranges globus-cli nextflow nf-core SQLAlchemy && \
    $CONDA_DIR/bin/conda clean -y --all && \
    ln -s /opt/conda/bin/* /usr/local/bin/

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
