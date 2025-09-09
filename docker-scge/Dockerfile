FROM ubuntu:focal-20221130

LABEL version="1.0"
LABEL description="Image with basic bioinformatics tools plus joblib for machine learning"

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
    groff \
    mailutils \
    && apt-get clean all

#set timezone to CDT
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime
RUN echo "America/Chicago" > /etc/timezone
RUN dpkg-reconfigure --frontend noninteractive tzdata

# Install conda (mamba version)
ENV CONDA_DIR /opt/conda
RUN cd /tmp && mkdir -p $CONDA_DIR && \
    wget "https://github.com/conda-forge/miniforge/releases/download/24.3.0-0/Mambaforge-24.3.0-0-Linux-x86_64.sh" && \
    bash Mambaforge-24.3.0-0-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm -f Mambaforge-24.3.0-0-Linux-x86_64.sh && \
    $CONDA_DIR/bin/mamba install -y --channel conda-forge --channel bioconda python=3.11 numpy joblib scikit-learn bcftools==1.21 htslib==1.21 samtools==1.21 pandas bedtools=2.31 pybedtools openpyxl scipy cython cyvcf2 pysam statsmodels pyranges pybigwig globus-cli nf-core nextflow==24.04.4 SQLAlchemy biotite pyodbc pybigtools bigtools vcftools && \
    $CONDA_DIR/bin/conda clean -y --all && \
    ln -s /opt/conda/bin/* /usr/local/bin/

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && ./aws/install

##################
# ucsc utilities #
RUN mkdir -p /tmp/ucsc && \
    cd /tmp/ucsc && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig && \
    chmod ugo+x * && \
    mv * /usr/local/bin/
    
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
