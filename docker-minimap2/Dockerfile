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
    less \
    libcurl4-openssl-dev \
    libncurses5-dev \
    libpng-dev \
    libssl-dev \
    libxml2-dev \
    make \
    ncurses-dev \
    nodejs \
    pkg-config \
    unzip \
    wget \
    zip \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    bc \
    curl \
    tzdata \
    gawk \
    openssh-client \
    grep \
    evince \
    libnss-sss \
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
    $CONDA_DIR/bin/mamba install -y --channel conda-forge --channel bioconda htslib==1.21 samtools==1.21 bedtools=2.31 && \
    $CONDA_DIR/bin/conda clean -y --all && \
    ln -s /opt/conda/bin/* /usr/local/bin/

RUN cd /opt/ && \
    git clone https://github.com/lh3/minimap2 && \
    cd minimap2 && make && \
    cp minimap2 /usr/local/bin/
    
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

ADD pear-0.9.6-bin-64 /usr/local/bin/

RUN cp /usr/local/bin/pear-0.9.6-bin-64 pear
