FROM ubuntu:20.04

LABEL version="1.0"
LABEL description="Image with basic bioinformatics tools and manta"

ENV DEBIAN_FRONTEND noninteractive

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
    groff \
    gnupg \
    unixodbc-dev \
    valgrind \
    && apt-get clean all

#set timezone to CDT
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime
RUN echo "America/Chicago" > /etc/timezone
RUN dpkg-reconfigure --frontend noninteractive tzdata

RUN apt-get update && apt-get install -y --no-install-recommends \
    python2 python2-dev wget ca-certificates && \
    ln -sf /usr/bin/python2.7 /usr/bin/python2 && \
    ln -sf /usr/bin/python2.7 /usr/bin/python && \
    update-alternatives --install /usr/bin/python python /usr/bin/python2.7 1 && \
    wget https://bootstrap.pypa.io/pip/2.7/get-pip.py && \
    python2 get-pip.py && \
    rm get-pip.py && \
    apt-get clean all

# Install htslib, samtools, and bcftools from source
ENV HTSLIB_VERSION 1.21
ENV SAMTOOLS_VERSION 1.21
ENV BCFTOOLS_VERSION 1.21
RUN cd /tmp && \
    wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/usr/local && make && make install && \
    cd /tmp && \
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/usr/local && make && make install && \
    cd /tmp && \
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    ./configure --prefix=/usr/local && make && make install && \
    rm -rf /tmp/htslib* /tmp/samtools* /tmp/bcftools*
    
ENV manta_version 1.6.0

WORKDIR /opt/

RUN wget https://github.com/Illumina/manta/releases/download/v${manta_version}/manta-${manta_version}.centos6_x86_64.tar.bz2 && \
    tar -jxvf manta-${manta_version}.centos6_x86_64.tar.bz2 && \
    mv manta-${manta_version}.centos6_x86_64 /usr/local/src/manta

COPY configManta.hg38.py.ini /usr/local/src/manta/bin/configManta.py.ini
