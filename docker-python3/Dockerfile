FROM ubuntu:20.04

LABEL David Spencer <dspencer@wustl.edu>

LABEL Base Spencerlab image with python3 and some standard bioinformatics modules

ENV DEBIAN_FRONTEND=noninteractive

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
    libdeflate-dev \
    make \
    ncurses-dev \
    pkg-config \
    python3 \
    python3-dev \
    virtualenv \
    python3-pip \
    zip \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    tzdata \
    gawk \
    openssh-client \
    grep \
    libnss-sss \
    && apt-get clean all && rm -rf /var/lib/apt/lists/*

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 1 && \
    update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1
    
RUN pip install numpy cython pyfaidx pandas pysam statsmodels pyranges gget openpyxl PyYAML

#set timezone to CDT
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime
RUN echo "America/Chicago" > /etc/timezone
RUN dpkg-reconfigure --frontend noninteractive tzdata
