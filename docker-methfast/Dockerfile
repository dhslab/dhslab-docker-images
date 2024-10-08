FROM ubuntu:focal-20221130

LABEL version="1.0"
LABEL description="Image for methfast"

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
    unzip \
    wget \
    zip \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    bc \
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

# get dhslab/methfast from github
WORKDIR /tmp/
RUN wget --no-check-certificate https://github.com/dhslab/methfast/archive/refs/tags/methfast-0.1.tar.gz && \
    tar -xvf methfast-0.1.tar.gz && \
    cd methfast-methfast-0.1 && \
    make && make install && \
    rm -rf /tmp/methfast-0.1.tar.gz /tmp/methfast-methfast-0.1

    
