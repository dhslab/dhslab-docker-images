FROM ubuntu:focal-20221130

LABEL version="2.1"
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
    ca-certificates \
    && apt-get clean all

#set timezone to CDT
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime
RUN echo "America/Chicago" > /etc/timezone
RUN dpkg-reconfigure --frontend noninteractive tzdata

WORKDIR /tmp/

RUN wget http://opengene.org/fastp/fastp && \
    chmod a+x ./fastp && \
    mv ./fastp /usr/local/bin/
