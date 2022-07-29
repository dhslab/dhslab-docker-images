FROM ubuntu:groovy
MAINTAINER David Spencer <dspencer@wustl.edu>

LABEL Image for glimpse

#some basic tools

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Chicago

#set timezone to CDT
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime
#LSF: Java bug that need to change the /etc/timezone.
#     The above /etc/localtime is not enough.
RUN echo "America/Chicago" > /etc/timezone

RUN apt-get update -y && apt-get install -y --no-install-recommends \
    build-essential \
    locate \
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
    autoconf

RUN dpkg-reconfigure --frontend noninteractive tzdata

# needed for MGI data mounts
RUN apt-get update && apt-get install -y libnss-sss && apt-get clean all

# some other utils
RUN apt-get update && apt-get install -y --no-install-recommends gawk openssh-client grep evince && apt-get clean all

RUN apt-get install -y --no-install-recommends samtools bedtools bcftools bio-eagle shapeit4 libbz2-dev libssl-dev liblzma-dev libgsl0-dev

WORKDIR /opt

RUN git config --global http.sslverify false && git clone --branch=develop --recurse-submodules https://github.com/samtools/htslib.git && \
    git clone --branch=develop https://github.com/samtools/bcftools.git && \
    /bin/rm -f bcftools/plugins/{{mocha,beta_binom,genome_rules}.h,{mocha,trio-phase,mochatools,extendFMT}.c} && \
    git clone https://github.com/freeseek/mocha.git && \
    mv -t bcftools/plugins/ mocha/mocha.h mocha/beta_binom.h mocha/genome_rules.h mocha/mocha.c mocha/trio-phase.c mocha/mochatools.c mocha/extendFMT.c && \
    cd htslib && autoheader && (autoconf || autoconf) && ./configure --disable-bz2 --disable-gcs --disable-lzma && make && cd .. && \
    cd bcftools && make && cd .. && \
    /bin/cp bcftools/bcftools bcftools/plugins/*.so /usr/local/bin/ && \
    export BCFTOOLS_PLUGINS="/usr/local/bin"

RUN wget --no-check-certificate -O impute5_v1.1.4.zip https://www.dropbox.com/sh/mwnceyhir8yze2j/AAAPJVJv3kI2glXGDdc3sSHga/impute5_v1.1.4.zip?dl=0 && \
    unzip impute5_v1.1.4.zip && \
    /bin/mv impute5_v1.1.4/impute5_1.1.4_static /usr/local/bin/impute5 && \
    /bin/mv impute5_v1.1.4/imp5Converter_1.1.4_static /usr/local/bin/imp5Converter && \
    wget --no-check-certificate -P /usr/local/bin http://faculty.washington.edu/browning/beagle/bref3.18May20.d20.jar

RUN mkdir -p /usr/local/libexec/bcftools && /bin/cp bcftools/plugins/*.so /usr/local/libexec/bcftools/

