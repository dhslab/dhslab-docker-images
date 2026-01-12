FROM python:3.11-slim
    
LABEL version="1.0"
LABEL description="Image for running SCGE off-target analysis with joblib for machine learning"

ENV DEBIAN_FRONTEND=noninteractive

# Install essential packages
#some basic tools
RUN apt-get update -y && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libcurl4-openssl-dev \
    unixodbc-dev \
    bzip2 \
    curl \
    g++ \
    less \
    libpng-dev \
    libssl-dev \
    libxml2-dev \
    make \
    ncurses-dev \
    pkg-config \
    unzip \
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
    && apt-get clean all \
    && rm -rf /var/lib/apt/lists/*

#set timezone to CDT
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime
RUN echo "America/Chicago" > /etc/timezone
RUN dpkg-reconfigure --frontend noninteractive tzdata

##############
#HTSlib 1.18#
##############
ENV HTSLIB_INSTALL_DIR=/opt/htslib

WORKDIR /tmp
RUN wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.21.tar.bz2 && \
    cd /tmp/htslib-1.21 && \
    ./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/ && \
    ln -s $HTSLIB_INSTALL_DIR/bin/tabix /usr/bin/tabix && \
    ln -s $HTSLIB_INSTALL_DIR/bin/bgzip /usr/bin/bgzip && \
    rm -Rf /tmp/htslib-1.21

################
#Samtools 1.18#
################
ENV SAMTOOLS_INSTALL_DIR=/opt/samtools

WORKDIR /tmp
RUN wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar --bzip2 -xf samtools-1.21.tar.bz2 && \
    cd /tmp/samtools-1.21 && \
    ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$SAMTOOLS_INSTALL_DIR && \
    make && \
    make install && \
    ln -s /opt/samtools/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/samtools-1.21

################
#bcftools 1.18#
################
ENV BCFTOOLS_INSTALL_DIR=/opt/bcftools
WORKDIR /tmp
RUN wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 && \
    tar --bzip2 -xf bcftools-1.21.tar.bz2 && \
    cd /tmp/bcftools-1.21 && \
    make prefix=$BCFTOOLS_INSTALL_DIR && \
    make prefix=$BCFTOOLS_INSTALL_DIR install && \
    ln -s /opt/bcftools/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/bcftools-1.21


# 3. Install Python Packages via PIP
# Note: We upgrade pip first to ensure binary wheels are found (speeding up build)
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
    numpy \
    joblib \
    scikit-learn \
    pandas \
    pyranges \
    openpyxl \
    scipy \
    cython \
    pysam \
    statsmodels \
    pyranges \
    biotite \
    pyodbc \
    ncls \
    edlib

