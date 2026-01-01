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

WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
    tar -xjvf htslib-1.21.tar.bz2 && \
    cd htslib-1.21 && \
    ./configure --prefix=/usr/local && make -j4 && make install && \
    cd .. && \
    \
    wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar -xjvf samtools-1.21.tar.bz2 && \
    cd samtools-1.21 && \
    ./configure --prefix=/usr/local --with-htslib=system && make -j4 && make install && \
    cd .. && \
    \
    wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 && \
    tar -xjvf bcftools-1.21.tar.bz2 && \
    cd bcftools-1.21 && \
    ./configure --prefix=/usr/local --with-htslib=system && make -j4 && make install && \
    cd .. && \
    rm -rf htslib-1.21* samtools-1.21* bcftools-1.21*

# 3. Install Python Packages via PIP
# Note: We upgrade pip first to ensure binary wheels are found (speeding up build)
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
    numpy \
    joblib \
    scikit-learn \
    pandas \
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

