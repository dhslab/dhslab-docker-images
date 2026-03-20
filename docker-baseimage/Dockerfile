FROM ubuntu:22.04

LABEL version="1.0"
LABEL description="Image with basic bioinformatics tools"

ENV DEBIAN_FRONTEND=noninteractive

# Install essential packages
# some basic tools
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
    nano \
    gawk \
    openssh-client \
    grep \
    evince \
    libnss-sss \
    openjdk-21-jdk \
    openjdk-21-jre \
    git-lfs \
    smbclient \
    groff \
    gnupg \
    unixodbc-dev \
    valgrind \
    && apt-get clean all


# 2. Add Microsoft's official repository
RUN curl https://packages.microsoft.com/keys/microsoft.asc | apt-key add - && \
    curl https://packages.microsoft.com/config/ubuntu/22.04/prod.list > /etc/apt/sources.list.d/mssql-release.list

# 3. Install the driver and clean up
RUN apt-get update && \
    ACCEPT_EULA=Y apt-get install -y msodbcsql18 && \
    rm -rf /var/lib/apt/lists/*

#set timezone to CDT
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime
RUN echo "America/Chicago" > /etc/timezone
RUN dpkg-reconfigure --frontend noninteractive tzdata

# Install Python 3.11 and pip
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.11 python3.11-venv python3.11-dev python3-pip && \
    ln -sf /usr/bin/python3.11 /usr/bin/python3 && \
    ln -sf /usr/bin/pip3 /usr/bin/pip && \
    update-alternatives --install /usr/bin/python python /usr/bin/python3 1 && \
    update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1 && \
    apt-get clean all

# Install the latest Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && chmod +x /usr/local/bin/nextflow

# Install Python packages via pip (remove bedtools, globus-cli, and nf-core, which are not pip packages)
RUN python3 -m pip install --upgrade pip && \
    pip install --no-cache-dir \
    nf-core pandas pybedtools openpyxl scipy cython cyvcf2 pysam statsmodels pyranges bamread pybigwig SQLAlchemy biotite pyodbc pybigtools bigtools sqlalchemy

# Install globus-cli (official method)
RUN apt-get update && apt-get install -y --no-install-recommends python3-pip && \
    pip install --upgrade pip && \
    pip install --no-cache-dir globus-cli

# Install bedtools from source
ENV BEDTOOLS_VERSION=2.31.1
RUN cd /tmp && \
    wget https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz && \
    tar -xzf bedtools-${BEDTOOLS_VERSION}.tar.gz && \
    cd bedtools2 && \
    make && \
    cp -r bin/* /usr/local/bin/ && \
    cd /tmp && rm -rf bedtools* 

# Install htslib, samtools, and bcftools from source
ENV HTSLIB_VERSION=1.21
ENV SAMTOOLS_VERSION=1.21
ENV BCFTOOLS_VERSION=1.21
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
